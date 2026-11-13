#!/usr/bin/env python
"""
BURNER - figure riassuntive del pilota (21/10/2026). Etichette delle figure in inglese.
=====================================================================================
Usa SOLO cartelle gia' prodotte: run del pilota (pilot_omics_background.py), tetti di rumore
(noise_ceiling.py), controllo dei driver (residual_vs_raw_drivers.py) e la matrice delle feature.

Figure
  1) fig1_core_vs_modules.(png|pdf)  gruppi di feature usati da tutte le omiche (core), da alcune
     (partially shared) o da una sola (module), dal run indicato con --modules_run (quello SENZA
     distanze da telomero/centromero). Eventuali gruppi di telomero/centromero rimasti sono esclusi.
  2) fig2_performance.(png|pdf)  R2 per omica e per run: raw e normalizzato per il tetto di rumore.
     Barra = R2 del modello genome-wide (run genome) o mediana sui cromosomi (run per cromosoma);
     punti = R2 per cromosoma; trattino nero = tetto. Le metriche sono RICALCOLATE dalle predizioni
     out-of-fold sulla stessa scala del tetto:
        CNA          R2 sulla frequenza di pazienti (bin di tutto il genoma)
        mutazioni    R2 sulla densita' (campioni mutati per kb di CDS), bin con cds_bp >= --min_cds,
                     SENZA bin driver (confronto con il tetto senza driver; la versione con i driver
                     e' nella tabella). Run glm (Poisson): conteggi / kb di CDS; run squared: expm1
                     del target log1p(densita').
        ATAC         R2 sul logFC; tetto non stimato -> solo raw
     Normalizzato = R2 / tetto (pavimento assunto 0; lo shift circolare non e' ancora stimato).
  3) fig3_mutations_drivers.(png|pdf)  mutazioni: AUROC con IC (osservato vs residuo), curve
     precision-recall nei bin con geni (con la prevalenza come linea del caso), arricchimento nel
     top 1% / 5%.
Tabelle: performance_summary.tsv, performance_per_chr.tsv.

Uso (dalla root della repo, env py_ml_env):
  python summary_figures_pilot.py --dataset $DS \
     --run pilot_chrom_full "Per-chromosome, with chr. coordinates" \
     --run pilot_genome_nopos "Genome-wide, no chr. coordinates" \
     --modules_run pilot_genome_nopos \
     --ceiling cna_ampl ceiling_ampl/noise_ceiling.tsv --ceiling cna_del ceiling_del/noise_ceiling.tsv \
     --ceiling mut ceiling_mut_nodrv/noise_ceiling.tsv \
     --drivers_tsv drivers_check_genome_nopos/residual_vs_raw.tsv --drivers_run pilot_genome_nopos \
     --out_dir summary_figures
"""
import argparse, os, sys, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--run", nargs=2, action="append", metavar=("DIR", "LABEL"), required=True)
p.add_argument("--modules_run", default=None, help="cartella del run senza coordinate per la figura 1")
p.add_argument("--ceiling", nargs=2, action="append", metavar=("OMIC", "TSV"), default=[])
p.add_argument("--drivers_tsv", default=None); p.add_argument("--drivers_run", default=None)
p.add_argument("--min_cds", type=float, default=500)
p.add_argument("--ymin_chr", type=float, default=-0.5, help="punti per cromosoma sotto questo valore disegnati al bordo")
p.add_argument("--out_dir", default="summary_figures")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)

OMIC_NAME = {"cna_ampl": "CNA amplifications", "cna_del": "CNA deletions",
             "mut": "Coding mutations", "atac": "Chromatin accessibility"}
OMIC_ORDER = ["cna_ampl", "cna_del", "mut", "atac"]
RUN_COLORS = ["#2c7fb8", "#d95f0e", "#31a354", "#756bb1", "#636363", "#e7298a"]
def save(fig, name):
    for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------- dataset: cds, driver, bin con geni ----------------
cols = pd.read_csv(a.dataset, sep="\t", nrows=0).columns
MUT_COLS = [c for c in ["mut_n_samples", "length_coding_agnostic"] if c in cols]
need = [c for c in ["binID", "cds_bp", "dist.to.closest.OG", "dist.to.closest.TSG",
                    "num.overlapping.genes_protein_coding"] + MUT_COLS if c in cols]
D = pd.read_csv(a.dataset, sep="\t", usecols=need, dtype={"binID": str})
D["is_driver"] = False
for c in ("dist.to.closest.OG", "dist.to.closest.TSG"):
    if c in D: D["is_driver"] |= D[c].fillna(np.inf) == 0
D["gene_bin"] = D["num.overlapping.genes_protein_coding"].fillna(0) > 0 if "num.overlapping.genes_protein_coding" in D else True
D = D.set_index("binID")

# ---------------- tetti ----------------
def load_ceilings():
    C = {}
    for om, f in a.ceiling:
        T = pd.read_csv(f, sep="\t")
        col = "split_half_SB" if "split_half_SB" in T else "alpha"
        if om == "mut":
            dens = T[T["target"].str.contains("densita|density", regex=True)]
            sel = {"nodrv": dens[dens["target"].str.contains("senza_driver|no_driver")],
                   "all": dens[~dens["target"].str.contains("senza_driver|no_driver")]}
        else:
            sel = {"all": T[~T["target"].str.contains("senza_driver|no_driver")]}
        for k, t in sel.items():
            if t.empty: continue
            C[(om, k)] = {str(s).replace("chr", ""): v for s, v in zip(t["scope"], t[col])}
    return C
CEIL = load_ceilings()
def ceil(om, key, scope):
    c = CEIL.get((om, key)); return np.nan if c is None else c.get(scope, np.nan)

# ---------------- R2 sulla scala del tetto, dalle predizioni OOF ----------------
def r2(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float); ok = np.isfinite(y) & np.isfinite(p)
    y, p = y[ok], p[ok]
    if len(y) < 20 or y.var() == 0: return np.nan
    return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)

def mut_scale(g, kb):
    """Scala del target mutazionale nelle predizioni OOF. I run piu' vecchi non hanno la colonna
    'family': si riconosce confrontando 'obs' con le trasformazioni candidate delle colonne del dataset."""
    if "family" in g and (g["family"] == "poisson").all(): return "poisson_count"
    obs = g["obs"].to_numpy(float)
    for c in MUT_COLS:
        cnt = g[c].to_numpy(float)
        for name, cand in [("log1p_density", np.log1p(cnt / kb)), ("log1p_count", np.log1p(cnt))]:
            ok = np.isfinite(cand) & np.isfinite(obs)
            if ok.sum() > 0 and np.nanmax(np.abs(cand[ok] - obs[ok])) < 1e-6:
                print(f"  mutazioni: target riconosciuto = {name} di {c}", file=sys.stderr)
                return name
    if "family" in g and (g["family"] == "gaussian").all():
        print("  mutazioni: scala non verificabile, assumo log1p(densita') (default del pilota)", file=sys.stderr)
        return "log1p_density"
    sys.exit("mutazioni: impossibile riconoscere la scala del target in oof_predictions.tsv "
             f"(colonne candidate nel dataset: {MUT_COLS}). Controlla la riga 'Mutazioni: target' nel log del run.")

def to_ceiling_scale(g, om, floor=False):
    """Ritorna dataframe con obs, pred sulla scala del tetto (+ chr, is_driver).
    floor=True: righe del pavimento (target spostato): esposizione e stato di driver vengono dal bin
    SORGENTE del valore (src_binID), la posizione (cromosoma) dal bin a cui il valore e' stato assegnato."""
    g = g.copy(); g["chr"] = g["binID"].str.split("_").str[0]
    key = "src_binID" if floor else "binID"
    if om == "mut":
        g = g.join(D[["cds_bp", "is_driver"] + MUT_COLS], on=key)
        g = g[g["cds_bp"] >= a.min_cds]
        kb = g["cds_bp"].to_numpy(float) / 1000.0
        if floor and "expo_kb" in g and np.isfinite(g["expo_kb"]).all(): kb = g["expo_kb"].to_numpy(float)
        scale = mut_scale(g, kb)
        if scale == "poisson_count":
            g["obs_s"], g["pred_s"] = g["obs"] / kb, g["pred"] / kb
        elif scale == "log1p_density":
            g["obs_s"], g["pred_s"] = np.expm1(g["obs"]), np.expm1(g["pred"])
        else:                                                  # log1p del conteggio -> densita'
            g["obs_s"], g["pred_s"] = np.expm1(g["obs"]) / kb, np.expm1(g["pred"]) / kb
    else:
        g["obs_s"], g["pred_s"] = g["obs"], g["pred"]
        g = g.join(D[["is_driver"]], on=key)
    return g

def floor_of_run(rdir, mode):
    """Pavimento sulla scala del tetto: {(omica, bins): (valore riassuntivo, {cromosoma: valore})}, media sulle ripetizioni."""
    f = f"{rdir}/floor_oof_predictions.tsv.gz"
    if not os.path.exists(f): return {}
    F = pd.read_csv(f, sep="\t", dtype={"binID": str, "src_binID": str}); out = {}
    for om, g0 in F.groupby("omic"):
        g0 = to_ceiling_scale(g0, om, floor=True)
        variants = [("nodrv", ~g0["is_driver"].fillna(False).astype(bool)), ("all", pd.Series(True, index=g0.index))] if om == "mut" \
                   else [("all", pd.Series(True, index=g0.index))]
        for key, m in variants:
            h = g0[m.to_numpy()]; summ, per = [], {}
            for rep, hr in h.groupby("rep"):
                pc = {c: r2(x["obs_s"], x["pred_s"]) for c, x in hr.groupby("chr")}
                for c, v in pc.items(): per.setdefault(c, []).append(v)
                summ.append(r2(hr["obs_s"], hr["pred_s"]) if mode == "genome" else np.nanmedian(list(pc.values())))
            out[(om, key)] = (float(np.nanmean(summ)), {c: float(np.nanmean(v)) for c, v in per.items()})
    return out

def norm(v, fl, ce):
    """Frazione della varianza spiegabile catturata OLTRE il pavimento: (R2 - pavimento) / (tetto - pavimento)."""
    fl = 0.0 if not np.isfinite(fl) else fl
    return (v - fl) / (ce - fl) if np.isfinite(ce) and np.isfinite(v) and ce > fl else np.nan

rows, rows_chr = [], []
for ri, (rdir, rlab) in enumerate(a.run):
    perf = pd.read_csv(f"{rdir}/performance.tsv", sep="\t")
    mode = "genome" if set(perf["unit"].astype(str)) == {"genome"} else "chrom"
    oof = pd.read_csv(f"{rdir}/oof_predictions.tsv", sep="\t", dtype={"binID": str})
    FLOOR = floor_of_run(rdir, mode)
    print(f"{rlab}: pavimento {'presente' if FLOOR else 'ASSENTE (assunto 0)'}", file=sys.stderr)
    for om, g in oof.groupby("omic"):
        g = to_ceiling_scale(g, om)
        variants = [("nodrv", ~g["is_driver"].fillna(False)), ("all", pd.Series(True, index=g.index))] if om == "mut" \
                   else [("all", pd.Series(True, index=g.index))]
        for key, m in variants:
            h = g[m.to_numpy()]
            per_chr = {c: r2(x["obs_s"], x["pred_s"]) for c, x in h.groupby("chr")}
            fl_sum, fl_chr = FLOOR.get((om, key), (np.nan, {}))
            for c, v in per_chr.items():
                cc = ceil(om, key, c); fc = fl_chr.get(c, np.nan)
                rows_chr.append(dict(run=rlab, run_dir=rdir, mode=mode, omic=om, bins=key, chr=c, R2=v, floor=fc,
                                     ceiling=cc, R2_norm=norm(v, fc, cc)))
            pc = pd.Series(per_chr, dtype=float)
            if mode == "genome":
                summ, cg = r2(h["obs_s"], h["pred_s"]), ceil(om, key, "genome")
            else:
                summ = float(pc.median()); cg = float(np.nanmedian([ceil(om, key, c) for c in pc.index])) if (om, key) in CEIL else np.nan
            rows.append(dict(run=rlab, run_dir=rdir, mode=mode, omic=om, bins=key, n_bins=len(h),
                             R2=summ, R2_per_chr_median=float(pc.median()), floor=fl_sum, ceiling=cg,
                             R2_norm=norm(summ, fl_sum, cg),
                             R2_norm_per_chr_median=float(np.nanmedian([norm(v, fl_chr.get(c, np.nan), ceil(om, key, c))
                                                                        for c, v in per_chr.items()]))
                             if (om, key) in CEIL else np.nan))
S = pd.DataFrame(rows); SC = pd.DataFrame(rows_chr)
S.to_csv(f"{a.out_dir}/performance_summary.tsv", sep="\t", index=False)
SC.to_csv(f"{a.out_dir}/performance_per_chr.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(S.round(3).to_string(index=False))

# ---------------- figura 2: performance raw e normalizzata ----------------
def plot_perf():
    P = S[(S["omic"] != "mut") | (S["bins"] == "nodrv")]
    PC = SC[(SC["omic"] != "mut") | (SC["bins"] == "nodrv")]
    oms = [o for o in OMIC_ORDER if o in set(P["omic"])]
    runs = [lab for _, lab in a.run]; nr = len(runs); w = 0.8 / nr
    fig, axes = plt.subplots(1, 2, figsize=(5 + 2.2 * len(oms), 4.6))
    for ax, (col, ccol, title) in zip(axes, [("R2", "R2", "Raw out-of-fold R²"),
                                             ("R2_norm", "R2_norm", "Normalized R²: (R² − floor) / (ceiling − floor)")]):
        for k, lab in enumerate(runs):
            color = RUN_COLORS[k % len(RUN_COLORS)]
            for i, om in enumerate(oms):
                x = i + (k - (nr - 1) / 2) * w
                r = P[(P["run"] == lab) & (P["omic"] == om)]
                if r.empty: continue
                v = r[col].iloc[0]
                if np.isfinite(v): ax.bar(x, v, w * 0.9, color=color, alpha=0.75, zorder=2)
                pts = PC[(PC["run"] == lab) & (PC["omic"] == om)][ccol].dropna().to_numpy()
                if len(pts):
                    jit = np.random.default_rng(k * 10 + i).uniform(-w * 0.3, w * 0.3, len(pts))
                    ax.scatter(x + jit, np.clip(pts, a.ymin_chr, None), s=7, color="k", alpha=0.45, zorder=3, lw=0)
                if col == "R2" and np.isfinite(r["ceiling"].iloc[0]):
                    ax.plot([x - w * 0.45, x + w * 0.45], [r["ceiling"].iloc[0]] * 2, color="k", lw=2, zorder=4)
                if col == "R2" and np.isfinite(r["floor"].iloc[0]):
                    ax.plot([x - w * 0.45, x + w * 0.45], [r["floor"].iloc[0]] * 2, color="k", lw=1.5, ls=":", zorder=4)
        ax.axhline(0, color="grey", lw=0.8)
        if col == "R2_norm":
            ax.axhline(1, color="k", lw=0.8, ls="--")
            for i, om in enumerate(oms):
                if not np.isfinite(P[P["omic"] == om]["R2_norm"].to_numpy(float)).any():
                    ax.text(i, 0.5, "no ceiling\nestimated", ha="center", va="center", fontsize=8, color="grey")
        ax.set_xticks(range(len(oms))); ax.set_xticklabels([OMIC_NAME.get(o, o) for o in oms], rotation=20, ha="right")
        ax.set_title(title, fontsize=10); ax.set_ylabel("R²" if col == "R2" else "fraction of explainable variance above floor")
        lo = min(a.ymin_chr, 0) if len(PC) else 0
        ax.set_ylim(lo - 0.05, 1.05 if col == "R2_norm" else max(1.0, np.nanmax(P["ceiling"].to_numpy(float), initial=0) + 0.05))
    handles = [Patch(color=RUN_COLORS[k % len(RUN_COLORS)], alpha=0.75, label=lab) for k, lab in enumerate(runs)]
    handles += [Line2D([], [], marker="o", ls="", color="k", alpha=0.5, ms=4, label="per chromosome"),
                Line2D([], [], color="k", lw=2, label="noise ceiling"),
                Line2D([], [], color="k", lw=1.5, ls=":", label="floor (circularly shifted target)")]
    fig.legend(handles=handles, loc="lower center", ncol=min(4, len(handles)), fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("Background models across omics (bars: genome-wide R² or median across chromosomes)", fontsize=11)
    fig.text(0.5, -0.17, f"Mutations: density per kb of CDS (bins with CDS ≥ {a.min_cds:g} bp), known driver bins excluded. "
             f"Chromatin accessibility: no ceiling estimated. Floor = same model fitted on the target circularly shifted along the genome "
             f"(0 if not computed). Per-chromosome points below {a.ymin_chr} drawn at the axis limit.",
             ha="center", fontsize=7, color="dimgrey")
    fig.tight_layout(); save(fig, "fig2_performance")

# ---------------- figura 1: core vs modules ----------------
def plot_modules():
    rdir = a.modules_run
    U = pd.read_csv(f"{rdir}/usage_per_omic.tsv", sep="\t")
    U = U[~U["group_label"].str.contains("telomer|centromer", case=False, regex=True)]
    oms = [o for o in OMIC_ORDER if o in set(U["omic"])]
    used = U.pivot(index="group_label", columns="omic", values="used").reindex(columns=oms).fillna(False).astype(bool)
    imp = U.pivot(index="group_label", columns="omic", values="mean_imp").reindex(columns=oms).fillna(0)
    n = used.sum(axis=1)
    cat = np.where(n == len(oms), "Core (all omics)",
          np.where(n == 1, "Module: " + used.idxmax(axis=1).map(lambda o: OMIC_NAME.get(o, o)),
          np.where(n == 0, "Not used", "Partially shared")))
    cat = pd.Series(cat, index=used.index)
    keep = cat != "Not used"
    used, imp, cat, n = used[keep], imp[keep], cat[keep], n[keep]
    cat_order = ["Core (all omics)", "Partially shared"] + [f"Module: {OMIC_NAME.get(o, o)}" for o in oms]
    order = sorted(used.index, key=lambda g: (cat_order.index(cat[g]) if cat[g] in cat_order else 99, -imp.loc[g].max()))
    used, imp, cat = used.loc[order], imp.loc[order], cat.loc[order]
    # importanza relativa per colonna (ogni omica somma a 1 sui gruppi mostrati) per confrontare le omiche
    M = imp / imp.sum(axis=0).replace(0, np.nan)
    cmap_cat = {c: col for c, col in zip(cat_order, ["#1b9e77", "#7570b3", "#d95f02", "#e7298a", "#66a61e", "#e6ab02"])}
    fig = plt.figure(figsize=(1.3 * len(oms) + 5.5, 0.3 * len(M) + 2.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.25, 1.3 * len(oms), 2.2], wspace=0.08)
    axc, ax, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    axc.imshow(np.array([[matplotlib.colors.to_rgb(cmap_cat.get(c, "grey"))] for c in cat]), aspect="auto")
    axc.set_xticks([]); axc.set_yticks(range(len(M))); axc.set_yticklabels(M.index, fontsize=7)
    im = ax.imshow(M.to_numpy(float), aspect="auto", cmap="viridis")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if used.iat[i, j]: ax.text(j, i, "●", ha="center", va="center", color="white", fontsize=7)
    ax.set_yticks([]); ax.set_xticks(range(len(oms))); ax.xaxis.tick_top()
    ax.set_xticklabels([OMIC_NAME.get(o, o) for o in oms], rotation=30, ha="left")
    cb = fig.colorbar(im, ax=[axc, ax], orientation="horizontal", fraction=0.05, pad=0.04, aspect=40)
    cb.set_label("relative importance within each omic (mean |SHAP|, out-of-fold)", fontsize=8)
    counts = cat.value_counts().reindex([c for c in cat_order if c in set(cat)])
    axb.barh(range(len(counts)), counts.values, color=[cmap_cat.get(c, "grey") for c in counts.index])
    axb.set_yticks(range(len(counts))); axb.set_yticklabels(counts.index, fontsize=8); axb.invert_yaxis()
    axb.yaxis.tick_right(); axb.set_xlabel("feature groups"); axb.spines[["top", "left"]].set_visible(False)
    fig.suptitle("Shared core and omic-specific modules of occurrence features\n"
                 "(● = above decoy features; telomere/centromere distances excluded)", fontsize=10, y=1.22)
    save(fig, "fig1_core_vs_modules")
    pd.DataFrame({"category": cat}).join(used.add_prefix("used_")).join(M.add_prefix("relimp_")).to_csv(
        f"{a.out_dir}/core_vs_modules_table.tsv", sep="\t")

# ---------------- figura 3: mutazioni e driver ----------------
def pr_curve(s, l):
    o = np.argsort(-s, kind="mergesort"); l = l[o]
    tp = np.cumsum(l); prec = tp / np.arange(1, len(l) + 1); rec = tp / max(l.sum(), 1)
    return rec, prec

def plot_mut_drivers():
    R = pd.read_csv(a.drivers_tsv, sep="\t"); R = R[R["omic"] == "mut"]
    name = {"osservato": "observed", "residuo": "residual", "observed": "observed", "residual": "residual"}
    R["score_en"] = R["score"].map(name)
    univ_name = {"all": "all bins", "gene_bins": "bins with genes"}
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    # (a) AUROC
    ax = axes[0]; unis = [u for u in ["all", "gene_bins"] if u in set(R["universe"])]
    for k, (sc, col) in enumerate([("observed", "#7f8c8d"), ("residual", "#c0392b")]):
        s = R[R["score_en"] == sc].set_index("universe").reindex(unis)
        x = np.arange(len(unis)) + (k - 0.5) * 0.38
        ax.bar(x, s["AUROC"], 0.38, color=col, label=sc,
               yerr=[s["AUROC"] - s["AUROC_lo"], s["AUROC_hi"] - s["AUROC"]], capsize=3)
    ax.axhline(0.5, color="k", lw=0.8, ls="--"); ax.set_ylim(0.3, 1); ax.set_xticks(range(len(unis)))
    ax.set_xticklabels([univ_name[u] for u in unis]); ax.set_ylabel("AUROC (known driver bins)")
    ax.set_title("a. Ranking of driver bins (AUROC)", fontsize=10); ax.legend(fontsize=8)
    # (b) precision-recall nei bin con geni
    ax = axes[1]
    if a.drivers_run:
        oof = pd.read_csv(f"{a.drivers_run}/oof_predictions.tsv", sep="\t", dtype={"binID": str})
        g = oof[oof["omic"] == "mut"].join(D[["is_driver", "gene_bin"]], on="binID")
        g = g[g["gene_bin"].fillna(False).astype(bool)]
        l = g["is_driver"].fillna(False).to_numpy(bool)
        for sc, col, colname in [("observed", "#7f8c8d", "obs"), ("residual", "#c0392b", "resid")]:
            s = g[colname].to_numpy(float); ok = np.isfinite(s)
            rec, prec = pr_curve(s[ok], l[ok])
            lo = l[ok][np.argsort(-s[ok], kind="mergesort")]
            ap = prec[lo].mean() if lo.any() else np.nan            # AP = media della precisione a ogni driver ritrovato
            ax.plot(rec, prec, color=col, lw=1.5, label=f"{sc} (AP = {ap:.3f})")
        ax.axhline(l.mean(), color="k", lw=0.8, ls="--", label=f"random (prevalence = {l.mean():.3f})")
        ax.set_xlabel("recall (fraction of driver bins found)"); ax.set_ylabel("precision (fraction of top bins that are drivers)")
        ax.set_ylim(0, 1); ax.legend(fontsize=7)
    ax.set_title("b. Precision-recall, bins with genes", fontsize=10)
    # (c) arricchimento nel top 1% / 5%
    ax = axes[2]; s = R[R["universe"] == "gene_bins"] if "gene_bins" in set(R["universe"]) else R
    for k, (sc, col) in enumerate([("observed", "#7f8c8d"), ("residual", "#c0392b")]):
        t = s[s["score_en"] == sc]
        if t.empty: continue
        ax.bar(np.arange(2) + (k - 0.5) * 0.38, [t["enrich_top1"].iloc[0], t["enrich_top5"].iloc[0]], 0.38, color=col, label=sc)
    ax.axhline(1, color="k", lw=0.8, ls="--"); ax.set_xticks(range(2)); ax.set_xticklabels(["top 1%", "top 5%"])
    ax.set_ylabel("enrichment of driver bins (× prevalence)"); ax.set_title("c. Enrichment at the top of the ranking", fontsize=10)
    ax.legend(fontsize=8)
    fig.suptitle("Coding mutations: does the residual of the occurrence model point to known drivers?", fontsize=11)
    fig.tight_layout(); save(fig, "fig3_mutations_drivers")

plot_perf()
if a.modules_run: plot_modules()
if a.drivers_tsv: plot_mut_drivers()
print(f"\nFigure e tabelle in {a.out_dir}/")