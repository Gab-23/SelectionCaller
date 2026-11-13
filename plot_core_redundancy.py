#!/usr/bin/env python
"""
BURNER - figure del core SEVERO (criterio used_strict di strict_usage.py) e della RIDONDANZA delle feature del core.
====================================================================================================================
Il pilota si allena su tutte le feature singole e raggruppa DOPO (|Spearman| >= 0.8, average linkage) solo per sommare
le importanze. Le feature correlate sotto quella soglia restano separate e si contendono l'importanza: questa figura
mostra, per ogni gruppo del core, quali altre feature portano informazione simile (la "ridondanza").

Figura 1  fig_strict_core_across_resolutions: righe = gruppi usati (criterio severo) da almeno una delle omiche scelte
          a qualche risoluzione; colonne = risoluzioni. Cerchio pieno = core (usato da TUTTE le omiche), mezzo cerchio =
          solo da una (lettera dell'omica); colore = importanza relativa media sulle omiche che lo usano.
Figura 2  fig_core_redundancy_<ris>: per ogni risoluzione, le feature dei gruppi ROBUSTI (core e moduli, con una striscia
          di colore per categoria) e le altre feature
          con |Spearman| >= --redundancy_rho con almeno un membro (vicine ridondanti): heatmap delle correlazioni
          (con segno) e, a destra, l'importanza di ogni feature per omica (|SHAP| relativo medio, modello vero).
          Autovettori Hi-C orientati come nel pilota (per cromosoma, correlazione con GC > 0).
Tabella   core_redundancy.tsv: per ogni risoluzione e gruppo del core, membri e vicine con il loro |rho|.

Uso: python plot_core_redundancy.py --omics atac mut_mis \\
       --run analyses/pilot_classes_floorimp/20000bp/genome_nopos 20kb DATASET_20kb ... --out_dir DIR
"""
import argparse, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.colors
from matplotlib.lines import Line2D

p = argparse.ArgumentParser()
p.add_argument("--run", nargs=3, action="append", metavar=("DIR", "LABEL", "DATASET"), required=True)
p.add_argument("--omics", nargs="+", default=["atac", "mut_mis"])
p.add_argument("--criterion", default="used_strict"); p.add_argument("--usage_file", default="usage_per_omic_strict.tsv")
p.add_argument("--redundancy_rho", type=float, default=0.5)
p.add_argument("--max_neighbours", type=int, default=8, help="vicine ridondanti mostrate nella heatmap (le piu' correlate); il totale e' nel titolo e nella tabella")
p.add_argument("--out_dir", default="core_redundancy")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
OMIC_NAME = {"cna_ampl": "CNA amplifications", "cna_del": "CNA deletions", "atac": "Chromatin accessibility",
             "mut_syn": "Synonymous mutations", "mut_mis": "Missense mutations", "mut_trunc": "Truncating mutations"}
OMIC_COL = {"atac": "#1b9e77", "mut_mis": "#d95f02", "cna_ampl": "#7570b3", "cna_del": "#e7298a", "mut_syn": "#66a61e", "mut_trunc": "#a6761d"}
NICE = {"replication.timing.smoothed": "Replication timing", "GC_content": "GC content"}
nice = lambda f: NICE.get(f, f.replace("eigen_", "Hi-C eigen "))
oms = a.omics

def load(rdir, dataset):
    U = pd.read_csv(f"{rdir}/{a.usage_file}", sep="\t"); U = U[U["omic"].isin(oms)]
    U = U[~U["group_label"].str.contains("telomer|centromer", case=False, regex=True)]
    use = U.pivot(index="group_label", columns="omic", values=a.criterion).reindex(columns=oms).fillna(False).astype(bool)
    imp = U.pivot(index="group_label", columns="omic", values="mean_imp").reindex(columns=oms).fillna(0)
    rel = imp / imp.sum(axis=0)
    groups = pd.read_csv(f"{rdir}/feature_groups.tsv", sep="\t")
    I = pd.read_csv(f"{rdir}/importance_per_chr.tsv", sep="\t")
    I = I[I["omic"].isin(oms) & ~I["feature"].str.startswith(("DECOY_", "SHADOW"))]
    fimp = I.groupby(["omic", "feature"])["imp"].mean().unstack(0).reindex(columns=oms).fillna(0)
    fimp = fimp / fimp.sum(axis=0)                                            # importanza relativa per feature, per omica
    feats = list(groups["feature"])
    D = pd.read_csv(dataset, sep="\t", usecols=lambda c: c in set(feats) | {"binID", "GC_content"})
    chrom = D["binID"].astype(str).str.split("_").str[0]
    for ch, idx in D.groupby(chrom).groups.items():                           # orientamento come nel pilota
        for c in [c for c in D.columns if c.startswith("eigen_")]:
            r = D.loc[idx, c].corr(D.loc[idx, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0: D.loc[idx, c] = -D.loc[idx, c]
    rho = D[[f for f in feats if f in D.columns]].rank().corr()
    return use, rel, groups, fimp, rho

data = {lab: load(d, ds) for d, lab, ds in a.run}
labels = [lab for _, lab, _ in a.run]

# ---------------- figura 1: core severo attraverso le risoluzioni ----------------
# NOTE (28/10/2026): i gruppi si calcolano a ogni risoluzione, quindi lo stesso gruppo puo' avere membri diversi
# ("3_PromD1 (+1)", "(+5)", "(+7)"). Una riga per feature RAPPRESENTATIVA (etichetta senza "(+n)"); il numero di
# altri membri a quella risoluzione e' scritto accanto al cerchio.
base_lab = lambda g: re.sub(r" \(\+\d+\)$", "", g)
n_extra = lambda g: int(re.search(r"\(\+(\d+)\)$", g).group(1)) if re.search(r"\(\+(\d+)\)$", g) else 0
rows = []
for lab in labels:
    use, rel, *_ = data[lab]
    for g in use.index[use.any(axis=1)]:
        n = int(use.loc[g].sum())
        rows.append(dict(group=base_lab(g), extra=n_extra(g), res=lab, n=n, who="+".join(o for o in oms if use.loc[g, o]),
                         imp=float(rel.loc[g, use.loc[g]].mean())))
R = pd.DataFrame(rows)
core_any = R[R["n"] == len(oms)]["group"].unique()
order = sorted(R["group"].unique(), key=lambda g: (-(R[(R.group == g) & (R.n == len(oms))].shape[0]), -R[R.group == g]["imp"].max()))
fig, ax = plt.subplots(figsize=(1.3 * len(labels) + 4.2, 0.32 * len(order) + 1.8))
vmax = R["imp"].max()
for i, g in enumerate(order):
    for j, lab in enumerate(labels):
        r = R[(R.group == g) & (R.res == lab)]
        if r.empty: continue
        r = r.iloc[0]; col = plt.cm.viridis(r["imp"] / vmax)
        if r["extra"] > 0:
            ax.text(j + 0.17, i, f"+{int(r['extra'])}", ha="left", va="center", fontsize=6.5, color="0.35", zorder=4)
        if r["n"] == len(oms):
            ax.scatter(j, i, s=260, color=col, edgecolor="k", lw=0.8, zorder=3)
        else:
            o = r["who"]; ax.scatter(j, i, s=150, color="white", edgecolor=OMIC_COL.get(o, "grey"), lw=2, zorder=3)
            ax.text(j, i, {"atac": "A", "mut_mis": "M"}.get(o, o[:1].upper()), ha="center", va="center", fontsize=6.5, color=OMIC_COL.get(o, "grey"), zorder=4)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels); ax.xaxis.tick_top()
ax.set_yticks(range(len(order))); ax.set_yticklabels([nice(g) for g in order], fontsize=8); ax.invert_yaxis()
ax.set_xlim(-0.6, len(labels) - 0.4); ax.grid(axis="x", color="0.9", zorder=0)
for k, g in enumerate(order):
    if g in core_any: ax.get_yticklabels()[k].set_fontweight("bold")
sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, vmax)); fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02,
    label="relative importance (mean |SHAP|, OOF)")
handles = [Line2D([], [], marker="o", ls="", ms=11, mfc="0.6", mec="k", label="core: robust in all omics")] + \
          [Line2D([], [], marker="o", ls="", ms=9, mfc="white", mec=OMIC_COL[o], mew=2, label=f"robust only in {OMIC_NAME[o].lower()}") for o in oms] + \
          [Line2D([], [], marker="$+n$", ls="", ms=12, color="0.35", label="n other features grouped with it at that resolution (|Spearman| ≥ 0.8)")]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.02 - 0.6 / (0.32 * len(order) + 1.8)), fontsize=7.5, frameon=False, ncol=1)
ax.set_title("Robust feature groups across resolutions (" + " + ".join(OMIC_NAME[o] for o in oms) + ")\n"
             "robust = above decoys AND above all importances with the shifted target, in most folds", fontsize=10, pad=28)
for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_strict_core_across_resolutions.{ext}", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---------------- figura 2: ridondanza delle feature del core, per risoluzione ----------------
tab = []
for lab in labels:
    use, rel, groups, fimp, rho = data[lab]
    # NOTE (28/10/2026): oltre al core, anche i MODULI (gruppi robusti per una sola omica), con una striscia di colore per
    # categoria; poi le vicine ridondanti NON robuste di tutti questi gruppi.
    robust = list(use.index[use.any(axis=1)])
    if not robust: print(f"{lab}: nessun gruppo robusto"); continue
    def cat_of(g):
        w = [o for o in oms if use.loc[g, o]]
        return "core" if len(w) == len(oms) else ("module:" + w[0] if len(w) == 1 else "shared:" + "+".join(w))
    cats = {g: cat_of(g) for g in robust}
    cat_order = ["core"] + [f"module:{o}" for o in oms] + sorted({c for c in cats.values() if c.startswith("shared:")})
    robust = sorted(robust, key=lambda g: (cat_order.index(cats[g]), -rel.loc[g, use.loc[g]].mean()))
    members, neigh, owner, fcat = [], [], {}, {}
    for g in robust:
        m = [f for f in groups.loc[groups["group_label"] == g, "feature"] if f in rho.index]
        members += m
        for f in m: owner[f] = g; fcat[f] = cats[g]
    for g in robust:
        m = [f for f in groups.loc[groups["group_label"] == g, "feature"] if f in rho.index]
        r = rho.loc[m].abs().max(axis=0).drop(index=members, errors="ignore")
        nb = r[r >= a.redundancy_rho].sort_values(ascending=False)
        for f, v in nb.items():
            if f not in neigh: neigh.append(f); owner.setdefault(f, f"~{g}")
        tab.append(dict(resolution=lab, group=g, category=cats[g], members=", ".join(m),
                        redundant_neighbours=", ".join(f"{f} ({v:.2f})" for f, v in nb.items())))
    n_tot = len(neigh)
    if n_tot > a.max_neighbours:                                  # le piu' correlate con le feature robuste
        best = rho.loc[members, neigh].abs().max(axis=0).sort_values(ascending=False)
        neigh = list(best.index[:a.max_neighbours])
    feats = members + neigh
    C = rho.loc[feats, feats]
    H = 0.38 * len(feats) + 2.4
    fig = plt.figure(figsize=(0.38 * len(feats) + 7.2, H))
    gs = fig.add_gridspec(1, 6, width_ratios=[0.25, 0.38 * len(feats), 0.25, 3.0, 0.6, 0.12], wspace=0.06)
    axs, ax, axb, cax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[3]), fig.add_subplot(gs[5])
    CAT_RGB = {"core": (0.15, 0.15, 0.15), **{f"module:{o}": matplotlib.colors.to_rgb(OMIC_COL.get(o, "grey")) for o in oms}}
    strip = [CAT_RGB.get(fcat[f], (0.45, 0.45, 0.45)) if f in fcat else (0.85, 0.85, 0.85) for f in feats]
    axs.imshow(np.array([[c] for c in strip]), aspect="auto")
    axs.set_xticks([]); axs.set_yticks(range(len(feats))); axs.set_yticklabels([nice(f) for f in feats], fontsize=7.5)
    cmap_ = plt.get_cmap("RdBu_r").copy(); cmap_.set_bad("0.85")
    im = ax.imshow(np.ma.masked_invalid(C.to_numpy(float)), vmin=-1, vmax=1, cmap=cmap_, aspect="auto")
    for i in range(len(feats)):
        for j in range(len(feats)):
            v = C.iat[i, j]
            if i != j and np.isfinite(v): ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6, color="white" if abs(v) > 0.6 else "black")
    ax.axhline(len(members) - 0.5, color="k", lw=1.2); ax.axvline(len(members) - 0.5, color="k", lw=1.2)
    ax.set_yticks([]); ax.set_xticks(range(len(feats))); ax.set_xticklabels([nice(f) for f in feats], rotation=60, ha="right", fontsize=7)
    y = np.arange(len(feats)); h = 0.8 / len(oms)
    for k, o in enumerate(oms):
        axb.barh(y - 0.4 + h * (k + 0.5), fimp.reindex(feats)[o].fillna(0).to_numpy(), height=h, color=OMIC_COL.get(o, "grey"), label=OMIC_NAME[o])
    axb.set_ylim(len(feats) - 0.5, -0.5); axb.set_yticks([]); axb.axhline(len(members) - 0.5, color="k", lw=1.2)
    axb.set_xlabel("relative importance of the single feature\n(mean |SHAP|, out-of-fold)", fontsize=8)
    axb.legend(fontsize=7, frameon=False, loc="lower right"); axb.spines[["top", "right"]].set_visible(False)
    cb = fig.colorbar(im, cax=cax); cb.set_label("Spearman correlation", fontsize=8); cb.ax.tick_params(labelsize=7)
    pos = cax.get_position(); cax.set_position([pos.x0 + 0.005, pos.y0 + pos.height * 0.25, pos.width, pos.height * 0.5])
    from matplotlib.patches import Patch
    leg = [Patch(color=CAT_RGB["core"], label="core (robust in all omics)")] + \
          [Patch(color=CAT_RGB[f"module:{o}"], label=f"module: robust only in {OMIC_NAME[o].lower()}") for o in oms] + \
          [Patch(color=(0.85, 0.85, 0.85), label=f"not robust, |Spearman| ≥ {a.redundancy_rho} with a robust feature")]
    fig.legend(handles=leg, loc="lower left", bbox_to_anchor=(0.01, 0.0), fontsize=7, frameon=False, ncol=2)
    fig.suptitle(f"Robust features and their redundancy (bin size {lab}): core and modules above the line, "
                 f"correlated non-robust features below ({len(neigh)} most correlated of {n_tot} shown)\nthe model is trained on single features: "
                 "correlated features compete, so their importance can move from one to another", fontsize=9.5, y=1.0 - 0.1 / H)
    fig.subplots_adjust(top=1 - 1.0 / H, bottom=1.6 / H)
    for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_core_redundancy_{lab}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
T = pd.DataFrame(tab); T.to_csv(f"{a.out_dir}/core_redundancy.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 160)
print(T.to_string(index=False)); print(f"\nOutput in {a.out_dir}/")