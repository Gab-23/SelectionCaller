#!/usr/bin/env python
"""
BURNER - figura "nucleo condiviso e moduli" (heatmap + barre laterali) per i run del pilota con le classi di mutazioni.
Stessa figura di summary_figures_pilot.py (fig1_core_vs_modules), aggiornata alle omiche attuali
(cna_ampl, cna_del, atac, mut_syn, mut_mis, mut_trunc). Legge solo usage_per_omic.tsv di ogni run.
Categorie: Core = gruppo usato (sopra le esche) da TUTTE le omiche del run; Module = da una sola; Partially shared =
da piu' di una ma non da tutte. Gruppi di distanze da telomero/centromero esclusi.
Colori = importanza relativa dentro ogni omica (|SHAP| medio OOF, ogni colonna somma a 1 sui gruppi mostrati).

Con --omics si scelgono le omiche: core / moduli sono calcolati SOLO su quelle (es. le mutazioni rappresentate dalle
sole missense, segnale osservato come per CNA e accessibilita').
RIDONDANZA (facoltativa, --datasets): una cella non robusta e' segnata ○ se quell'omica usa in modo robusto un altro gruppo
correlato (|Spearman| >= --redundancy_rho); un gruppo robusto in almeno un'omica e "coperto" cosi' in tutte le altre va nella
categoria "Core through correlated features" (stesso asse, rappresentato da feature diverse in omiche diverse).
Uso: python plot_core_modules.py --run analyses/pilot_classes/20000bp/genome_nopos 20kb \\
                                 --run analyses/pilot_classes/100000bp/genome_nopos 100kb ... --out_dir DIR
"""
import argparse, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.colors

p = argparse.ArgumentParser()
p.add_argument("--run", nargs=2, action="append", metavar=("DIR", "LABEL"), required=True)
p.add_argument("--omics", nargs="+", default=None,
               help="omiche da usare per le categorie e nella figura (default: tutte quelle del run), es. cna_ampl cna_del atac mut_mis")
p.add_argument("--usage_file", default="usage_per_omic.tsv",
               help="file di uso nella cartella del run (usage_per_omic_strict.tsv per il criterio used_strict)")
p.add_argument("--criterion", default="used", choices=["used", "used_robust", "used_shadow", "used_shadow_robust", "used_strict"],
               help="colonna di usage_per_omic.tsv che decide se un gruppo e' usato: used = esche (come prima); "
                    "used_shadow = copie spostate dei membri (--shadow_mode per_feature); *_robust = usato e NON usato nel pavimento "
                    "(--floor_importance); used_strict = usato rispetto alle esche E sopra il massimo del pavimento nella "
                    "maggioranza delle fold (strict_usage.py, --usage_file usage_per_omic_strict.tsv)")
p.add_argument("--datasets", nargs="+", default=None,
               help="(facoltativo) un dataset per run, nello stesso ordine dei --run: attiva la RIDONDANZA (vedi --redundancy_rho)")
p.add_argument("--redundancy_rho", type=float, default=0.5,
               help="con --datasets: una cella non robusta e' segnata ○ se l'omica usa in modo robusto un ALTRO gruppo con |Spearman| >= soglia "
                    "(massimo tra i membri dei due gruppi; autovettori orientati come nel pilota)")
p.add_argument("--out_dir", default="core_modules_figures")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
OMIC_NAME = {"cna_ampl": "CNA amplifications", "cna_del": "CNA deletions", "atac": "Chromatin accessibility",
             "mut_syn": "Synonymous mutations", "mut_mis": "Missense mutations", "mut_trunc": "Truncating mutations",
             "mut": "Coding mutations"}
OMIC_ORDER = ["cna_ampl", "cna_del", "atac", "mut_syn", "mut_mis", "mut_trunc", "mut"]
PALETTE = ["#1b9e77", "#7570b3", "#d95f02", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#1f78b4"]

def group_corr(rdir, dataset, groups_needed):
    """|Spearman| massimo tra i membri di due gruppi, per i gruppi richiesti (autovettori orientati: corr con GC > 0 per cromosoma)."""
    G = pd.read_csv(f"{rdir}/feature_groups.tsv", sep="\t")[["feature", "group_label"]]
    G = G[G["group_label"].isin(groups_needed)]
    feats = sorted(set(G["feature"]))
    D = pd.read_csv(dataset, sep="\t", usecols=lambda c: c in set(feats) | {"binID", "GC_content"}, dtype={"binID": str})
    chrom = D["binID"].str.split("_").str[0]
    for ch, idx in D.groupby(chrom).groups.items():
        for c in [c for c in D.columns if c.startswith("eigen_")]:
            r = D.loc[idx, c].corr(D.loc[idx, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0: D.loc[idx, c] = -D.loc[idx, c]
    R = D[[f for f in feats if f in D.columns]].rank().corr().abs()
    mem = {g: [f for f in G.loc[G["group_label"] == g, "feature"] if f in R.index] for g in groups_needed}
    out = pd.DataFrame(np.nan, index=groups_needed, columns=groups_needed)
    for g in groups_needed:
        for h in groups_needed:
            if mem[g] and mem[h]: out.loc[g, h] = np.nanmax(R.loc[mem[g], mem[h]].to_numpy())
    return out

def plot(rdir, label, dataset=None):
    U = pd.read_csv(f"{rdir}/{a.usage_file}", sep="\t")
    U = U[~U["group_label"].str.contains("telomer|centromer", case=False, regex=True)]
    oms = [o for o in OMIC_ORDER if o in set(U["omic"]) and (a.omics is None or o in a.omics)]
    U = U[U["omic"].isin(oms)]
    if a.criterion not in U.columns:
        raise SystemExit(f"{rdir}: colonna {a.criterion} assente in {a.usage_file} (rilanciare il pilota con le opzioni giuste, "
                         "o strict_usage.py per used_strict)")
    used = U.pivot(index="group_label", columns="omic", values=a.criterion).reindex(columns=oms).fillna(False).astype(bool)
    imp = U.pivot(index="group_label", columns="omic", values="mean_imp").reindex(columns=oms).fillna(0)
    # ---- ridondanza (facoltativa): una cella non robusta e' "coperta" se la stessa omica usa in modo robusto un altro gruppo
    #      correlato (|rho| >= soglia). Serve a distinguere "l'omica non usa questo asse" da "lo usa tramite una feature correlata".
    covered = pd.DataFrame(False, index=used.index, columns=oms); subst = {}
    if dataset is not None:
        shown = list(used.index[used.any(axis=1)])
        GC_ = group_corr(rdir, dataset, shown)
        for g in shown:
            for o in oms:
                if used.loc[g, o]: continue
                cands = [h for h in shown if h != g and used.loc[h, o] and np.isfinite(GC_.loc[g, h]) and GC_.loc[g, h] >= a.redundancy_rho]
                if cands:
                    best = max(cands, key=lambda h: GC_.loc[g, h]); covered.loc[g, o] = True
                    subst[(g, o)] = f"{best} ({GC_.loc[g, best]:.2f})"
    n = used.sum(axis=1); n_cov = (used | covered).sum(axis=1)
    cat = pd.Series(np.where(n == len(oms), "Core (all omics)",
                    np.where((n_cov == len(oms)) & (n >= 1), "Core through correlated features",
                    np.where(n == 1, "Module: " + used.idxmax(axis=1).map(lambda o: OMIC_NAME.get(o, o)),
                    np.where(n == 0, "Not used", "Partially shared")))), index=used.index)
    keep = cat != "Not used"; used, imp, cat, covered = used[keep], imp[keep], cat[keep], covered[keep]
    cat_order = ["Core (all omics)", "Core through correlated features", "Partially shared"] + [f"Module: {OMIC_NAME.get(o, o)}" for o in oms]
    order = sorted(used.index, key=lambda g: (cat_order.index(cat[g]) if cat[g] in cat_order else 99, -imp.loc[g].max()))
    used, imp, cat, covered = used.loc[order], imp.loc[order], cat.loc[order], covered.loc[order]
    M = imp / imp.sum(axis=0).replace(0, np.nan)
    cmap_cat = dict(zip(cat_order, ["#1b9e77", "#66c2a5"] + PALETTE[1:]))
    H = 0.3 * len(M) + 2.6
    fig = plt.figure(figsize=(1.15 * len(oms) + 6.0, H))
    fig.subplots_adjust(top=1 - 1.9 / H)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.25, 1.15 * len(oms), 2.4], wspace=0.08)
    axc, ax, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    axc.imshow(np.array([[matplotlib.colors.to_rgb(cmap_cat.get(c, "grey"))] for c in cat]), aspect="auto")
    axc.set_xticks([]); axc.set_yticks(range(len(M))); axc.set_yticklabels(M.index, fontsize=7)
    im = ax.imshow(M.to_numpy(float), aspect="auto", cmap="viridis")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if used.iat[i, j]: ax.text(j, i, "●", ha="center", va="center", color="white", fontsize=7)
            elif covered.iat[i, j]: ax.text(j, i, "○", ha="center", va="center", color="white", fontsize=8)
    ax.set_yticks([]); ax.set_xticks(range(len(oms))); ax.xaxis.tick_top()
    ax.set_xticklabels([OMIC_NAME.get(o, o) for o in oms], rotation=30, ha="left", fontsize=8)
    cb = fig.colorbar(im, ax=[axc, ax], orientation="horizontal", fraction=0.05, pad=0.04, aspect=40)
    cb.set_label("relative importance within each omic (mean |SHAP|, out-of-fold)", fontsize=8)
    counts = cat.value_counts().reindex([c for c in cat_order if c in set(cat)])
    axb.barh(range(len(counts)), counts.values, color=[cmap_cat.get(c, "grey") for c in counts.index])
    for i, v in enumerate(counts.values): axb.text(v, i, f" {v}", va="center", fontsize=7)
    axb.set_xlim(0, counts.max() * 1.25)
    axb.set_yticks(range(len(counts))); axb.set_yticklabels(counts.index, fontsize=8); axb.invert_yaxis()
    axb.yaxis.tick_right(); axb.set_xlabel("feature groups"); axb.spines[["top", "left"]].set_visible(False)
    fig.suptitle(f"Shared core and omic-specific modules of occurrence features (bin size {label})\n"
                 + {"used": "(● = above decoy features", "used_robust": "(● = above decoys and NOT above decoys with the shifted target",
                    "used_shadow": "(● = above shifted copies of its own features",
                    "used_shadow_robust": "(● = above shifted copies of its own features, and NOT with the shifted target",
                    "used_strict": "(● = above decoys, and above ALL its importances with the shifted target in most folds"}[a.criterion]
                 + "; telomere/centromere distances excluded)"
                 + (f"\n○ = not robust, but this omic robustly uses a correlated group (|Spearman| ≥ {a.redundancy_rho})" if dataset is not None else ""),
                 fontsize=10, y=1 - 0.25 / H)
    for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_core_vs_modules_{label}_{a.criterion}{'_redundancy' if dataset is not None else ''}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    T = pd.DataFrame({"category": cat}).join(used.add_prefix("used_")).join(M.add_prefix("relimp_"))
    if dataset is not None:
        for o in oms: T[f"correlated_substitute_{o}"] = [subst.get((g, o), "") for g in T.index]
    T.to_csv(f"{a.out_dir}/core_vs_modules_table_{label}_{a.criterion}{'_redundancy' if dataset is not None else ''}.tsv", sep="\t")
    return T, oms

if a.datasets is not None and len(a.datasets) != len(a.run): raise SystemExit("--datasets: uno per ogni --run, nello stesso ordine")
summary = []
for k, (rdir, label) in enumerate(a.run):
    T, oms = plot(rdir, label, a.datasets[k] if a.datasets else None)
    c = T["category"].value_counts()
    summary.append(dict(resolution=label, omics=len(oms), **{k: int(v) for k, v in c.items()},
                        core=", ".join(T.index[T["category"] == "Core (all omics)"]),
                        core_through_correlated=", ".join(T.index[T["category"] == "Core through correlated features"])))
S = pd.DataFrame(summary).fillna(0); S.to_csv(f"{a.out_dir}/core_vs_modules_summary_{a.criterion}{'_redundancy' if a.datasets else ''}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 200)
print(S.to_string(index=False)); print(f"\nOutput in {a.out_dir}/")