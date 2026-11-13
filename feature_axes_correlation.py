#!/usr/bin/env python
"""
BURNER - asse cromatina/replicazione: correlazioni tra feature e soglia di raggruppamento del pilota
==================================================================================================
Domanda: tempo di replicazione e autovettori Hi-C portano la stessa informazione? Se si', il nucleo "solo replicazione"
in un run e "solo autovettori" in un altro e' lo stesso segnale attribuito diversamente (feature collineari).
Il pilota raggruppa le feature con average linkage su 1 - |Spearman|, soglia --group_corr (default 0.8).

Per ogni risoluzione (dataset in --datasets, etichette in --labels):
  - Spearman tra le feature dell'asse (--axis) e tabella completa sulle feature del pilota (--feature_list, opzionale:
    il feature_groups.tsv di un run del pilota; altrimenti tutte le colonne numeriche dell'asse);
  - con le stesse regole del pilota, il gruppo di ciascuna feature dell'asse a diverse soglie (--thresholds).
Figure (etichette in inglese): heatmap dell'asse per risoluzione; |rho| tra replicazione e ogni autovettore al variare
della risoluzione, con le soglie; con quali feature finisce la replicazione a ogni soglia.
Gli autovettori Hi-C vengono orientati come nel pilota (per cromosoma, correlazione con GC > 0) prima di tutto: il loro
segno e' arbitrario per cromosoma, e senza orientamento le correlazioni su tutto il genoma si annullano.
Solo lettura dei dataset: si puo' lanciare mentre altri job sono in corso.
"""
import argparse, os
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--datasets", nargs="+", required=True); p.add_argument("--labels", nargs="+", required=True)
p.add_argument("--feature_lists", nargs="+", default=None,
               help="un feature_groups.tsv del pilota per ogni dataset (stesse feature del pilota a quella risoluzione)")
p.add_argument("--axis", nargs="+", default=["replication.timing.smoothed", "GC_content"] + [f"eigen_{i}" for i in range(1, 7)])
p.add_argument("--thresholds", nargs="+", type=float, default=[0.8, 0.7, 0.6, 0.5])
p.add_argument("--rt", default="replication.timing.smoothed"); p.add_argument("--out_dir", default="feature_axes")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
assert len(a.datasets) == len(a.labels)
assert a.feature_lists is None or len(a.feature_lists) == len(a.datasets)
NICE = {"replication.timing.smoothed": "Replication timing", "GC_content": "GC content"}
nice = lambda f: NICE.get(f, f.replace("eigen_", "Hi-C eigen "))

R, rows, membership = {}, [], []
for k_, (ds, lab) in enumerate(zip(a.datasets, a.labels)):
    feats_all = list(pd.read_csv(a.feature_lists[k_], sep="\t")["feature"]) if a.feature_lists else None
    head = pd.read_csv(ds, sep="\t", nrows=0).columns
    feats = [f for f in (feats_all or a.axis) if f in head]
    axis = [f for f in a.axis if f in head]
    D = pd.read_csv(ds, sep="\t", usecols=sorted(set(feats) | set(axis) | {"binID", "GC_content"}))
    # stesso orientamento degli autovettori Hi-C del pilota (segno arbitrario per cromosoma): corr con GC > 0
    chrom = D["binID"].astype(str).str.split("_").str[0]; n_flip = 0
    eig = [c for c in D.columns if c.startswith("eigen_")]
    for ch, idx in D.groupby(chrom).groups.items():
        for c in eig:
            r = D.loc[idx, c].corr(D.loc[idx, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0: D.loc[idx, c] = -D.loc[idx, c]; n_flip += 1
    print(f"{lab}: autovettori orientati come nel pilota, {n_flip} segni invertiti su {len(eig) * chrom.nunique()}")
    rho = D[feats].rank().corr()                                          # Spearman = Pearson sui ranghi (come il pilota)
    rho.round(4).to_csv(f"{a.out_dir}/spearman_all_{lab}.tsv", sep="\t")
    R[lab] = rho.loc[axis, axis]
    for f in axis:
        if f != a.rt and a.rt in axis:
            rows.append(dict(resolution=lab, feature=f, abs_rho_with_RT=abs(rho.loc[a.rt, f])))
    A = rho.abs().fillna(0).to_numpy(copy=True); np.fill_diagonal(A, 1)
    Z = linkage(squareform(1 - A, checks=False), "average")
    for t in a.thresholds:
        g = pd.Series(fcluster(Z, t=1 - t, criterion="distance"), index=feats)
        if a.rt in g.index:
            mates = [f for f in g.index[g == g[a.rt]] if f != a.rt]
            membership.append(dict(resolution=lab, threshold=t, n_groups=g.nunique(), n_features=len(feats),
                                   RT_group=", ".join(mates) if mates else "(alone)"))
C = pd.DataFrame(rows); M = pd.DataFrame(membership)
C.to_csv(f"{a.out_dir}/abs_spearman_with_replication_timing.tsv", sep="\t", index=False)
M.to_csv(f"{a.out_dir}/replication_timing_group_by_threshold.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_colwidth", 120)
print("|Spearman| con il tempo di replicazione:\n", C.pivot(index="feature", columns="resolution", values="abs_rho_with_RT")[a.labels].round(2).to_string())
print("\nCon chi finisce il tempo di replicazione, per soglia (regole del pilota: average linkage, 1 - |rho|):\n", M.to_string(index=False))

# ---------------- figura 1: heatmap dell'asse per risoluzione (griglia, etichette solo sui bordi) ----------------
n = len(a.labels)
ncol = 2 if n > 1 else 1; nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol + 1.2, 4.9 * nrow), squeeze=False)
for k, lab in enumerate(a.labels):
    ax = axes[k // ncol, k % ncol]; r = R[lab]
    im = ax.imshow(r.to_numpy(), vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(r))); ax.set_yticks(range(len(r)))
    last_row = (k // ncol == nrow - 1) or (k + ncol >= n)
    ax.set_xticklabels([nice(f) for f in r.columns] if last_row else [], rotation=55, ha="right", fontsize=8)
    ax.set_yticklabels([nice(f) for f in r.index] if k % ncol == 0 else [], fontsize=8)
    for i in range(len(r)):
        for j in range(len(r)):
            v = r.iloc[i, j]; ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="white" if abs(v) > 0.6 else "black")
    ax.set_title(f"bin size {lab}", fontsize=10)
for k in range(n, nrow * ncol): axes[k // ncol, k % ncol].axis("off")
fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, label="Spearman correlation")
fig.suptitle("Chromatin / replication axis: correlation between features, by resolution\n"
             "(Hi-C eigenvectors oriented per chromosome: positive correlation with GC)", fontsize=11)
for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_axis_correlation_heatmaps.{ext}", dpi=200, bbox_inches="tight")

# ---------------- figura 2: |rho| con la replicazione al variare della risoluzione, e soglie ----------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.1, 1]})
x = np.arange(n)
for f, g in C.groupby("feature"):
    g = g.set_index("resolution").reindex(a.labels)
    ax1.plot(x, g["abs_rho_with_RT"], marker="o", lw=1.6 if f in ("eigen_1", "GC_content") else 1,
             label=nice(f), alpha=1 if f in ("eigen_1", "GC_content") else 0.6)
for t in a.thresholds:
    ax1.axhline(t, color="k", ls="--" if t == 0.8 else ":", lw=0.8)
    ax1.text(n - 1 + 0.08, t, f"{t:.1f}" + (" (pilot default)" if t == 0.8 else ""), va="center", fontsize=7)
ax1.set_xticks(x); ax1.set_xticklabels(a.labels); ax1.set_ylim(0, 1); ax1.set_xlim(-0.2, n - 0.2 + 0.9)
ax1.set_ylabel("|Spearman| with replication timing"); ax1.set_xlabel("bin size")
ax1.set_title("How strongly each feature tracks replication timing", fontsize=10)
ax1.legend(fontsize=7, ncol=4, loc="upper left", bbox_to_anchor=(0, -0.14), frameon=False)
# pannello 2: tabella di appartenenza
ax2.axis("off")
tab = M.pivot(index="threshold", columns="resolution", values="RT_group")[a.labels].sort_index(ascending=False)
def short(v, kmax=3):
    if v == "(alone)": return "alone"
    items = [nice(x) for x in str(v).split(", ")]
    return "\n".join(items) if len(items) <= kmax else f"{len(items)} features:\n" + "\n".join(items[:kmax]) + "\n…"
cell = [[short(v) for v in row] for row in tab.to_numpy()]
tb = ax2.table(cellText=cell, rowLabels=[f"|rho| ≥ {t:.1f}" for t in tab.index], colLabels=a.labels, loc="center", cellLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1, 2.6)
ax2.set_title("Features grouped with replication timing\n(pilot rule: average linkage on 1 − |Spearman|)", fontsize=10)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_replication_timing_grouping.{ext}", dpi=200, bbox_inches="tight")
print(f"\nOutput in {a.out_dir}/")