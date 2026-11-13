#!/usr/bin/env python
"""
BURNER - confronto per gene: il nostro residuo (gene_level_drivers_classes.py) contro dNdScv sugli stessi dati.
Stessa valutazione per tutti: AP (pari merito come un unico gradino), AUROC, driver nei primi 20/100, IC con bootstrap
sui geni; stesso universo (geni presenti in entrambi) e stessa lista di driver (colonna is_driver dei nostri file).
Corrispondenze: missense -> pmis_cv; troncanti -> ptrunc_cv; non sinonime -> pallsubs_cv (tutte le sostituzioni).
Punteggi: nostro = residuo (-log10 coda binomiale negativa), osservato = conteggio; dNdScv = -log10 p.
Uso: python compare_with_dndscv.py --ours analyses/pilot_classes/100000bp/drivers_gene_level --label 100kb \
       --dndscv analyses/dndscv_BRCA/sel_cv_hg19.tsv [--dndscv_label dNdScv] --out_dir DIR
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
p = argparse.ArgumentParser()
p.add_argument("--ours", nargs="+", required=True); p.add_argument("--label", nargs="+", required=True)
p.add_argument("--dndscv", nargs="+", required=True); p.add_argument("--dndscv_label", nargs="+", default=None)
p.add_argument("--ours_variants", nargs="*", default=[],
               help="altri attesi nostri da gene_scores_variants_<classe>.tsv (es. sites length length+syn model+syn)")
p.add_argument("--n_boot", type=int, default=1000); p.add_argument("--seed", type=int, default=1)
p.add_argument("--out_dir", default="compare_dndscv")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True); rng = np.random.default_rng(a.seed)
dl = a.dndscv_label or [os.path.basename(f).replace(".tsv", "") for f in a.dndscv]
MAP = {"mis": "pmis_cv", "trunc": "ptrunc_cv", "nonsyn": "pallsubs_cv"}
def pr_steps(s, l):
    o = np.argsort(-s, kind="mergesort"); ss, lo = s[o], l[o]
    last = np.r_[np.where(np.diff(ss) != 0)[0], len(ss) - 1]; tp = np.cumsum(lo)[last]
    return tp / max(l.sum(), 1), tp / (last + 1)
def ap(s, l):
    if not l.any(): return np.nan
    rec, prec = pr_steps(s, l); return float(np.sum(np.diff(np.r_[0, rec]) * prec))
def auroc(s, l):
    r = pd.Series(s).rank().to_numpy(); n1 = l.sum(); n0 = len(l) - n1
    return (r[l].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
def topk(s, l, k):
    o = np.lexsort((rng.random(len(s)), -s)); return int(l[o[:k]].sum())
rows = []
D = {lab: pd.read_csv(f, sep="\t").drop_duplicates("gene_name").set_index("gene_name") for f, lab in zip(a.dndscv, dl)}
for odir, olab in zip(a.ours, a.label):
    for cls, pcol in MAP.items():
        G = pd.read_csv(f"{odir}/gene_scores_{cls}.tsv", sep="\t", index_col=0)
        common = G.index
        for lab, Dd in D.items(): common = common.intersection(Dd.index)
        G = G.loc[common]; l = G["is_driver"].to_numpy(bool)
        bottom = lambda s: np.where(np.isfinite(s), s, np.nanmin(s[np.isfinite(s)]) - 1)
        scores = {"observed (count)": G["O"].to_numpy(float), "ours (residual)": bottom(G["residual"].to_numpy(float))}
        if a.ours_variants:
            V = pd.read_csv(f"{odir}/gene_scores_variants_{cls}.tsv", sep="\t", index_col=0).reindex(common)
            for v in a.ours_variants:
                if v not in V.columns: raise SystemExit(f"{odir}: variante {v} assente (disponibili: {list(V.columns)})")
                scores[f"ours: {v}"] = bottom(V[v].to_numpy(float))
        for lab, Dd in D.items():
            pv = Dd.loc[common, pcol].to_numpy(float); scores[lab] = bottom(-np.log10(np.clip(pv, 1e-300, 1)))
        pos, neg = np.where(l)[0], np.where(~l)[0]
        for name, s in scores.items():
            bs = []
            for _ in range(a.n_boot):
                i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]; bs.append(ap(s[i], l[i]))
            rows.append(dict(resolution=olab, cls=cls, method=name, n_genes=len(s), n_drivers=int(l.sum()), AP=ap(s, l),
                             AP_lo=np.percentile(bs, 2.5), AP_hi=np.percentile(bs, 97.5), AUROC=auroc(s, l),
                             drivers_top20=topk(s, l, 20), drivers_top100=topk(s, l, 100)))
R = pd.DataFrame(rows); R.to_csv(f"{a.out_dir}/comparison_with_dndscv.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(R.round(3).to_string(index=False))
# figura: AP con IC per classe e metodo, una colonna per risoluzione
NAME = {"mis": "Missense", "trunc": "Truncating", "nonsyn": "Non-synonymous"}
meths = list(dict.fromkeys(R["method"]))
ours_v = [m for m in meths if m.startswith("ours: ")]; dn = [m for m in meths if m not in ("observed (count)", "ours (residual)") and m not in ours_v]
meths = ["observed (count)", "ours (residual)"] + ours_v + dn
cols = ["#7f8c8d", "#c0392b"] + list(plt.cm.Oranges(np.linspace(0.35, 0.75, max(len(ours_v), 1))))[:len(ours_v)] + ["#2c3e50", "#2980b9", "#8e44ad"][:len(dn)]
fig, axes = plt.subplots(len(MAP), len(a.label), figsize=((2.0 + 0.45 * len(meths)) * len(a.label), 3.2 * len(MAP)), squeeze=False)
for i, cls in enumerate(MAP):
    for j, olab in enumerate(a.label):
        ax = axes[i, j]; r = R[(R.cls == cls) & (R.resolution == olab)].set_index("method").reindex(meths)
        for k, mth in enumerate(meths):
            v = r.loc[mth, "AP"]; ax.bar(k, v, 0.7, color=cols[k], yerr=[[v - r.loc[mth, "AP_lo"]], [r.loc[mth, "AP_hi"] - v]], capsize=2)
            ax.text(k, r.loc[mth, "AP_hi"], f"{v:.2f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(range(len(meths))); ax.set_xticklabels(meths, rotation=35, ha="right", fontsize=7)
        ax.set_ylim(0, R[R.cls == cls]["AP_hi"].max() * 1.25)
        if i == 0: ax.set_title(f"our background: {olab}", fontsize=9)
        if j == 0: ax.set_ylabel(f"{NAME[cls]}\naverage precision (95% CI)", fontsize=8)
fig.suptitle(f"Known breast cancer drivers (IntOGen) at gene level: our residual vs dNdScv on the same mutations\n"
             f"same gene universe and driver list; {int(R['n_drivers'].iloc[0])} drivers", fontsize=10)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_comparison_with_dndscv.{ext}", dpi=200, bbox_inches="tight")
print(f"Output in {a.out_dir}/")