#!/usr/bin/env python
"""
BURNER - Venn dei driver noti ritrovati da noi e da dNdScv, dentro l'universo dei driver (rettangolo).
======================================================================================================
Rettangolo = tutti i driver IntOGen del nostro universo di geni (colonna is_driver dei nostri gene_scores_<classe>.tsv).
Cerchi = driver "ritrovati" da noi e da dNdScv. Fuori dai cerchi = driver non ritrovati da nessuno.
Criterio "ritrovato" (--mode):
  q    (default): q-value < --q per entrambi. dNdScv: la sua colonna q (qallsubs_cv / qmis_cv / qtrunc_cv).
                  Noi: Benjamini-Hochberg sui nostri p-value, p = 10^(-punteggio) (coda della binomiale negativa),
                  su tutti i geni dell'universo (geni senza punteggio: p = 1).
  top           : i primi N geni di ciascuna classifica, con N = numero di geni significativi di dNdScv (q < --q):
                  confronta le classifiche a parita' di numero di geni chiamati, senza dipendere dalla calibrazione.
In basso a destra: quanti NON driver ciascun metodo chiama (possibili scoperte o falsi positivi).
Uso: python venn_drivers_vs_dndscv.py --ours_dir analyses/v2/20000bp/pilot_classes/drivers_gene_level --cls nonsyn \
       --dndscv analyses/dndscv_BRCA/sel_cv_hg19_hg38_epigenome_pcawg.tsv [--score residual|<variante>] [--mode q|top] --out_dir DIR
"""
import argparse, os, textwrap
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
p = argparse.ArgumentParser()
p.add_argument("--ours_dir", required=True); p.add_argument("--cls", default="nonsyn", choices=["mis", "trunc", "nonsyn"])
p.add_argument("--dndscv", required=True); p.add_argument("--dndscv_label", default="dNdScv")
p.add_argument("--score", default="residual", help="residual (nostro modello) o una colonna di gene_scores_variants_<classe>.tsv (es. model+syn)")
p.add_argument("--mode", default="q", choices=["q", "top"]); p.add_argument("--q", type=float, default=0.1)
p.add_argument("--label", default=""); p.add_argument("--out_dir", default="venn_drivers")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
QCOL = {"mis": "qmis_cv", "trunc": "qtrunc_cv", "nonsyn": "qallsubs_cv"}
PCOL = {"mis": "pmis_cv", "trunc": "ptrunc_cv", "nonsyn": "pallsubs_cv"}
G = pd.read_csv(f"{a.ours_dir}/gene_scores_{a.cls}.tsv", sep="\t", index_col=0)
if a.score == "residual": s = G["residual"]
else: s = pd.read_csv(f"{a.ours_dir}/gene_scores_variants_{a.cls}.tsv", sep="\t", index_col=0)[a.score].reindex(G.index)
D = pd.read_csv(a.dndscv, sep="\t").drop_duplicates("gene_name").set_index("gene_name")
U = G.index.intersection(D.index); drivers = set(U[G.loc[U, "is_driver"].astype(bool)])
pv = np.power(10.0, -s.reindex(U).astype(float).fillna(0.0).to_numpy())          # senza punteggio -> p = 1
def bh(p_):
    n = len(p_); o = np.argsort(p_); q = np.empty(n); q[o] = np.minimum.accumulate((p_[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)
q_ours = pd.Series(bh(pv), index=U); q_dn = D.loc[U, QCOL[a.cls]].astype(float)
n_dn = int((q_dn < a.q).sum())
if a.mode == "q":
    ours_hit = set(U[q_ours < a.q]); dn_hit = set(U[q_dn < a.q]); crit = f"q < {a.q} (ours: Benjamini-Hochberg on our p-values)"
else:
    ours_hit = set(pd.Series(pv, index=U).sort_values().index[:n_dn]); dn_hit = set(q_dn.sort_values().index[:n_dn])
    crit = f"top {n_dn} genes of each ranking (N = genes with dNdScv q < {a.q})"
both, only_o, only_d = drivers & ours_hit & dn_hit, (drivers & ours_hit) - dn_hit, (drivers & dn_hit) - ours_hit
none_ = drivers - ours_hit - dn_hit
nd_o, nd_d = sorted(ours_hit - drivers), sorted(dn_hit - drivers)
rank = lambda S, key: sorted(S, key=key)
k_o = lambda g: -s.get(g, -np.inf) if pd.notna(s.get(g)) else np.inf
k_d = lambda g: D.loc[g, PCOL[a.cls]]
fig, ax = plt.subplots(figsize=(9.5, 6.2))
ax.add_patch(Rectangle((-3.0, -2.15), 6.0, 4.3, facecolor="#f4f4f4", edgecolor="0.4", lw=1.2))
ax.add_patch(Circle((-0.75, 0.0), 1.4, facecolor="#c0392b", alpha=0.18, edgecolor="#c0392b", lw=1.6))
ax.add_patch(Circle((0.75, 0.0), 1.4, facecolor="#2c3e50", alpha=0.15, edgecolor="#2c3e50", lw=1.6))
def put(x, y, genes, bold=False, fs=7.4):
    txt = "\n".join(genes) if genes else "—"
    ax.text(x, y, txt, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", linespacing=1.2)
put(-1.45, 0.0, rank(only_o, k_o)); put(1.45, 0.0, rank(only_d, k_d)); put(0, 0.0, rank(both, k_d), bold=True)
ax.text(-1.55, 1.25, f"ours ({a.score})\n{len(drivers & ours_hit)} drivers", ha="right", va="bottom", fontsize=9, color="#c0392b", fontweight="bold")
ax.text(1.55, 1.25, f"{a.dndscv_label}\n{len(drivers & dn_hit)} drivers", ha="left", va="bottom", fontsize=9, color="#2c3e50", fontweight="bold")
ax.text(-3.0, 2.22, f"rectangle = all known drivers in the gene universe (IntOGen): {len(drivers)}", ha="left", va="bottom", fontsize=8.5)
nn = sorted(none_)
ax.text(-2.9, -2.05, textwrap.fill(f"found by neither ({len(nn)}): " + ", ".join(nn), 105), ha="left", va="bottom", fontsize=6.8)
ax.text(2.9, 1.95, f"non-driver genes called\nours: {len(nd_o)}   {a.dndscv_label}: {len(nd_d)}", ha="right", va="top", fontsize=7.5, color="0.3")
ax.set_xlim(-3.1, 3.1); ax.set_ylim(-2.3, 2.6); ax.set_aspect("equal"); ax.axis("off")
cn = {"mis": "missense", "trunc": "truncating", "nonsyn": "non-synonymous"}[a.cls]
ax.set_title(f"Known breast cancer drivers found by our method and by {a.dndscv_label} ({cn}{', ' + a.label if a.label else ''})\n"
             f"called = {crit}; within each region ordered by significance", fontsize=9.5)
tag = f"{a.cls}_{a.score.replace('+', 'plus')}_{a.mode}"
for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_venn_drivers_{tag}.{ext}", dpi=200, bbox_inches="tight")
rows = [dict(gene=g, is_driver=g in drivers, ours_called=g in ours_hit, dndscv_called=g in dn_hit, ours_q=q_ours.get(g), dndscv_q=q_dn.get(g))
        for g in sorted(drivers | ours_hit | dn_hit)]
pd.DataFrame(rows).to_csv(f"{a.out_dir}/venn_drivers_{tag}.tsv", sep="\t", index=False)
print(f"driver: {len(drivers)} | entrambi {len(both)} | solo nostri {len(only_o)} | solo {a.dndscv_label} {len(only_d)} | nessuno {len(none_)}")
print(f"non driver chiamati: nostri {len(nd_o)} ({', '.join(nd_o[:15])}{'...' if len(nd_o) > 15 else ''}) | {a.dndscv_label} {len(nd_d)} ({', '.join(nd_d[:15])})")
print(f"Output in {a.out_dir}/")