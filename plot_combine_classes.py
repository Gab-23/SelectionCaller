#!/usr/bin/env python
"""
BURNER - Figura dei risultati di combine_classes_test.py: AP sui driver per risoluzione, una colonna per coorte.
Linee: somma delle classi (oggi) e combinazioni (fisher, minp, lrt); riferimento tratteggiato: dNdScv.
Legge le tabelle combine_classes_<COORTE>_<risoluzione>.tsv scritte da combine_classes_test.py (--label "COORTE risoluzione").
Uso: python plot_combine_classes.py --in_dir analyses/combine_classes --cohorts BRCA LUAD --out_dir analyses/combine_classes
"""
import argparse, glob, os, re
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
p = argparse.ArgumentParser()
p.add_argument("--in_dir", required=True); p.add_argument("--cohorts", nargs="+", default=["BRCA", "LUAD"])
p.add_argument("--out_dir", default=None)
a = p.parse_args(); out = a.out_dir or a.in_dir; os.makedirs(out, exist_ok=True)
STYLE = {  # metodo: (etichetta, colore, stile)
    "sum of classes (model)":     ("sum of classes (model)",     "#9e9e9e", "-"),
    "sum of classes (model+syn)": ("sum of classes (model+syn)", "#616161", "-"),
    "fisher (model+syn)":         ("Fisher (model+syn)",         "#1f77b4", "-"),
    "minp (model+syn)":           ("min-p (model+syn)",          "#d62728", "-"),
    "lrt (model)":                ("LRT, dNdScv-like (model)",   "#2ca02c", "-"),
}
RES_LABEL = {20000: "20 kb", 100000: "100 kb", 500000: "500 kb", 1000000: "1 Mb"}
fig, axes = plt.subplots(1, len(a.cohorts), figsize=(5.2 * len(a.cohorts), 4.3), sharey=True)
axes = [axes] if len(a.cohorts) == 1 else axes
for ax, coh in zip(axes, a.cohorts):
    T = []
    for f in glob.glob(f"{a.in_dir}/combine_classes_{coh}_*.tsv"):
        W = int(re.search(rf"combine_classes_{coh}_(\d+)\.tsv$", f).group(1)); T.append(pd.read_csv(f, sep="\t").assign(W=W))
    if not T: ax.set_title(f"{coh}: no tables"); continue
    T = pd.concat(T); Ws = sorted(T.W.unique()); x = range(len(Ws))
    for m, (lab, col, ls) in STYLE.items():
        s = T[T.method == m].set_index("W").reindex(Ws)["AP"]
        ax.plot(x, s.values, ls, color=col, marker="o", lw=2 if "sum" not in m else 1.5, label=lab)
    d = T[T.method == "dNdScv"].set_index("W").reindex(Ws)["AP"]
    ax.plot(x, d.values, "--", color="black", lw=1.5, label="dNdScv")
    ax.set_xticks(list(x)); ax.set_xticklabels([RES_LABEL.get(w, str(w)) for w in Ws])
    n_drv = None
    ax.set_title(f"{coh}"); ax.set_xlabel("background resolution"); ax.grid(alpha=0.3)
axes[0].set_ylabel("average precision (known drivers)")
fig.legend(*axes[0].get_legend_handles_labels(), fontsize=8.5, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
fig.suptitle("Non-synonymous driver detection: summing vs combining missense and truncating evidence", fontsize=11)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(f"{out}/fig_combine_classes.{ext}", dpi=200, bbox_inches="tight")
print(f"Figura in {out}/fig_combine_classes.png")