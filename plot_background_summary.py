#!/usr/bin/env python
"""
BURNER - riassunto delle prestazioni del background per omica, alle diverse risoluzioni.
=========================================================================================
Legge, per ogni run del pilota, performance.tsv (modello XGBoost e lineare) e floor.tsv (pavimento: stesso modello
allenato sul target spostato circolarmente, che conserva l'autocorrelazione ma rompe il legame con le feature).
Figura: un pannello per omica; per ogni risoluzione la devianza spiegata out-of-fold del modello (punto pieno), del
modello lineare (grigio) e del pavimento (media e intervallo tra le ripetizioni). Quello che il modello impara oltre la
struttura spaziale e' la distanza tra il punto pieno e il pavimento. Ogni pannello ha la sua scala (le omiche hanno
livelli molto diversi). Valori del lineare sotto il limite inferiore del pannello: segnati come triangolo sul bordo.
Uso: python plot_background_summary.py --run analyses/v2/20000bp/pilot_classes/genome_nopos 20kb ... --out_dir DIR
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
p = argparse.ArgumentParser()
p.add_argument("--run", nargs=2, action="append", metavar=("DIR", "LABEL"), required=True)
p.add_argument("--metric", default="dev_expl", choices=["dev_expl", "R2"])
p.add_argument("--out_dir", default="background_summary")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
NAME = {"cna_ampl": "CNA amplifications", "cna_del": "CNA deletions", "atac": "Chromatin accessibility",
        "mut_syn": "Synonymous mutations", "mut_mis": "Missense mutations", "mut_trunc": "Truncating mutations", "mut": "Mutations (old target)"}
lin_col = {"dev_expl": "dev_expl_lin", "R2": "R2_ridge"}[a.metric]
rows = []
for d, lab in a.run:
    P = pd.read_csv(f"{d}/performance.tsv", sep="\t"); P = P[P["unit"] == "genome"] if "unit" in P else P
    F = pd.read_csv(f"{d}/floor.tsv", sep="\t") if os.path.exists(f"{d}/floor.tsv") else pd.DataFrame(columns=["omic", a.metric])
    for _, r in P.iterrows():
        f = F[F["omic"] == r["omic"]][a.metric].astype(float)
        rows.append(dict(resolution=lab, omic=r["omic"], family=r.get("family", ""), n_bins=r.get("n_bins", np.nan), model=r[a.metric], linear=r[lin_col],
                         floor_mean=f.mean() if len(f) else np.nan, floor_min=f.min() if len(f) else np.nan,
                         floor_max=f.max() if len(f) else np.nan))
T = pd.DataFrame(rows); T["model_minus_floor"] = T["model"] - T["floor_mean"]
T.to_csv(f"{a.out_dir}/background_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(T.round(3).to_string(index=False))
labels = [lab for _, lab in a.run]; oms = [o for o in NAME if o in set(T["omic"])] + sorted(set(T["omic"]) - set(NAME))
nc = 3; nr = int(np.ceil(len(oms) / nc))
fig, axes = plt.subplots(nr, nc, figsize=(4.3 * nc, 3.3 * nr), squeeze=False)
x = np.arange(len(labels))
for k, om in enumerate(oms):
    ax = axes[k // nc, k % nc]; t = T[T["omic"] == om].set_index("resolution").reindex(labels)
    lo_vals = np.r_[t["model"].to_numpy(float), t["floor_min"].to_numpy(float), 0.0]
    hi_vals = np.r_[t["model"].to_numpy(float), t["floor_max"].to_numpy(float), t["linear"].clip(upper=None).to_numpy(float)]
    lo, hi = np.nanmin(lo_vals), np.nanmax(hi_vals); pad = 0.12 * (hi - lo if hi > lo else 0.05); lo, hi = lo - pad, hi + pad
    ax.fill_between(x, t["floor_min"], t["floor_max"], color="#e6b0aa", alpha=0.6, lw=0, label="floor (shifted target), range")
    ax.plot(x, t["floor_mean"], "--", color="#c0392b", lw=1.2, label="floor, mean")
    lin = t["linear"].to_numpy(float); inside = lin >= lo
    ax.plot(x[inside], lin[inside], "s", color="0.55", ms=5, label="linear model")
    if (~inside).any(): ax.plot(x[~inside], np.full((~inside).sum(), lo), "v", color="0.55", ms=6, clip_on=False, label="linear (below axis)")
    ax.plot(x, t["model"], "o-", color="#1f618d", ms=7, lw=1.5, label="XGBoost (out-of-fold)")
    ax.axhline(0, color="k", lw=0.6); ax.set_ylim(lo, hi)
    ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_title(NAME.get(om, om), fontsize=10)
    # NOTE (28/10/2026): con la perdita gaussiana (accessibilita') la devianza spiegata E' l'R2 (1 - RSS/TSS): lo si scrive.
    #      Per le altre perdite: devianza spiegata rispetto al modello di riferimento (solo intercetta o, per le mutazioni,
    #      solo l'offset dell'esposizione).
    fam = str(t["family"].dropna().iloc[0]) if t["family"].notna().any() else ""
    if a.metric == "R2" or fam == "gaussian": ylab = "R² (out-of-fold)"
    elif om.startswith("mut"): ylab = "deviance explained\n(vs exposure-only model)"
    else: ylab = "deviance explained\n(vs intercept-only model)"
    ax.set_ylabel(ylab, fontsize=8)
    if k == 0: ax.legend(fontsize=6.5, frameon=False, loc="best")
for k in range(len(oms), nr * nc): axes[k // nc, k % nc].axis("off")
fig.suptitle("Agnostic background model: same features and procedure for every alteration type\n"
             "what the model learns beyond spatial structure = distance between the blue points and the red floor", fontsize=10)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_background_summary.{ext}", dpi=200, bbox_inches="tight")
print(f"Output in {a.out_dir}/")