#!/usr/bin/env python
"""
BURNER - criterio d'uso SEVERO rispetto al pavimento, calcolato a posteriori sui run del pilota (nessun nuovo modello).
=================================================================================================================
Richiede un run lanciato con --floor_importance. Legge importance_per_chr.tsv (modello vero, una riga per fold e feature),
importance_floor_per_rep.tsv (modelli del pavimento: target spostato, una riga per ripetizione x fold e feature),
feature_groups.tsv e usage_per_omic.tsv.

Per ogni omica e gruppo di feature (importanza del gruppo = somma delle importanze dei suoi membri, come nel pilota):
  floor_max          = la MASSIMA importanza del gruppo su TUTTE le repliche del pavimento (ripetizioni x fold)
  beyond_floor (fold) = importanza del gruppo nel modello vero, in quella fold, > floor_max
  frac_beyond_floor  = frazione delle fold del modello vero in cui vale beyond_floor
  used_beyond_floor  = frac_beyond_floor >= --min_frac (default 0.5, come --used_frac_chroms del pilota)
  used_strict        = used (sopra le esche, criterio del pilota) E used_beyond_floor
Nessuna soglia da scegliere sul rapporto: il riferimento e' l'intera distribuzione del pavimento (il suo massimo).
Con 3 ripetizioni x 5 fold = 15 repliche del pavimento: sotto l'ipotesi che il gruppo non sia legato all'occorrenza
(importanza vera scambiabile con quelle del pavimento), la probabilita' di superarle tutte in una fold e' 1/16.
Uscita: usage_per_omic_strict.tsv nella cartella del run (le colonne di usage_per_omic.tsv invariate + le nuove).
"""
import argparse, sys
import numpy as np, pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--run", nargs="+", required=True, help="cartelle dei run del pilota (con --floor_importance)")
p.add_argument("--min_frac", type=float, default=0.5)
a = p.parse_args()
IS_DECOY = lambda f: f.str.startswith("DECOY_") | f.str.startswith("SHADOW")

for rdir in a.run:
    try:
        I = pd.read_csv(f"{rdir}/importance_per_chr.tsv", sep="\t", dtype={"chr": str})
        F = pd.read_csv(f"{rdir}/importance_floor_per_rep.tsv", sep="\t", dtype={"chr": str})
    except FileNotFoundError as e:
        sys.exit(f"{rdir}: {e.filename} assente (il run deve essere lanciato con --floor_importance)")
    groups = pd.read_csv(f"{rdir}/feature_groups.tsv", sep="\t")[["feature", "group_label"]]
    U = pd.read_csv(f"{rdir}/usage_per_omic.tsv", sep="\t")
    gsum = lambda X: (X[~IS_DECOY(X["feature"])].merge(groups, on="feature")
                      .groupby(["omic", "chr", "group_label"], as_index=False)["imp"].sum())
    R, FL = gsum(I), gsum(F)
    fmax = FL.groupby(["omic", "group_label"])["imp"].agg(floor_max="max", floor_mean="mean", n_floor_reps="size").reset_index()
    R = R.merge(fmax, on=["omic", "group_label"], how="left")
    R["beyond_floor"] = R["imp"] > R["floor_max"]
    S = R.groupby(["omic", "group_label"]).agg(frac_beyond_floor=("beyond_floor", "mean"),
                                               real_mean_imp=("imp", "mean"), floor_max=("floor_max", "first"),
                                               floor_mean_imp=("floor_mean", "first"), n_floor_reps=("n_floor_reps", "first")).reset_index()
    S["used_beyond_floor"] = S["frac_beyond_floor"] >= a.min_frac
    out = U.merge(S, on=["omic", "group_label"], how="left")
    out["used_beyond_floor"] = out["used_beyond_floor"].fillna(False).astype(bool)
    out["used_strict"] = out["used"].astype(bool) & out["used_beyond_floor"]
    out.to_csv(f"{rdir}/usage_per_omic_strict.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    show = out[out["used"]].sort_values(["omic", "frac_beyond_floor"], ascending=[True, False])
    print(f"\n=== {rdir}\nrepliche del pavimento per gruppo: {int(S['n_floor_reps'].median())}\n"
          + show[["omic", "group_label", "used", "frac_beyond_floor", "real_mean_imp", "floor_max", "used_strict"]].round(3).to_string(index=False))
    print("gruppi usati (criterio del pilota) -> severi, per omica:\n" +
          out.groupby("omic")[["used", "used_strict"]].sum().astype(int).to_string())