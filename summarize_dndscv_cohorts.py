#!/usr/bin/env python
"""
BURNER - riassunto di dNdScv con e senza covariate per piu' coorti: quanto contano le covariate epigenetiche.
Per ogni coorte (cartella con sel_cv_none.tsv, sel_cv_<cv>.tsv e i log log_none.txt, log_cv.txt di run_dndscv.R):
campioni e sostituzioni, theta della regressione senza e con covariate (e guadagno relativo), geni con q < soglia
(tutte le sostituzioni) nei due casi, geni guadagnati e persi con le covariate.
Uso: python summarize_dndscv_cohorts.py --cohort SKCM analyses/dndscv_cohorts/SKCM --cohort LUSC ... [--cv hg19_hg38_epigenome_pcawg]
"""
import argparse, re, os
import pandas as pd
p = argparse.ArgumentParser()
p.add_argument("--cohort", nargs=2, action="append", metavar=("NOME", "CARTELLA"), required=True)
p.add_argument("--cv", default="hg19_hg38_epigenome_pcawg"); p.add_argument("--q", type=float, default=0.1)
a = p.parse_args()
def from_log(path):
    if not os.path.exists(path): return {}
    t = open(path).read(); r = {}
    m = re.search(r"theta = ([\d.]+)", t); r["theta"] = float(m.group(1)) if m else None
    m = re.search(r"sostituzioni singole uniche \(campione x mutazione\): (\d+), campioni: (\d+)", t)
    if m: r["subs"], r["samples"] = int(m.group(1)), int(m.group(2))
    m = re.search(r"(\d+) samples excluded", t); r["hyper_excluded"] = int(m.group(1)) if m else 0
    return r
rows, diffs = [], []
for name, d in a.cohort:
    L0, L1 = from_log(f"{d}/log_none.txt"), from_log(f"{d}/log_cv.txt")
    S0 = pd.read_csv(f"{d}/sel_cv_none.tsv", sep="\t").drop_duplicates("gene_name").set_index("gene_name")
    S1 = pd.read_csv(f"{d}/sel_cv_{a.cv}.tsv", sep="\t").drop_duplicates("gene_name").set_index("gene_name")
    g0, g1 = set(S0.index[S0["qallsubs_cv"] < a.q]), set(S1.index[S1["qallsubs_cv"] < a.q])
    th0, th1 = L0.get("theta"), L1.get("theta")
    rows.append(dict(coorte=name, campioni=L0.get("samples"), ipermutati_esclusi=L0.get("hyper_excluded"), sostituzioni=L0.get("subs"),
                     theta_senza=th0, theta_con=th1, guadagno_theta=(th1 / th0 - 1) if th0 and th1 else None,
                     geni_q_senza=len(g0), geni_q_con=len(g1), in_comune=len(g0 & g1), solo_senza=len(g0 - g1), solo_con=len(g1 - g0)))
    diffs.append(f"{name}: solo SENZA covariate: {', '.join(sorted(g0 - g1)) or '-'} | solo CON covariate: {', '.join(sorted(g1 - g0)) or '-'}")
T = pd.DataFrame(rows)
pd.set_option("display.width", 220); print(T.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\n" + "\n".join(diffs))
print("\nLettura: guadagno_theta = quanto le covariate riducono la variabilita' tra geni nella regressione sulle sinonime; "
      "solo_senza = geni significativi solo senza covariate (candidati falsi positivi che le covariate eliminano).")