"""
Gather the label randomizations of pilot_predictor_identity.py (run_identity_randomizations.sh) into single tables,
with a 'shuffle' column (Gab, 20/11/2026). Only plumbing: no statistic and no test is computed here.

Inputs
  --rand_dir   folder with shuffle_<k>/ subfolders (one per randomization)
  --truth_dir  optional: the truth run (same dataset / labels / seeds / block_bins / NA rule); added with shuffle = 0
Outputs (in --rand_dir)
  all_base_models.tsv, all_removal_deltas.tsv, all_addition_deltas.tsv, all_swap_deltas_summary.tsv,
  all_cluster_labels.tsv (cluster, representative, cluster label per omic, per shuffle)
Checks (fail loudly): incomplete shuffles are listed and left out; every shuffle must have the same clusters as the
first one (same features, same cluster ids), otherwise the runs are not comparable and the script stops.

Usage (cluster, pilot env)
  python collect_identity_randomizations.py --rand_dir analyses/v4/BRCA/1000000bp/predictor_identity_na0_rand \
      --truth_dir analyses/v4/BRCA/1000000bp/predictor_identity_na0/<omic>
  (the truth runs of 20/11/2026 have one folder per omic: pass them all, --truth_dir takes several folders)
"""
import argparse, glob, os, re, sys
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--rand_dir", required=True)
p.add_argument("--truth_dir", nargs="*", default=[], help="truth run folder(s), added as shuffle 0")
a = p.parse_args()

TABLES = ["base_models", "removal_deltas", "addition_deltas", "swap_deltas_summary"]

runs = []                                                       # (shuffle, folder)
for d in sorted(glob.glob(f"{a.rand_dir}/shuffle_*")):
    m = re.fullmatch(r"shuffle_(\d+)", os.path.basename(d))
    if m: runs.append((int(m.group(1)), d))
if not runs:
    sys.exit(f"no shuffle_<k> folders in {a.rand_dir}: did run_identity_randomizations.sh write there?")
complete = [(k, d) for k, d in runs if all(os.path.isfile(f"{d}/{t}.tsv") for t in TABLES + ["clusters"])]
incomplete = sorted(k for k, d in runs if (k, d) not in complete)
if incomplete:
    print(f"WARNING: {len(incomplete)} incomplete shuffles left out (still running or failed; check their logs): {incomplete}")
runs = [(0, d) for d in a.truth_dir] + complete
for d in a.truth_dir:
    missing = [t for t in TABLES + ["clusters"] if not os.path.isfile(f"{d}/{t}.tsv")]
    if missing: sys.exit(f"truth folder {d}: missing {missing}")

# the clusters must be identical in every run (they do not depend on the labels; a difference means another dataset or
# another NA rule)
ref_k, ref_d = runs[0]
ref = pd.read_csv(f"{ref_d}/clusters.tsv", sep="\t")[["feature", "cluster"]].sort_values("feature").reset_index(drop=True)
out = {t: [] for t in TABLES}
labels = []
for k, d in runs:
    C = pd.read_csv(f"{d}/clusters.tsv", sep="\t")
    same = C[["feature", "cluster"]].sort_values("feature").reset_index(drop=True)
    if not same.equals(ref):
        sys.exit(f"{d}: clusters differ from {ref_d}: not the same dataset / NA rule; runs are not comparable")
    label_cols = [c for c in C.columns if c.startswith("cluster_label_")]
    L = C.groupby("cluster").first()[["representative"] + label_cols].reset_index()
    L.insert(0, "shuffle", k); labels.append(L)
    for t in TABLES:
        T = pd.read_csv(f"{d}/{t}.tsv", sep="\t")
        T.insert(0, "shuffle", k); out[t].append(T)

for t in TABLES:
    pd.concat(out[t], ignore_index=True).to_csv(f"{a.rand_dir}/all_{t}.tsv", sep="\t", index=False)
# truth folders split by omic share the same clusters -> keep one label row per (shuffle, cluster), merging omic columns
LAB = pd.concat(labels, ignore_index=True).groupby(["shuffle", "cluster"], as_index=False).first()
LAB.to_csv(f"{a.rand_dir}/all_cluster_labels.tsv", sep="\t", index=False)
n_rand = len([k for k, _ in runs if k > 0])
print(f"gathered {n_rand} randomizations" + (f" + truth ({len(a.truth_dir)} folders)" if a.truth_dir else "")
      + f" into {a.rand_dir}/all_*.tsv")
