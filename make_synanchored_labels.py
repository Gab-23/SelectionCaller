"""
Predictor identity: labels anchored on a neutral reference (Gab, 23/11/2026).

Idea (Gab, 23/11/2026): synonymous mutations are (close to) neutral, so the labels are DEFINED by behaviour on them
alone, with no prior: a cluster that predicts synonymous occurrence is O, every other cluster is S. The labels are then
transferred unchanged to every omic (tune on one omic, generalize) and validated by the identity pilot and its label
randomizations on the OTHER targets (missense, truncating, CNA, ATAC), never on the anchor.
S therefore means "no detectable occurrence signal on synonymous mutations": it includes features without any signal
(noise), which this rule cannot separate from selection.

Rule (Gab's definition; threshold and importance chosen by Claude, reported 23/11/2026):
  importance  = dev_single of the anchor omic (pilot_predictor_identity.py --single_cluster_importance): held-out
                deviance explained of the model with that cluster alone, mean over seeds (NOTE 23/11/2026: not the
                drop-cluster importance, which misses occurrence predictors that are redundant with each other; Gab:
                the label is "does it predict synonymous mutations", so the direct measure)
  noise       = SE = pooled seed sd over all clusters (sqrt of the mean per-cluster variance) / sqrt(n seeds)
  cluster O if mean dev_single > 2 SE, otherwise S; every feature takes the label of its cluster, for every omic
  (the clusters are then never mixed). The prior labels (--labels) are used only to report what changed.

Inputs
  --importance  single_cluster_importance.tsv of the anchor run (same dataset, NA rule and labels as the --clusters run)
  --clusters    clusters.tsv of the same run
  --labels      prior labels (config/predictor_identity_labels.tsv): feature + one column per omic (O / S / B); only for
                the report and for the list of omic columns
  --anchor      anchor omic (default mut_syn)
Outputs
  --out         new labels, same format as --labels (feature + one column per omic, only O / S)
  <out>.report.tsv  one row per cluster: representative, members, dev_single mean / sd, SE, occurrence-like, and the
                    features whose label changed, per omic
Usage (cluster, pilot env)
  python make_synanchored_labels.py --importance <run>/single_cluster_importance.tsv --clusters <run>/clusters.tsv \
      --labels config/predictor_identity_labels.tsv --out <run>/labels_synanchored.tsv
"""
import argparse, sys
import numpy as np, pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--importance", required=True)
p.add_argument("--clusters", required=True)
p.add_argument("--labels", required=True)
p.add_argument("--anchor", default="mut_syn")
p.add_argument("--n_se", type=float, default=2.0, help="occurrence-like if mean dev_single > n_se * SE")
p.add_argument("--out", required=True)
a = p.parse_args()

I = pd.read_csv(a.importance, sep="\t")
I = I[I.omic == a.anchor]
if I.empty:
    sys.exit(f"{a.importance}: no rows for the anchor omic {a.anchor}: was the run made with --single_cluster_importance "
             f"and --omics including {a.anchor}?")
C = pd.read_csv(a.clusters, sep="\t")
L = pd.read_csv(a.labels, sep="\t", dtype=str).set_index("feature")
omics = [c for c in L.columns]

# importance per cluster (mean and sd over seeds) and the pooled noise of the anchor model
per_cluster = I.groupby(["cluster", "representative"]).dev_single.agg(["mean", "std", "count"]).reset_index()
n_seeds = int(per_cluster["count"].iloc[0])
if (per_cluster["count"] != n_seeds).any():
    sys.exit(f"{a.importance}: clusters with a different number of seeds: {per_cluster[per_cluster['count'] != n_seeds]}")
se = np.sqrt((per_cluster["std"] ** 2).mean()) / np.sqrt(n_seeds)
per_cluster["se"] = se
per_cluster["occurrence_like"] = per_cluster["mean"] > a.n_se * se
missing = sorted(set(C.cluster) - set(per_cluster.cluster))
if missing:
    sys.exit(f"clusters of {a.clusters} without importance rows: {missing}: importance and clusters from different runs?")
occ_like = dict(zip(per_cluster.cluster, per_cluster.occurrence_like))
cluster_of = dict(zip(C.feature, C.cluster))
not_in_run = [f for f in L.index if f not in cluster_of]
print(f"anchor {a.anchor}: {len(per_cluster)} clusters, SE {se:.5f} ({n_seeds} seeds), threshold {a.n_se * se:+.5f}; "
      f"occurrence-like: {int(per_cluster.occurrence_like.sum())} clusters "
      f"({', '.join(per_cluster.loc[per_cluster.occurrence_like, 'representative'])})")
if not_in_run:
    print(f"NOTE: {len(not_in_run)} features of {a.labels} are not in the run (dropped features): kept with their prior "
          f"label, the pilot ignores them: {not_in_run}")

# labels from the anchor only: O if the cluster is occurrence-like, S otherwise; the same for every omic column
NEW = L.copy()
changes = {om: [] for om in omics}
for om in omics:
    for f in L.index:
        if f not in cluster_of: continue
        prior = L.loc[f, om].strip().upper()
        new = "O" if occ_like[cluster_of[f]] else "S"
        NEW.loc[f, om] = new
        if new != prior: changes[om].append(f"{f}:{prior}->{new}")
    print(f"{om}: {len(changes[om])} labels differ from the prior: {changes[om]}")
NEW.reset_index().to_csv(a.out, sep="\t", index=False)

members = C.groupby("cluster").feature.apply(lambda s: ", ".join(sorted(s)))
per_cluster["members"] = per_cluster.cluster.map(members)
for om in omics:
    per_cluster[f"changed_{om}"] = per_cluster.cluster.map(
        lambda k: ", ".join(c for c in changes[om] if cluster_of.get(c.split(":")[0]) == k))
per_cluster.rename(columns={"mean": "dev_single_mean", "std": "dev_single_sd", "count": "n_seeds"}) \
    .sort_values("dev_single_mean", ascending=False).to_csv(f"{a.out}.report.tsv", sep="\t", index=False)
print(f"written {a.out} and {a.out}.report.tsv")
