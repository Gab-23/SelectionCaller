"""
PILOT - identity of the predictors: occurrence (O) vs selection (S), per omic, with the pilot's XGBoost models.

Idea (Gab, 19/11/2026):
  1. Prior labels O / S per feature and per omic from domain knowledge (--labels), taken as ground truth.
  2. Redundancy: hierarchical clustering of the features in Spearman space (distance 1 - |rho|, complete linkage,
     cut at --redundancy_rho: every pair inside a cluster has |rho| >= that value). Each cluster is represented by its
     medoid (the member with the highest mean |rho| to the other members). A cluster whose members have different
     labels for an omic ("mixed") is labelled O for that omic ("occurrence until proven selection").
  3. Swap test: for each omic, an O-only model (all O clusters) and an S-only model (all S clusters). For each cluster A
     of the model, A is removed, and then replaced by each cluster B of the other class.
  4. Every swap is decomposed as   d_swap = d_remove(A) + d_add(B | A removed)   (deltas of held-out deviance
     explained), so that the distortion due to losing A is separated from the one due to adding B.
     How the deltas become a confidence score is NOT implemented here: this script only writes the deltas.
  5. NOTE (19/11/2026): each cluster B of the other class is also added to the FULL model of the class (A present):
     d_add_full(B) = dev(class + B) - dev(class). Then  shared(A, B) = d_add(B | A removed) - d_add_full(B)
     separates "B stands in for A" (B adds value only when A is missing: A and B share information about the target)
     from "B brings its own signal" (B adds the same with or without A).
  6. NOTE (20/11/2026): --label_shuffle_seed permutes the cluster labels (same permutation for every omic, class sizes
     in clusters kept): one randomization of the null for the concordance test. Many seeds are run as a SLURM array by
     run_identity_randomizations.sh and gathered by collect_identity_randomizations.py. The test itself is not here.
  7. NOTE (20/11/2026): O+S model (all clusters) and its null: the same model with the S block circularly shifted
     within each chromosome (--n_shift repeats): does S add to O in combination, beyond what any S-like block would?
     --full_models_only runs only the base models and this block (no swap test), for setups whose swaps already exist.

Models: the same XGBoost regressors as pilot_omics_background.py --loss glm (fixed hyperparameters, 300 trees), with
the family of each omic: CNA binomial (fraction of patients), ATAC gaussian (logFC), mutation classes Poisson with
offset log(possible sites of the class). NOTE (19/11/2026): first sketched with a penalized GLM; switched to XGBoost
on Gab's request. Every fit is repeated with --seeds (XGBoost subsampling is random: on toy data a single seed moved
d_swap by up to 0.02).
CV: interleaved blocks with buffer (INTERLEAVED_KFOLD_IDXS), one model for the whole genome, same folds for every fit of
an omic. Metric: deviance explained on the pooled out-of-fold predictions, against the null fitted on the training folds.

Feature set: the same as pilot_omics_background.py (no alteration columns, no exposure columns, no cds_bp, no
dist.to.closest.OG / TSG which stay for validation; distance NAs capped at 2 Mb; optionally (--gene_na_to_zero,
20/11/2026) NA of gene-level features set to 0; features with >= 50% NA dropped;
Hi-C eigenvectors oriented by GC per chromosome). Constant features are dropped (logged).

Inputs
  --dataset   feature matrix (data/feature_matrix/<dataset_name>)
  --labels    TSV: column 'feature' + one column per omic (values O, S or B = both; B is treated as O, "occurrence until
              proven selection", and kept as B in clusters.tsv). Every feature of the feature set needs a
              label for every requested omic; unknown features in the file are reported and ignored.
Outputs (in --out_dir)
  clusters.tsv            feature, cluster, representative, mean |rho| to the cluster, per-omic feature and cluster labels
  base_models.tsv         per omic, class and seed: number of clusters, held-out deviance explained of the full model
  swap_deltas.tsv         one row per (omic, model class, removed cluster A, added cluster B, seed):
                          dev_base, dev_removed, dev_swapped, d_remove, d_add, d_swap, d_add_full, shared
                          (paired: same seed in every fit)
  swap_deltas_summary.tsv the same without the seed: mean and sd over seeds of d_remove, d_add, d_swap, d_add_full, shared
  addition_deltas.tsv     one row per (omic, model class, added cluster B, seed): dev_base, dev_full_plus_B, d_add_full
  removal_deltas.tsv      one row per (omic, model class, removed cluster A, seed): dev_base, dev_removed, d_remove
  full_models.tsv         NOTE (20/11/2026): one row per (omic, model, shift_rep, seed): model "O+S" (all clusters,
                          shift_rep -1) or "O+S_shifted" (S block circularly shifted, --n_shift repeats); dev_O, dev_S,
                          dev_expl, gain_over_O = dev_expl - dev_O (paired by seed)
  label_template.tsv      only with --write_label_template: the feature set with its clusters and empty label columns

Usage (cluster, pilot env; start with the 1 Mb dataset: the number of fits is 2 * n_O * n_S + 2 * (n_O + n_S) + 2 per omic)
  python pilot_predictor_identity.py --dataset <feature_matrix> --write_label_template --out_dir identity_out
  python pilot_predictor_identity.py --dataset <feature_matrix> --labels <labels.tsv> --out_dir identity_out \
      --omics cna_ampl cna_del mut_mis --bin_size 1000000
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import xgboost as xgb
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

sys.path.insert(0, os.getcwd())
try:
    from workflow.scripts.helpers.INTERLEAVED_KFOLD_IDXS import INTERLEAVED_KFOLD_IDXS
except ImportError:
    sys.exit("workflow/scripts/helpers/INTERLEAVED_KFOLD_IDXS.py is needed: run from the root of the repository.")

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--labels", default=None, help="TSV: feature + one column per omic with O / S / B (B treated as O)")
p.add_argument("--write_label_template", action="store_true", help="write label_template.tsv (feature set + clusters) and stop")
p.add_argument("--out_dir", default="identity_out")
p.add_argument("--omics", nargs="+", default=["cna_ampl", "cna_del", "atac", "mut_syn", "mut_mis", "mut_trunc"])
p.add_argument("--redundancy_rho", type=float, default=0.5, help="cut of the tree: |rho| >= this inside a cluster (complete linkage)")
p.add_argument("--n_splits", type=int, default=5)
p.add_argument("--block_bins", type=int, default=100, help="size of the interleaved blocks (bins)")
p.add_argument("--max_buffer_bins", type=int, default=100, help="cap of the buffer estimated from the target autocorrelation (bins)")
p.add_argument("--min_exonic_bp", type=float, default=500, help="mutation bins: at least this many coding bp (as in the pilot)")
p.add_argument("--n_trees", type=int, default=300, help="boosting rounds (as in the pilot)")
p.add_argument("--n_jobs", type=int, default=4)
p.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3],
               help="XGBoost seeds (subsample / colsample): every fit is repeated with each seed and the deltas are computed "
                    "seed by seed (paired: base, removal and swap with the same seed)")
p.add_argument("--gene_na_to_zero", action="store_true",
               help="NOTE (20/11/2026): set NA to 0 in the gene-level features (no gene in the bin) before the >= 50%% NA "
                    "filter (default off = previous behaviour)")
p.add_argument("--label_shuffle_seed", type=int, default=None,
               help="NOTE (20/11/2026): permute the cluster labels with this seed (null of the concordance test; one "
                    "permutation shared by all omics, class sizes in clusters kept). Default: true labels")
p.add_argument("--n_shift", type=int, default=20,
               help="NOTE (20/11/2026): circular shifts of the S block for the null of the O+S model (0 = O+S model only)")
p.add_argument("--shift_seed", type=int, default=1, help="seed of the circular shifts (same shifts for every omic run)")
p.add_argument("--full_models_only", action="store_true",
               help="NOTE (20/11/2026): fit only the base models (O-only, S-only) and the O+S model with its shifted-S "
                    "null; skip the swap test (removals, additions, swaps) and write only base_models.tsv, "
                    "full_models.tsv and clusters.tsv. Use a separate --out_dir from the full runs")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
t0 = time.time()
def log(msg): print(f"[{time.time() - t0:.0f}s] {msg}", flush=True)

# ============================ data and feature set (same rules as pilot_omics_background.py) ============================
D = pd.read_csv(a.dataset, sep="\t")
D["chr"] = D["binID"].str.split("_").str[0]
D["idx"] = D["binID"].str.split("_").str[1].astype(int)
chrom_order = lambda v: int(v) if str(v).isdigit() else 100 + ord(str(v)[0])
D = D.sort_values(["chr", "idx"], key=lambda s: s.map(chrom_order) if s.name == "chr" else s).reset_index(drop=True)

ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
VALIDATION_ONLY = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_"))]
EXONIC = ["cds_bp"] if "cds_bp" in D.columns else []
FEATS = [c for c in D.columns if c not in ALT + VALIDATION_ONLY + EXPOSURE + EXONIC + ["binID", "chr", "idx"]]
# NOTE (20/11/2026): with --gene_na_to_zero, NA of the gene-level features (= no gene, or no gene with that annotation, in
# the bin) becomes 0 BEFORE the >= 50% NA filter, so that at 20 kb they are not dropped (Gab, 20/11/2026). Post-hoc here
# for time reasons; TODO: move it to a preprocessing rule of the ML pipeline (the dataset files are never modified).
# Checked on BRCA 20 kb / 1 Mb and LUAD 1 Mb: the count features (n_*_bin, complex_density, n_*_trans) are NA only in
# bins without protein-coding genes (0 = exact); ohnolog/paralog means and totals are NA also where the bin has genes but
# none with ohnologs/paralogs (observed minimum 1, so 0 = "none" is exact); ess_score (0 = non-essential in DepMap),
# haploinsufficiency_score (observed minimum 0.2) and GTEx_* (centred values: 0 = average, NOT "absent") are
# approximations. ChromHMM, Hi-C eigenvectors and replication timing are NOT filled: their NA = missing measurement.
GENE_LEVEL_NA_TO_ZERO = [c for c in FEATS
                         if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                                  "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
                         or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs",
                                          "total_n_paralogs"))]
if a.gene_na_to_zero:
    n_filled = {c: int(D[c].isna().sum()) for c in GENE_LEVEL_NA_TO_ZERO}
    for c in GENE_LEVEL_NA_TO_ZERO:
        D[c] = D[c].fillna(0)
    log(f"--gene_na_to_zero: NA set to 0 in {len(GENE_LEVEL_NA_TO_ZERO)} gene-level features (NA count per feature): {n_filled}")
for c in [c for c in FEATS if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)                                   # NA on a whole chromosome = no gene of that list: cap
too_many_na = [c for c in FEATS if D[c].isna().mean() >= 0.5]
constant = [c for c in FEATS if c not in too_many_na and D[c].nunique() <= 1]
FEATS = [c for c in FEATS if c not in too_many_na + constant]
log(f"{len(D)} bins, {len(FEATS)} features; dropped: {len(too_many_na)} with >= 50% NA {too_many_na}, "
    f"{len(constant)} constant {constant}")

# Hi-C eigenvectors: arbitrary sign per chromosome -> oriented so that their correlation with GC is positive (as the pilot)
eig = [c for c in FEATS if c.startswith("eigen_")]
if eig and "GC_content" in FEATS:
    n_flip = 0
    for ch, rows in D.groupby("chr").groups.items():
        for c in eig:
            r = D.loc[rows, c].corr(D.loc[rows, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0:
                D.loc[rows, c] = -D.loc[rows, c]; n_flip += 1
    log(f"Hi-C eigenvectors oriented by GC: {n_flip} signs flipped (eigenvector x chromosome)")

# ============================ step 2: redundancy clusters (shared by all omics) ============================
# Spearman on all bins of the dataset (pairwise complete), so the clusters are the same for every omic; only the labels
# change between omics.
rho = D[FEATS].corr(method="spearman").abs()
# NOTE (19/11/2026): an undefined |rho| (e.g. at 1 Mb dist.to.closest.gene_protein_coding is constant = 0 on every bin
# where the gene-level features are non-NA) is set to 0 = not redundant, the same rule as pilot_omics_background.py
# (.corr().abs().fillna(0)); with complete linkage such a pair can never share a cluster. Logged pair by pair.
undefined = [(f, g) for i, f in enumerate(FEATS) for g in FEATS[i + 1:] if pd.isna(rho.loc[f, g])]
if undefined:
    log(f"{len(undefined)} feature pairs with undefined Spearman correlation, set to |rho| = 0: "
        + ", ".join(f"{f} ~ {g}" for f, g in undefined))
    rho = rho.fillna(0)
dist = 1 - rho.to_numpy(copy=True); np.fill_diagonal(dist, 0)
tree = linkage(squareform(dist, checks=False), method="complete")
cluster_id = fcluster(tree, t=1 - a.redundancy_rho, criterion="distance")
CL = pd.DataFrame({"feature": FEATS, "cluster": cluster_id})
# medoid: the member with the highest mean |rho| to the OTHER members (scale-free, same space as the tree)
mean_rho = {}
for k, members in CL.groupby("cluster")["feature"]:
    m = list(members)
    for f in m:
        others = [g for g in m if g != f]
        mean_rho[f] = rho.loc[f, others].mean() if others else 1.0
CL["mean_abs_rho_in_cluster"] = CL["feature"].map(mean_rho)
rep = CL.loc[CL.groupby("cluster")["mean_abs_rho_in_cluster"].idxmax(), ["cluster", "feature"]]
CL["representative"] = CL["cluster"].map(dict(zip(rep["cluster"], rep["feature"])))
CL["cluster_size"] = CL["cluster"].map(CL["cluster"].value_counts())
log(f"redundancy clusters (complete linkage, |rho| >= {a.redundancy_rho}): {CL['cluster'].nunique()} clusters from "
    f"{len(FEATS)} features; largest: {CL['cluster_size'].max()} features")

if a.write_label_template:
    T = CL.sort_values(["cluster", "feature"])[["feature", "cluster", "representative"]].copy()
    for om in a.omics: T[om] = ""
    T.to_csv(f"{a.out_dir}/label_template.tsv", sep="\t", index=False)
    log(f"written {a.out_dir}/label_template.tsv: fill the omic columns with O / S / B and pass it with --labels")
    sys.exit(0)

# ============================ step 1: prior labels ============================
if not a.labels:
    sys.exit("--labels is required (or use --write_label_template to create the file to fill)")
L = pd.read_csv(a.labels, sep="\t", dtype=str)
missing_cols = [c for c in ["feature"] + a.omics if c not in L.columns]
if missing_cols:
    sys.exit(f"{a.labels}: columns missing: {missing_cols} (needed: feature + one column per requested omic)")
L = L.set_index("feature")
unknown = sorted(set(L.index) - set(FEATS))
if unknown:
    print(f"NOTE: {len(unknown)} features of {a.labels} are not in the feature set and are ignored: {unknown}")
not_labelled = [f for f in FEATS if f not in L.index]
if not_labelled:
    sys.exit(f"{a.labels}: {len(not_labelled)} features of the feature set have no row: {not_labelled}")
for om in a.omics:
    vals = L.loc[FEATS, om].str.strip().str.upper()
    bad = vals[~vals.isin(["O", "S", "B"])]
    if len(bad):
        sys.exit(f"{a.labels}, column {om}: labels must be O, S or B; wrong or empty for {len(bad)} features: "
                 f"{bad.index.tolist()[:20]}")
    # NOTE (19/11/2026): B (both) is treated as O, as the mixed clusters ("occurrence until proven selection", Gab)
    CL[f"prior_{om}"] = CL["feature"].map(vals)
    n_B = int((vals == "B").sum())
    if n_B:
        log(f"{om}: {n_B} features labelled B, treated as O: {vals.index[vals == 'B'].tolist()}")
    vals = vals.replace("B", "O")
    CL[f"label_{om}"] = CL["feature"].map(vals)
    # cluster label: unanimous label, otherwise O ("occurrence until proven selection")
    per_cluster = CL.groupby("cluster")[f"label_{om}"].agg(lambda s: s.iloc[0] if s.nunique() == 1 else "MIXED")
    mixed = per_cluster[per_cluster == "MIXED"].index
    CL[f"cluster_label_{om}"] = CL["cluster"].map(per_cluster.replace("MIXED", "O"))
    CL[f"mixed_{om}"] = CL["cluster"].isin(mixed)
    if len(mixed):
        parts = []
        for k in mixed:
            members = CL[CL["cluster"] == k]
            parts.append(f"cluster {k}: " + ", ".join(f"{f}={lab}" for f, lab in zip(members["feature"], members[f"label_{om}"])))
        detail = "; ".join(parts)
        log(f"{om}: {len(mixed)} mixed clusters labelled O -> {detail}")

# NOTE (20/11/2026): label randomization (Gab, 20/11/2026), the null of the concordance test. The CLUSTER labels (after
# the mixed -> O rule) are permuted among clusters: the class sizes in clusters stay those of the truth, the clusters
# themselves (label-free) do not change. ONE permutation shared by all omics: it depends only on --label_shuffle_seed and
# on the cluster ids, so separate jobs (one per omic) with the same seed get the same shuffle. The true labels stay in
# clusters.tsv as true_cluster_label_<omic>; label_<omic> / prior_<omic> (feature level) are left as they were.
if a.label_shuffle_seed is not None:
    cluster_ids = np.sort(CL["cluster"].unique())
    permuted_ids = np.random.default_rng(a.label_shuffle_seed).permutation(cluster_ids)
    source_of = dict(zip(cluster_ids, permuted_ids))           # cluster k takes the true label of cluster source_of[k]
    for om in a.omics:
        true_label = CL.groupby("cluster")[f"cluster_label_{om}"].first()
        CL[f"true_cluster_label_{om}"] = CL[f"cluster_label_{om}"]
        CL[f"cluster_label_{om}"] = CL["cluster"].map(lambda k: true_label[source_of[k]])
        per_cluster = CL.groupby("cluster")[[f"true_cluster_label_{om}", f"cluster_label_{om}"]].first()
        n_changed = int((per_cluster.iloc[:, 0] != per_cluster.iloc[:, 1]).sum())
        log(f"{om}: cluster labels shuffled (seed {a.label_shuffle_seed}): {n_changed} of {len(per_cluster)} clusters "
            f"changed class; sizes O={int((per_cluster.iloc[:, 1] == 'O').sum())} S={int((per_cluster.iloc[:, 1] == 'S').sum())} "
            f"(true: O={int((per_cluster.iloc[:, 0] == 'O').sum())} S={int((per_cluster.iloc[:, 0] == 'S').sum())})")
CL.to_csv(f"{a.out_dir}/clusters.tsv", sep="\t", index=False)

# ============================ targets, families, models ============================
TARGETS = {"cna_ampl": ("ampl.frequency", "binomial", None), "cna_del": ("del.frequency", "binomial", None),
           "atac": ("logFC", "gaussian", None), "mut_syn": ("mut_n_syn", "poisson", "sites_syn"),
           "mut_mis": ("mut_n_mis", "poisson", "sites_mis"), "mut_trunc": ("mut_n_trunc", "poisson", "sites_trunc")}
XGB_OBJ = {"gaussian": "reg:squarederror", "poisson": "count:poisson", "binomial": "binary:logistic"}

def deviance(y, mu, family):
    """Deviance (up to multiplicative constants, which cancel in the ratio) - same as pilot_omics_background.py."""
    if family == "gaussian":
        return np.sum((y - mu) ** 2)
    if family == "poisson":
        mu = np.clip(mu, 1e-12, None)
        return 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0) - (y - mu))
    mu = np.clip(mu, 1e-9, 1 - 1e-9)                          # binomial on a fraction of patients
    t1 = np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)
    t0 = np.where(y < 1, (1 - y) * np.log(np.where(y < 1, 1 - y, 1) / (1 - mu)), 0.0)
    return 2 * np.sum(t1 + t0)

def corr_len(y, max_lag=3000):
    """Lag (bins) at which the autocorrelation of y drops below 0.5 - buffer of the CV, as in the pilot."""
    y = y - y.mean(); v = np.dot(y, y)
    if v == 0: return 1
    for lag in range(1, min(max_lag, len(y) - 1)):
        if np.dot(y[:-lag], y[lag:]) / v < 0.5: return lag
    return max_lag

# NOTE (19/11/2026): models = the same XGBoost regressors as pilot_omics_background.py (--loss glm): same objective per
# family and same fixed hyperparameters, 300 trees; NaN handled natively by XGBoost (no imputation, no scaling).
XGB_PARAMS = dict(max_depth=5, eta=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=10, tree_method="hist")

class OmicCV:
    """Folds, base margins and null model of one omic, computed once and reused by every fit."""
    def __init__(self, d, y, family, offset, reps):
        self.y, self.family, self.reps = y, family, reps
        chroms = d["chr"].to_numpy()
        ell = int(np.median([corr_len(y[chroms == c]) for c in np.unique(chroms)]))
        self.buffer = min(ell, a.max_buffer_bins)
        self.folds = INTERLEAVED_KFOLD_IDXS(d["binID"].to_numpy(), a.n_splits, a.block_bins, self.buffer, MIN_TRAIN_FRACTION=0.05)
        self.X = d[reps].to_numpy(float)
        off = np.zeros(len(y)) if offset is None else offset
        self.parts, null = [], np.full(len(y), np.nan)
        for k, (tr, te) in self.folds.items():
            margin = None
            if family == "poisson":                          # start from the mean rate of the training fold (see the pilot)
                margin = off + np.log(max(y[tr].sum(), 1e-12) / np.exp(off[tr]).sum())
                null[te] = np.exp(off[te]) * y[tr].sum() / np.exp(off[tr]).sum()
            else:
                null[te] = y[tr].mean()
            self.parts.append((tr, te, None if margin is None else margin[tr], None if margin is None else margin[te]))
        self.d_null = deviance(y, null, family)

    def dev_expl(self, features, X=None):
        """Held-out deviance explained of XGBoost on these representatives, one value per seed (array aligned with
        a.seeds). Empty list = null model = 0 for every seed.
        NOTE (20/11/2026): X = optional feature matrix replacing self.X (same columns: self.reps), used by the
        circular-shift null of the O+S model; same folds, margins and null model."""
        if not features:
            return np.zeros(len(a.seeds))
        X_used = self.X if X is None else X
        cols = [self.reps.index(f) for f in features]
        out = []
        for seed in a.seeds:
            pred = np.full(len(self.y), np.nan)
            params = dict(XGB_PARAMS, objective=XGB_OBJ[self.family], nthread=a.n_jobs, seed=seed)
            for tr, te, mtr, mte in self.parts:
                dtr = xgb.DMatrix(X_used[np.ix_(tr, cols)], label=self.y[tr], base_margin=mtr)
                dte = xgb.DMatrix(X_used[np.ix_(te, cols)], base_margin=mte)
                booster = xgb.train(params, dtr, num_boost_round=a.n_trees)
                pred[te] = booster.predict(dte)
            out.append(1 - deviance(self.y, pred, self.family) / self.d_null)
        return np.array(out)

# ============================ steps 3-4: base models, removals, swaps ============================
base_rows, removal_rows, swap_rows, addition_rows = [], [], [], []
full_rows = []                                                  # NOTE (20/11/2026): O+S model and its shifted-S null
for om in a.omics:
    col, family, expo = TARGETS[om]
    if col not in D.columns or (expo and expo not in D.columns):
        sys.exit(f"{om}: column {col if col not in D.columns else expo} missing in {a.dataset}")
    d = D
    if family == "poisson":                                   # same bins as the pilot: enough CDS and sites of the class
        if not EXONIC: sys.exit(f"{om}: column cds_bp missing in {a.dataset}")
        n0 = len(d); d = d[(d["cds_bp"] >= a.min_exonic_bp) & (d[expo] > 0)]
        log(f"{om}: {n0 - len(d)} bins removed (cds_bp < {a.min_exonic_bp:g} or no {expo})")
    n0 = len(d); d = d[d[col].notna()]
    if len(d) < n0: log(f"{om}: {n0 - len(d)} bins removed (NA in {col})")
    y = d[col].to_numpy(float)
    if family == "binomial" and (y.min() < 0 or y.max() > 1):
        sys.exit(f"{om}: {col} must be a fraction in [0, 1] for the binomial family (range {y.min()} .. {y.max()})")
    offset = np.log(d[expo].to_numpy(float)) if family == "poisson" else None

    reps_by_class = {}
    for cls in ("O", "S"):
        reps_by_class[cls] = sorted(CL.loc[CL[f"cluster_label_{om}"] == cls, "representative"].unique())
    if not reps_by_class["O"] or not reps_by_class["S"]:
        log(f"{om}: SKIPPED, one class has no cluster (O: {len(reps_by_class['O'])}, S: {len(reps_by_class['S'])})")
        continue
    reps = reps_by_class["O"] + reps_by_class["S"]
    nO, nS = len(reps_by_class["O"]), len(reps_by_class["S"])
    cv = OmicCV(d, y, family, offset, reps)
    log(f"{om} [{family}]: {len(y)} bins, buffer {cv.buffer} bins, clusters O={nO} S={nS}, "
        f"fits to run: {(0 if a.full_models_only else 2 * nO * nS + 2 * (nO + nS)) + 2 + 1 + a.n_shift} "
        f"(x {a.n_splits} folds x {len(a.seeds)} seeds)")
    cluster_of = dict(zip(CL["feature"], CL["cluster"]))

    dev_base_of = {}                                            # per class, one value per seed (reused by the O+S block)
    for cls, other in (("O", "S"), ("S", "O")):
        base = reps_by_class[cls]
        dev_base = cv.dev_expl(base)                             # one value per seed
        dev_base_of[cls] = dev_base
        for i, seed in enumerate(a.seeds):
            base_rows.append(dict(omic=om, model_class=cls, seed=seed, n_clusters=len(base), n_bins=len(y),
                                  buffer_bins=cv.buffer, dev_expl=dev_base[i]))
        log(f"  {om} {cls}-only model: deviance explained {dev_base.mean():.4f} (sd over seeds {dev_base.std(ddof=1) if len(a.seeds) > 1 else 0:.4f}; "
            f"{len(base)} clusters)")
        if a.full_models_only:                                   # NOTE (20/11/2026): base models only, no swap test
            continue
        # NOTE (19/11/2026): every cluster B of the other class added to the FULL model (nothing removed): its own added
        # value, the reference for d_add(B | A removed) -> shared(A, B)
        dev_full_plus = {}
        for B in reps_by_class[other]:
            dev_full_plus[B] = cv.dev_expl(base + [B])
            for i, seed in enumerate(a.seeds):
                addition_rows.append(dict(omic=om, model_class=cls, in_cluster=cluster_of[B], in_rep=B, seed=seed,
                                          dev_base=dev_base[i], dev_full_plus_B=dev_full_plus[B][i],
                                          d_add_full=dev_full_plus[B][i] - dev_base[i]))
        for A in base:
            without_A = [f for f in base if f != A]
            dev_removed = cv.dev_expl(without_A)
            for i, seed in enumerate(a.seeds):
                removal_rows.append(dict(omic=om, model_class=cls, out_cluster=cluster_of[A], out_rep=A, seed=seed,
                                         dev_base=dev_base[i], dev_removed=dev_removed[i], d_remove=dev_removed[i] - dev_base[i]))
            for B in reps_by_class[other]:
                dev_swapped = cv.dev_expl(without_A + [B])
                for i, seed in enumerate(a.seeds):             # paired deltas: same seed for base, removal and swap
                    swap_rows.append(dict(omic=om, model_class=cls, out_cluster=cluster_of[A], out_rep=A,
                                          in_cluster=cluster_of[B], in_rep=B, seed=seed, dev_base=dev_base[i],
                                          dev_removed=dev_removed[i], dev_swapped=dev_swapped[i],
                                          d_remove=dev_removed[i] - dev_base[i], d_add=dev_swapped[i] - dev_removed[i],
                                          d_swap=dev_swapped[i] - dev_base[i],
                                          d_add_full=dev_full_plus[B][i] - dev_base[i],
                                          shared=(dev_swapped[i] - dev_removed[i]) - (dev_full_plus[B][i] - dev_base[i])))
        log(f"  {om} {cls}-only model: {len(reps_by_class[other])} additions to the full model, {len(base)} removals and "
            f"{len(base) * len(reps_by_class[other])} swaps done")

    # NOTE (20/11/2026): O+S full model (Gab, 20/11/2026): does S add to O only in combination (S clusters together, or
    # S with O) although single S clusters add ~nothing? Null: the same O+S model with the S columns circularly shifted
    # TOGETHER (one offset per chromosome shared by all S columns, 10-90% of the chromosome, as the pilot's floor and
    # decoys): keeps the values, the correlations inside S, the spatial autocorrelation and the number of columns, and
    # breaks only the alignment of S with the target and with O. Within-chromosome shift: what S explains through
    # chromosome identity alone stays in the null. Same folds, margins and seeds as every other fit (paired).
    all_reps = reps_by_class["O"] + reps_by_class["S"]
    dev_OS = cv.dev_expl(all_reps)
    for i, seed in enumerate(a.seeds):
        full_rows.append(dict(omic=om, model="O+S", shift_rep=-1, seed=seed, dev_O=dev_base_of["O"][i],
                              dev_S=dev_base_of["S"][i], dev_expl=dev_OS[i], gain_over_O=dev_OS[i] - dev_base_of["O"][i]))
    s_cols = [cv.reps.index(f) for f in reps_by_class["S"]]
    shift_rng = np.random.default_rng(a.shift_seed)            # per omic: the same shifts whatever --omics lists
    chroms_om = d["chr"].to_numpy()
    dev_shifted_all = []
    for r in range(a.n_shift):
        idx = np.arange(len(chroms_om))
        for c in np.unique(chroms_om):                         # bins are sorted by position inside each chromosome
            m = np.where(chroms_om == c)[0]; n = len(m)
            idx[m] = np.roll(m, int(shift_rng.integers(int(.1 * n), max(int(.1 * n) + 1, int(.9 * n)))))
        X_shifted = cv.X.copy()
        X_shifted[:, s_cols] = cv.X[idx][:, s_cols]            # the whole S block moves together; O untouched
        dev_shifted = cv.dev_expl(all_reps, X=X_shifted)
        dev_shifted_all.append(dev_shifted)
        for i, seed in enumerate(a.seeds):
            full_rows.append(dict(omic=om, model="O+S_shifted", shift_rep=r, seed=seed, dev_O=dev_base_of["O"][i],
                                  dev_S=dev_base_of["S"][i], dev_expl=dev_shifted[i],
                                  gain_over_O=dev_shifted[i] - dev_base_of["O"][i]))
    msg = (f"  {om} O+S model: deviance explained {dev_OS.mean():.4f} (gain over O-only "
           f"{(dev_OS - dev_base_of['O']).mean():+.4f})")
    if a.n_shift:
        null_gain = np.array(dev_shifted_all).mean(axis=1) - dev_base_of["O"].mean()   # one value per shift (mean over seeds)
        n_ge = int((null_gain >= (dev_OS - dev_base_of["O"]).mean()).sum())
        msg += (f"; S shifted ({a.n_shift} shifts): gain {null_gain.mean():+.4f} (range {null_gain.min():+.4f} .. "
                f"{null_gain.max():+.4f}); shifts with gain >= observed: {n_ge}")
    log(msg)

pd.DataFrame(base_rows).to_csv(f"{a.out_dir}/base_models.tsv", sep="\t", index=False)
pd.DataFrame(full_rows).to_csv(f"{a.out_dir}/full_models.tsv", sep="\t", index=False)
if a.full_models_only:                                          # NOTE (20/11/2026): no swap tables (none computed)
    log(f"--full_models_only: written base_models.tsv, full_models.tsv, clusters.tsv in {a.out_dir} (no swap test)")
    sys.exit(0)
pd.DataFrame(removal_rows).to_csv(f"{a.out_dir}/removal_deltas.tsv", sep="\t", index=False)
pd.DataFrame(addition_rows).to_csv(f"{a.out_dir}/addition_deltas.tsv", sep="\t", index=False)
SW = pd.DataFrame(swap_rows); SW.to_csv(f"{a.out_dir}/swap_deltas.tsv", sep="\t", index=False)
# summary over seeds: mean and sd of every delta (sd = how much of a delta is XGBoost randomness)
keys = ["omic", "model_class", "out_cluster", "out_rep", "in_cluster", "in_rep"]
SUM = SW.groupby(keys)[["d_remove", "d_add", "d_swap", "d_add_full", "shared"]].agg(["mean", "std"])
SUM.columns = [f"{c}_{s}" for c, s in SUM.columns]
SUM = SUM.reset_index(); SUM["n_seeds"] = len(a.seeds)
SUM.to_csv(f"{a.out_dir}/swap_deltas_summary.tsv", sep="\t", index=False)
log(f"written base_models.tsv, removal_deltas.tsv, addition_deltas.tsv, swap_deltas.tsv (per seed), swap_deltas_summary.tsv (mean/sd over "
    f"{len(a.seeds)} seeds), clusters.tsv in {a.out_dir}")
