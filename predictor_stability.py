#!/usr/bin/env python
"""
PILOT - stability of the predictors across contexts (Gab, 23/11/2026).

Assumptions tested (Gab, 23/11/2026):
  - O predictors are always relevant: their correlation with the target varies little across contexts (low sd);
  - S predictors vary across contexts (high sd);
  - non-important predictors are always less important than average (low mean |rho|).
Contexts = cohort x chromosome, pooled, one panel per resolution. Measure = Spearman correlation only (no model).

For every dataset (cohort, resolution), mutation class and chromosome: Spearman rho between each feature and the
mutation RATE of the class (mut_n_<class> / sites_<class>; raw counts would mostly follow the number of sites), on the
same bins as the pilot (cds_bp >= --min_exonic_bp and sites of the class > 0). Then, per resolution, omic and feature,
over all (cohort, chromosome) values pooled:
  mean_rho, mean_abs_rho, sd_rho (signed: a sign flip counts as variability), sd_abs_rho,
  norm_sd = sd_rho / mean_abs_rho (sd relative to the average strength; large when mean_abs_rho is close to 0).
Noise reference: the same statistics with the target circularly shifted inside each chromosome (--n_shift repeats, shift
between 10% and 90% of the chromosome, as the pilot floor). It keeps the autocorrelation of target and features and
breaks their alignment, so its sd is the variability expected from chance alone with the same number of bins per
chromosome (small chromosomes give noisier rho).

Feature set: the same rules as pilot_predictor_identity.py (no alteration / exposure / cds_bp / OG-TSG distance columns;
distance NAs capped at 2 Mb; optional --gene_na_to_zero; features with >= 50% NA or constant dropped, per dataset;
Hi-C eigenvectors oriented by GC per chromosome). rho uses the bins where the feature is not NA; a chromosome with fewer
than --min_bins usable bins, or where the feature is constant, gives no value for that feature (logged).
Labels: --labels (config/predictor_identity_labels.tsv), column of each omic; O, S and B shown separately.

Inputs
  --dataset COHORT BIN_SIZE PATH   repeated, one per feature matrix (e.g. BRCA 1000000 data/feature_matrix/<name>)
  --labels                         feature + one column per omic (O / S / B)
Outputs (in --out_dir)
  correlations.tsv          resolution, cohort, chrom, omic, feature, n_bins, rho
  shifted_correlations.tsv  the same with shift_rep (circularly shifted target)
  stability.tsv             per resolution, omic, feature: label, n_values, mean_rho, mean_abs_rho, sd_rho, sd_abs_rho,
                            norm_sd, and the noise reference (mean over shifts of the same statistics, *_shift)
  stability_shifted.tsv     the same statistics per shift repeat (the grey reference of the figures)
  predictor_stability_<omic>.pdf/.png  one row per resolution: left, mean |rho| vs sd(rho) per feature, coloured by
                            label (grey: shifted target); right, distribution of norm_sd per label (log scale)
Usage (cluster, pilot env, from the root of the repository)
  python predictor_stability.py --labels config/predictor_identity_labels.tsv --gene_na_to_zero --out_dir analyses/v4/predictor_stability \
      --dataset BRCA 1000000 data/feature_matrix/<BRCA 1 Mb> --dataset LUAD 1000000 data/feature_matrix/<LUAD 1 Mb> \
      --dataset BRCA 20000 data/feature_matrix/<BRCA 20 kb> ...
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import rankdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--dataset", nargs=3, action="append", metavar=("COHORT", "BIN_SIZE", "PATH"), required=True)
p.add_argument("--labels", required=True)
p.add_argument("--omics", nargs="+", default=["mut_syn", "mut_mis", "mut_trunc"])
p.add_argument("--min_exonic_bp", type=float, default=500, help="bins with at least this many coding bp (as in the pilot)")
p.add_argument("--min_bins", type=int, default=20, help="minimum usable bins of a chromosome for a rho value")
p.add_argument("--min_values_frac", type=float, default=0.5,
               help="figures only: features with fewer (cohort x chromosome) values than this fraction of the maximum at "
                    "that resolution are left out (their sd is unreliable); they stay in the tables")
p.add_argument("--n_shift", type=int, default=20, help="circular shifts of the target (noise reference)")
p.add_argument("--shift_seed", type=int, default=1)
p.add_argument("--gene_na_to_zero", action="store_true",
               help="NA of the gene-level features set to 0 before the >= 50%% NA filter (as in pilot_predictor_identity.py)")
p.add_argument("--out_dir", default="predictor_stability")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
t0 = time.time()
def log(msg): print(f"[{time.time() - t0:.0f}s] {msg}", flush=True)

TARGETS = {"mut_syn": ("mut_n_syn", "sites_syn"), "mut_mis": ("mut_n_mis", "sites_mis"),
           "mut_trunc": ("mut_n_trunc", "sites_trunc")}
unknown_omics = [om for om in a.omics if om not in TARGETS]
if unknown_omics:
    sys.exit(f"--omics: only mutation classes are supported ({list(TARGETS)}); got {unknown_omics}")
L = pd.read_csv(a.labels, sep="\t", dtype=str).set_index("feature")
missing_cols = [om for om in a.omics if om not in L.columns]
if missing_cols:
    sys.exit(f"{a.labels}: no label column for {missing_cols}")
chrom_order = lambda v: int(v) if str(v).isdigit() else 100 + ord(str(v)[0])
res_name = lambda b: f"{b // 1000000} Mb" if b >= 1000000 else f"{b // 1000} kb"


def feature_set(D):
    """Same rules as pilot_predictor_identity.py; returns the feature list (D is modified in place: NA fills, Hi-C signs)."""
    alt = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
           or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
    validation_only = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]
    exposure = [c for c in D.columns if c.startswith(("sites_", "wsites_"))]
    feats = [c for c in D.columns if c not in alt + validation_only + exposure + ["cds_bp", "binID", "chr", "idx"]]
    if a.gene_na_to_zero:
        gene_level = [c for c in feats
                      if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                               "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
                      or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs",
                                       "total_n_paralogs"))]
        for c in gene_level:
            D[c] = D[c].fillna(0)
        log(f"  --gene_na_to_zero: NA set to 0 in {len(gene_level)} gene-level features")
    for c in [c for c in feats if c.startswith("dist")]:
        D[c] = D[c].fillna(2e6)                                 # NA on a whole chromosome = no gene of that list: cap
    too_many_na = [c for c in feats if D[c].isna().mean() >= 0.5]
    constant = [c for c in feats if c not in too_many_na and D[c].nunique() <= 1]
    feats = [c for c in feats if c not in too_many_na + constant]
    log(f"  {len(feats)} features; dropped {len(too_many_na)} with >= 50% NA {too_many_na}, {len(constant)} constant {constant}")
    # Hi-C eigenvectors: arbitrary sign per chromosome -> positive correlation with GC (as the pilot); without this the
    # sign of rho would flip at random between chromosomes and inflate sd_rho
    eig = [c for c in feats if c.startswith("eigen_")]
    if eig and "GC_content" in feats:
        for ch, rows in D.groupby("chr").groups.items():
            for c in eig:
                r = D.loc[rows, c].corr(D.loc[rows, "GC_content"], method="spearman")
                if pd.notna(r) and r < 0:
                    D.loc[rows, c] = -D.loc[rows, c]
    return feats


def spearman(x, y):
    """Spearman rho on the bins where x is not NA (y has no NA); NaN if too few bins or x / y constant there."""
    ok = ~np.isnan(x)
    if ok.sum() < a.min_bins:
        return np.nan, int(ok.sum())
    rx, ry = rankdata(x[ok]), rankdata(y[ok])
    if rx.std() == 0 or ry.std() == 0:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(rx, ry)[0, 1]), int(ok.sum())


# ============================ correlations per dataset, omic, chromosome ============================
corr_rows, shift_rows = [], []
for cohort, bin_size, path in a.dataset:
    bin_size = int(bin_size)
    log(f"{cohort} {res_name(bin_size)}: {path}")
    D = pd.read_csv(path, sep="\t")
    D["chr"] = D["binID"].str.split("_").str[0]
    D["idx"] = D["binID"].str.split("_").str[1].astype(int)
    D = D.sort_values(["chr", "idx"], key=lambda s: s.map(chrom_order) if s.name == "chr" else s).reset_index(drop=True)
    feats = feature_set(D)
    not_labelled = [f for f in feats if f not in L.index]
    if not_labelled:
        sys.exit(f"{a.labels}: {len(not_labelled)} features of {path} have no row: {not_labelled}")
    for om in a.omics:
        count_col, sites_col = TARGETS[om]
        for c in (count_col, sites_col, "cds_bp"):
            if c not in D.columns:
                sys.exit(f"{om}: column {c} missing in {path}")
        keep = (D["cds_bp"] >= a.min_exonic_bp) & (D[sites_col] > 0) & D[count_col].notna()
        d = D[keep]
        log(f"  {om}: {len(D) - len(d)} of {len(D)} bins removed (cds_bp < {a.min_exonic_bp:g}, no {sites_col} or NA count)")
        rng = np.random.default_rng(a.shift_seed)               # same shifts for every dataset and omic with the same bins
        skipped_chroms = []
        for ch, dc in d.groupby("chr", sort=False):
            n = len(dc)
            if n < a.min_bins:
                skipped_chroms.append(f"{ch} ({n})"); continue
            rate = (dc[count_col] / dc[sites_col]).to_numpy(float)
            shifts = [int(rng.integers(int(.1 * n), max(int(.1 * n) + 1, int(.9 * n)))) for _ in range(a.n_shift)]
            for f in feats:
                x = dc[f].to_numpy(float)
                rho, n_used = spearman(x, rate)
                corr_rows.append(dict(resolution=bin_size, cohort=cohort, chrom=ch, omic=om, feature=f, n_bins=n_used, rho=rho))
                for r, k in enumerate(shifts):
                    rho_s, _ = spearman(x, np.roll(rate, k))
                    shift_rows.append(dict(resolution=bin_size, cohort=cohort, chrom=ch, omic=om, feature=f,
                                           shift_rep=r, rho=rho_s))
        if skipped_chroms:
            log(f"  {om}: chromosomes skipped (< {a.min_bins} usable bins): {skipped_chroms}")
    del D

C = pd.DataFrame(corr_rows)
CS = pd.DataFrame(shift_rows)
C.to_csv(f"{a.out_dir}/correlations.tsv", sep="\t", index=False)
CS.to_csv(f"{a.out_dir}/shifted_correlations.tsv", sep="\t", index=False)

# report missing values: a feature without rho on some chromosome (constant there, or too many NA) has fewer values
na_rho = C[C["rho"].isna()].groupby(["resolution", "cohort", "omic", "feature"]).size()
if len(na_rho):
    log(f"rho undefined (feature constant on the chromosome or < {a.min_bins} non-NA bins) for "
        f"{int(na_rho.sum())} (feature x chromosome) values; per feature:\n{na_rho.to_string()}")


# ============================ stability per resolution, omic, feature ============================
def summarise(values):
    v = values.dropna().to_numpy()
    mean_abs = np.abs(v).mean() if len(v) else np.nan
    sd = v.std(ddof=1) if len(v) > 1 else np.nan
    return pd.Series(dict(n_values=len(v), mean_rho=v.mean() if len(v) else np.nan, mean_abs_rho=mean_abs, sd_rho=sd,
                          sd_abs_rho=np.abs(v).std(ddof=1) if len(v) > 1 else np.nan,
                          norm_sd=sd / mean_abs if mean_abs > 0 else np.nan))

keys = ["resolution", "omic", "feature"]
S = C.groupby(keys)["rho"].apply(summarise).unstack().reset_index()
SS = CS.groupby(keys + ["shift_rep"])["rho"].apply(summarise).unstack().reset_index()
shift_mean = SS.groupby(keys)[["mean_abs_rho", "sd_rho", "norm_sd"]].mean().add_suffix("_shift").reset_index()
S = S.merge(shift_mean, on=keys, how="left")
S["label"] = [L.loc[f, om].strip().upper() for f, om in zip(S["feature"], S["omic"])]
SS["label"] = [L.loc[f, om].strip().upper() for f, om in zip(SS["feature"], SS["omic"])]
n_max = S.groupby(["resolution", "omic"])["n_values"].transform("max")
fewer = S[S["n_values"] < n_max]
if len(fewer):
    log(f"features with fewer (cohort x chromosome) values than the others at the same resolution (absent from a "
        f"dataset or rho undefined on some chromosomes):\n"
        f"{fewer[keys + ['n_values']].assign(of=n_max[fewer.index]).to_string(index=False)}")
S.to_csv(f"{a.out_dir}/stability.tsv", sep="\t", index=False)
SS.to_csv(f"{a.out_dir}/stability_shifted.tsv", sep="\t", index=False)

# ============================ figures: one per omic, one row per resolution ============================
LABEL_COLOR = {"O": "#2a78d6", "S": "#eb6834", "B": "#1baf7a"}
SHIFT_COLOR = "#a0a0a0"
OMIC_NAME = {"mut_syn": "Synonymous mutations", "mut_mis": "Missense mutations", "mut_trunc": "Truncating mutations"}
resolutions = sorted(S["resolution"].unique())
# features with few values (e.g. constant on most chromosomes, like a distance to a gene type present on one chromosome)
# have an unreliable sd: left out of the figures, listed here, kept in stability.tsv
few = S["n_values"] < a.min_values_frac * n_max
if few.any():
    log(f"left out of the figures ({int(few.sum())} feature x omic x resolution, < {a.min_values_frac:g} of the maximum "
        f"number of values):\n{S.loc[few, keys + ['n_values']].to_string(index=False)}")
plot_keys = S.loc[~few, keys]
S_plot = S[~few]
SS_plot = SS.merge(plot_keys, on=keys, how="inner")
for om in a.omics:
    fig, axes = plt.subplots(len(resolutions), 2, figsize=(11, 3.6 * len(resolutions)), squeeze=False)
    for i, res in enumerate(resolutions):
        s = S_plot[(S_plot["omic"] == om) & (S_plot["resolution"] == res)]
        ss = SS_plot[(SS_plot["omic"] == om) & (SS_plot["resolution"] == res)]
        n_ctx = int(s["n_values"].max()) if len(s) else 0
        ax = axes[i, 0]
        ax.scatter(ss["mean_abs_rho"], ss["sd_rho"], s=6, color=SHIFT_COLOR, alpha=.25, linewidths=0,
                   label=f"shifted target ({a.n_shift} shifts)")
        for lab, col in LABEL_COLOR.items():
            sl = s[s["label"] == lab]
            if len(sl):
                ax.scatter(sl["mean_abs_rho"], sl["sd_rho"], s=22, color=col, edgecolors="white", linewidths=.5,
                           label=f"{lab} (n={len(sl)})")
        ax.set_xlabel("mean |Spearman rho| over cohort x chromosome")
        ax.set_ylabel("sd of rho")
        ax.set_title(f"{res_name(res)}: {n_ctx} cohort x chromosome values per feature", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax.grid(alpha=.2)
        ax = axes[i, 1]
        all_norm = np.concatenate([s["norm_sd"].dropna().to_numpy(), ss["norm_sd"].dropna().to_numpy()])
        all_norm = all_norm[all_norm > 0]
        if len(all_norm):
            bins = np.logspace(np.log10(all_norm.min()), np.log10(all_norm.max()), 30)
            ax.hist(ss["norm_sd"].dropna(), bins=bins, density=True, color=SHIFT_COLOR, alpha=.5, label="shifted target")
            for lab, col in LABEL_COLOR.items():
                v = s.loc[s["label"] == lab, "norm_sd"].dropna()
                if len(v):
                    ax.hist(v, bins=bins, density=True, histtype="step", linewidth=2, color=col,
                            label=f"{lab} (median {v.median():.2f})")
            ax.set_xscale("log")
        ax.set_xlabel("sd(rho) / mean |rho|  (log scale)")
        ax.set_ylabel("density")
        ax.legend(fontsize=8, frameon=False)
        ax.grid(alpha=.2)
    fig.suptitle(f"{OMIC_NAME[om]}: stability of the predictor correlations (rate = count / sites)", fontsize=12)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out_dir}/predictor_stability_{om}.{ext}", dpi=150)
    plt.close(fig)
    log(f"written {a.out_dir}/predictor_stability_{om}.pdf/.png")
log(f"written correlations.tsv, shifted_correlations.tsv, stability.tsv, stability_shifted.tsv in {a.out_dir}")
