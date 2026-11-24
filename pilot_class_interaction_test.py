"""
PILOT - per-feature supply vs selection test on mutation classes (a "dN/dS slope" per feature).
Gab, 24/11/2026.

Idea: supply (mutation rate, repair, ...) moves synonymous and non-synonymous counts TOGETHER; selection moves only the
non-synonymous ones. For every feature (one at a time, no background model):
  * supply effect beta  : Poisson GLM   syn ~ offset(log sites_syn) + a + beta * z
  * selection effect gamma (per contrast class C = mis or trunc), conditional on the bin total:
        C | (C + syn) ~ Binomial(C + syn, p),   logit p = log(sites_C / sites_syn) + a + gamma * z
    Any per-bin supply multiplier common to both classes cancels in this ratio, so gamma does NOT depend on how well a
    background is modelled. gamma = change of log dN/dS per SD of the feature.
  z = rank-based normal score of the feature (rank -> normal quantile): effects are per SD, comparable across features
  and robust to skewed features (distances, counts).

Inference: moving-block bootstrap over CONTIGUOUS bins within chromosomes (bins are autocorrelated: bin-level resampling
would give CIs that are too narrow). Reported: estimate, bootstrap SE, 95% percentile CI, p from z = estimate / SE_boot
(two-sided normal), BH q-value per contrast over the tested features (controls excluded).

Controls:
  * positive: dist.to.closest.OG / dist.to.closest.TSG (validation-only columns, never model features). Expected
    gamma < 0 (closer to a driver -> more non-synonymous). Optimistic: driver lists were built partly from TCGA.
  * negative (placebo): synonymous counts split at random in two halves (Binomial(syn, 0.5), equal sites) and tested
    with the SAME procedure (gamma of half A vs half B). Every feature should be null: the fraction with p < 0.05 over
    --n_placebo splits is the false-positive rate of the procedure. LIMIT: both halves share the synonymous
    trinucleotide composition, so the placebo does NOT catch a spectrum confound (a supply feature that changes the
    mutation spectrum, e.g. replication timing / GC, can give gamma != 0 without selection).

Same preprocessing as pilot_predictor_identity.py --gene_na_to_zero (gene-level NA -> 0, distance NA -> 2 Mb, drop >= 50%
NA and constant features, Hi-C eigenvectors oriented by GC per chromosome). Bins with NA in a feature are dropped for that
feature only (logged); bins with 0 sites of a class, or with C + syn = 0 (no information in the binomial), are dropped
per contrast (logged).

Inputs : --dataset  feature matrix (TSV, binID + features + mut_n_* + sites_*)
         --labels   config/predictor_identity_labels.tsv (columns feature, group, mut_syn, mut_mis, mut_trunc), for annotation
Outputs (--out_dir):
  effects.tsv          feature, role (feature / positive_control), prior labels, effect (supply / selection_mis /
                       selection_trunc), estimate, se_boot, ci_low, ci_high, p, q, n_bins, n_dropped_na
  placebo.tsv          feature, split, estimate, se_boot, p   (one row per placebo split)
  placebo_summary.tsv  feature, n_placebo, fpr_05 (fraction of splits with p < 0.05), mean and SD of the placebo gamma
Usage:
  python pilot_class_interaction_test.py --dataset <feature_matrix> --labels config/predictor_identity_labels.tsv \
      --out_dir analyses/v4/LUAD/1000000bp/class_interaction --block_bins 10 --n_boot 1000 --n_placebo 100 --n_boot_placebo 200
"""
import argparse, os, time
import numpy as np, pandas as pd
from scipy.stats import norm

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--labels", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--block_bins", type=int, default=10, help="bins per bootstrap block (contiguous, within a chromosome)")
p.add_argument("--n_boot", type=int, default=1000)
p.add_argument("--n_placebo", type=int, default=100)
p.add_argument("--n_boot_placebo", type=int, default=200)
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
t0 = time.time()
def log(msg): print(f"[{time.time() - t0:.0f}s] {msg}", flush=True)
rng = np.random.default_rng(a.seed)

# ============================ data and feature set (same rules as pilot_predictor_identity.py) ============================
D = pd.read_csv(a.dataset, sep="\t")
D["chr"] = D["binID"].str.split("_").str[0]
D["idx"] = D["binID"].str.split("_").str[1].astype(int)
chrom_order = lambda v: int(v) if str(v).isdigit() else 100 + ord(str(v)[0])
D = D.sort_values(["chr", "idx"], key=lambda s: s.map(chrom_order) if s.name == "chr" else s).reset_index(drop=True)

needed = ["mut_n_syn", "mut_n_mis", "mut_n_trunc", "sites_syn", "sites_mis", "sites_trunc"]
missing = [c for c in needed if c not in D.columns]
if missing:
    raise SystemExit(f"ERROR: {a.dataset} lacks the mutation columns {missing}: use a feature matrix with COSMIC classes.")

ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
CONTROLS = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]
if len(CONTROLS) < 2:
    log(f"WARNING: positive controls found: {CONTROLS} (expected dist.to.closest.OG and dist.to.closest.TSG)")
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_"))]
EXONIC = ["cds_bp"] if "cds_bp" in D.columns else []
FEATS = [c for c in D.columns if c not in ALT + CONTROLS + EXPOSURE + EXONIC + ["binID", "chr", "idx"]]
GENE_LEVEL_NA_TO_ZERO = [c for c in FEATS
                         if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                                  "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
                         or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs",
                                          "total_n_paralogs"))]
for c in GENE_LEVEL_NA_TO_ZERO:
    D[c] = D[c].fillna(0)
for c in [c for c in FEATS + CONTROLS if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)                                   # NA on a whole chromosome = no gene of that list: cap
too_many_na = [c for c in FEATS if D[c].isna().mean() >= 0.5]
constant = [c for c in FEATS if c not in too_many_na and D[c].nunique() <= 1]
FEATS = [c for c in FEATS if c not in too_many_na + constant]
log(f"{len(D)} bins, {len(FEATS)} features + {len(CONTROLS)} positive controls; dropped: {len(too_many_na)} with >= 50% "
    f"NA {too_many_na}, {len(constant)} constant {constant}")
eig = [c for c in FEATS if c.startswith("eigen_")]
if eig and "GC_content" in FEATS:
    for ch, rows in D.groupby("chr").groups.items():
        for c in eig:
            r = D.loc[rows, c].corr(D.loc[rows, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0:
                D.loc[rows, c] = -D.loc[rows, c]

L = pd.read_csv(a.labels, sep="\t").set_index("feature")
not_labelled = [f for f in FEATS if f not in L.index]
if not_labelled:
    raise SystemExit(f"ERROR: {a.labels} has no row for {len(not_labelled)} features: {not_labelled}")

# ============================ GLM with one covariate (intercept + slope, offset), Newton-Raphson ============================
def fit_glm_slope(z, offset, family, y, n=None, max_iter=50, tol=1e-10):
    """Slope of a GLM with intercept, one covariate z and an offset. family 'poisson' (y counts) or 'binomial'
    (y successes out of n trials, logit link). Returns the slope (np.nan if Newton does not converge)."""
    X = np.column_stack([np.ones_like(z), z])
    # start from the intercept-only solution (z slope 0)
    if family == "poisson":
        b = np.array([np.log(y.sum() / np.exp(offset).sum()), 0.0])
    else:
        frac = y.sum() / n.sum()
        b = np.array([np.log(frac / (1 - frac)) - np.mean(offset), 0.0])
    for _ in range(max_iter):
        eta = offset + X @ b
        if family == "poisson":
            mu = np.exp(eta); grad_w = y - mu; w = mu
        else:
            prob = 1 / (1 + np.exp(-eta)); grad_w = y - n * prob; w = n * prob * (1 - prob)
        score = X.T @ grad_w
        info = X.T @ (X * w[:, None])
        step = np.linalg.solve(info, score)
        b = b + step
        if np.max(np.abs(step)) < tol:
            return b[1]
    return np.nan

# ============================ moving-block bootstrap (blocks within chromosomes) ============================
def make_blocks(chrom, block_bins):
    """List of arrays of row positions: consecutive bins of the same chromosome, block_bins per block (last one shorter)."""
    blocks = []
    for ch in pd.unique(chrom):
        pos = np.flatnonzero(chrom == ch)
        for start in range(0, len(pos), block_bins):
            blocks.append(pos[start:start + block_bins])
    return blocks

def bootstrap_slope(fit_fun, blocks, n_boot, rng):
    """fit_fun(rows) -> slope on those rows. Resample blocks with replacement (as many as there are). Returns array."""
    out = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, len(blocks), len(blocks))
        rows = np.concatenate([blocks[j] for j in pick])
        out[i] = fit_fun(rows)
    return out

def summarise(est, boot):
    """estimate, SE, 95% percentile CI, two-sided p from est / SE (normal). NaN replicates (no convergence) are logged."""
    good = boot[np.isfinite(boot)]
    n_bad = len(boot) - len(good)
    se = good.std(ddof=1)
    pval = 2 * norm.sf(abs(est) / se) if se > 0 else np.nan
    return dict(estimate=est, se_boot=se, ci_low=np.quantile(good, 0.025), ci_high=np.quantile(good, 0.975),
                p=pval, n_boot_failed=n_bad)

def rank_normal(x):
    """Rank-based normal score (average ranks for ties)."""
    r = pd.Series(x).rank(method="average").to_numpy()
    return norm.ppf((r - 0.5) / len(r))

def bh(pvals):
    p_arr = np.asarray(pvals, float); m = len(p_arr)
    order = np.argsort(p_arr)
    q = p_arr[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m); out[order] = np.minimum(q, 1)
    return out

# ============================ tests ============================
CONTRASTS = {"mis": ("mut_n_mis", "sites_mis"), "trunc": ("mut_n_trunc", "sites_trunc")}
rows_out, placebo_rows = [], []
for feat in FEATS + CONTROLS:
    role = "positive_control" if feat in CONTROLS else "feature"
    ok_feat = D[feat].notna().to_numpy()
    n_na = int((~ok_feat).sum())
    if n_na:
        log(f"{feat}: {n_na} bins with NA dropped for this feature")

    # ---- supply: Poisson on synonymous ----
    keep = ok_feat & (D["sites_syn"].to_numpy() > 0)
    sub = D.loc[keep]
    z = rank_normal(sub[feat].to_numpy())
    y = sub["mut_n_syn"].to_numpy(float); off = np.log(sub["sites_syn"].to_numpy(float))
    blocks = make_blocks(sub["chr"].to_numpy(), a.block_bins)
    est = fit_glm_slope(z, off, "poisson", y)
    boot = bootstrap_slope(lambda r: fit_glm_slope(z[r], off[r], "poisson", y[r]), blocks, a.n_boot, rng)
    rows_out.append(dict(feature=feat, role=role, effect="supply", n_bins=int(keep.sum()), n_dropped_na=n_na,
                         **summarise(est, boot)))

    # ---- selection: binomial of class C vs synonymous ----
    for cname, (ycol, scol) in CONTRASTS.items():
        keep = ok_feat & (D["sites_syn"].to_numpy() > 0) & (D[scol].to_numpy() > 0) \
               & ((D[ycol] + D["mut_n_syn"]).to_numpy() > 0)
        sub = D.loc[keep]
        z = rank_normal(sub[feat].to_numpy())
        yc = sub[ycol].to_numpy(float); ntot = yc + sub["mut_n_syn"].to_numpy(float)
        off = np.log(sub[scol].to_numpy(float) / sub["sites_syn"].to_numpy(float))
        blocks = make_blocks(sub["chr"].to_numpy(), a.block_bins)
        est = fit_glm_slope(z, off, "binomial", yc, ntot)
        boot = bootstrap_slope(lambda r: fit_glm_slope(z[r], off[r], "binomial", yc[r], ntot[r]), blocks, a.n_boot, rng)
        rows_out.append(dict(feature=feat, role=role, effect=f"selection_{cname}", n_bins=int(keep.sum()),
                             n_dropped_na=n_na, **summarise(est, boot)))

    # ---- placebo: random half of synonymous vs the other half (same sites -> offset 0) ----
    keep = ok_feat & ((D["mut_n_syn"]).to_numpy() > 0)
    sub = D.loc[keep]
    z = rank_normal(sub[feat].to_numpy())
    syn = sub["mut_n_syn"].to_numpy().astype(int)
    off = np.zeros(len(sub))
    blocks = make_blocks(sub["chr"].to_numpy(), a.block_bins)
    for split in range(a.n_placebo):
        half_a = rng.binomial(syn, 0.5).astype(float)
        ntot = syn.astype(float)
        est = fit_glm_slope(z, off, "binomial", half_a, ntot)
        boot = bootstrap_slope(lambda r: fit_glm_slope(z[r], off[r], "binomial", half_a[r], ntot[r]),
                               blocks, a.n_boot_placebo, rng)
        s = summarise(est, boot)
        placebo_rows.append(dict(feature=feat, split=split, estimate=s["estimate"], se_boot=s["se_boot"], p=s["p"]))
    log(f"done {feat}")

E = pd.DataFrame(rows_out)
for col in ["group", "mut_syn", "mut_mis", "mut_trunc"]:
    E[f"prior_{col}" if col != "group" else "group"] = E["feature"].map(L[col]) if col in L.columns else np.nan
# BH per effect, over the model features only (controls are reported but not part of the family)
E["q"] = np.nan
for eff, g in E[E.role == "feature"].groupby("effect"):
    E.loc[g.index, "q"] = bh(g["p"].to_numpy())
n_failed = int(E["n_boot_failed"].sum())
if n_failed:
    log(f"WARNING: {n_failed} bootstrap fits did not converge (excluded from SE/CI), see column n_boot_failed")
E.to_csv(f"{a.out_dir}/effects.tsv", sep="\t", index=False)

P = pd.DataFrame(placebo_rows)
P.to_csv(f"{a.out_dir}/placebo.tsv", sep="\t", index=False)
PS = P.groupby("feature").agg(n_placebo=("p", "size"), fpr_05=("p", lambda s: float((s < 0.05).mean())),
                              mean_gamma=("estimate", "mean"), sd_gamma=("estimate", "std")).reset_index()
PS.to_csv(f"{a.out_dir}/placebo_summary.tsv", sep="\t", index=False)
log(f"placebo: overall false-positive rate at p < 0.05 = {(P['p'] < 0.05).mean():.3f} (expected 0.05)")
log(f"written {a.out_dir}/effects.tsv, placebo.tsv, placebo_summary.tsv")
