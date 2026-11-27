"""
PILOT - does a feature's selection effect differ INSIDE driver genes vs non-driver genic bins? (Gab, 27/11/2026)
Extension of pilot_class_interaction_test.py with a driver-bin indicator and its interaction with the feature.

Question: is the per-feature selection signal carried only by drivers, or do drivers show a DIFFERENT feature effect
than the rest of the coding genome?

Model, one feature at a time (z = rank-based normal score of the feature over the bins of the fit), coding bins only:
  * supply   : Poisson    syn ~ offset(log expo_syn) + a + b_z z + b_D D + b_zD z*D
  * selection: Binomial   C | (C + syn), logit p = log(expo_C / expo_syn) + a + g_z z + g_D D + g_zD z*D,  C = mis, trunc
  D = 1 if the bin holds coding sites of a driver gene (gene -> bin link from gene_bin_weighted_sites.tsv of
  context_weighted_sites.py; same rule as haplo_minus_drivers.py), 0 otherwise.
  Terms reported per (feature, effect):
    z        effect of the feature in NON-driver bins
    D        driver excess at z = 0 (selection: positive control, expected > 0, mostly truncating)
    zD       DIFFERENCE of the feature effect inside vs outside drivers  <- the quantity of interest
    z_driver effect of the feature inside drivers (z + zD)
  dist.to.closest.OG / dist.to.closest.TSG define D: never tested as features (circular).

Inference:
  * stratified bootstrap: non-driver bins resampled in contiguous blocks of --block_bins (within chromosomes, as the
    original test); driver bins resampled as CLUSTERS (driver bins of the same chromosome at most --driver_gap bins apart
    = one cluster, so one long gene or two neighbouring drivers are one unit). Each replicate draws as many blocks and
    clusters as there are. SE, 95% percentile CI, p = 2 * P(N(0,1) > |est| / SE).
  * random gene-set null for zD (calibration): D is replaced by --n_null random sets of NON-driver genes, matched one
    by one to the driver genes on total coding sites (strata = --n_strata quantiles of log total sites). Candidate genes
    exclude every IntOGen compendium gene (any cancer type) and the cohort list. Real driver bins are REMOVED from the
    null fits (the reference group is then non-driver bins, as in the real fit). p_null = (1 + #(|null| >= |obs|)) /
    (1 + n_null). Any gene subset differs a little from the rest: the real zD must stand out from this distribution.
    The same gene sets are used for every feature. Point estimates only (no bootstrap inside the null).
  * placebo: synonymous counts split at random (Binomial(syn, 0.5)), same model with offset 0; the zD term must be null.
    fpr_05 = share of splits with bootstrap p < 0.05.
  * power: mde80 = (1.96 + 0.84) * se_boot of zD (Wald approximation: smallest |zD| detected at p < 0.05 with 80% power;
    pilot_spikein_power.py validated this approximation on toy data, 25/11/2026). Interpret a null zD only if the
    expected effect is above mde80.

BH q-values per (effect, term) over the tested features, separately for p (bootstrap) and p_null.

Same preprocessing as pilot_class_interaction_test.py (gene-level NA -> 0, distance NA -> 2 Mb, drop >= 50% NA and
constant features, Hi-C eigenvectors oriented by GC per chromosome). A feature with fewer than --min_driver_levels
distinct values in the driver bins of a fit is 'not_testable' for that fit (zD not identifiable; logged).

Inputs : --dataset        feature matrix (TSV: binID, features, mut_n_syn/mis/trunc, sites_syn/mis/trunc)
         --labels         config/predictor_identity_labels.tsv (annotation only)
         --exposure_file  bin_weighted_sites.tsv (context_weighted_sites.py; wsites_* used in every offset). Optional:
                          without it, sites_* of the dataset.
         --gene_bin_sites gene_bin_weighted_sites.tsv (gene, binID, cls, sites, wsites)
         --cohort_bed     resources/IntOGen/processed/hg38_v45_<C>.sorted.bed (column 5 = gene symbol)
         --compendium     resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv
         --driver_list    cohort (cohort bed, default) | cancer_type (+ compendium genes with CANCER_TYPE == --cohort)
                          | pancancer (+ every compendium gene)
         --exclude_genes  optional gene symbols whose bins are removed from every fit (sensitivity, e.g. TP53 PIK3CA)
Outputs (--out_dir):
  effects.tsv          feature, effect (supply / selection_mis / selection_trunc), term (z / D / zD / z_driver), estimate,
                       se_boot, ci_low, ci_high, p, q, n_boot_failed; for zD also null_mean, null_sd, p_null, q_null,
                       mde80; n_bins, n_driver_bins, sd_z_driver (SD of z inside driver bins: low = little to learn from),
                       prior labels
  null_draws.tsv       feature, effect, draw, estimate (zD with random gene sets)
  null_gene_sets.tsv   draw, gene (the random sets), plus driver_units.tsv (driver genes, bins, cluster, total sites)
  placebo.tsv / placebo_summary.tsv   zD of the synonymous split per feature and split / fpr_05 per feature
Usage (pilot env, repo root):
  python -u pilot_driver_interaction_test.py --dataset <20 kb matrix> --labels config/predictor_identity_labels.tsv \
      --exposure_file analyses/v4/<C>/20000bp/context_weighted_sites_pooled/bin_weighted_sites.tsv \
      --gene_bin_sites analyses/v4/<C>/20000bp/context_weighted_sites_pooled/gene_bin_weighted_sites.tsv \
      --cohort_bed resources/IntOGen/processed/hg38_v45_<C>.sorted.bed \
      --compendium resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv --cohort <C> \
      --out_dir analyses/v4/<C>/20000bp/driver_interaction_cohort --block_bins 500
"""
import argparse, os, time
import numpy as np, pandas as pd
from scipy.stats import norm

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--labels", required=True)
p.add_argument("--exposure_file", default=None)
p.add_argument("--gene_bin_sites", required=True)
p.add_argument("--cohort_bed", required=True)
p.add_argument("--compendium", required=True)
p.add_argument("--cohort", required=True, help="cohort / IntOGen CANCER_TYPE code, e.g. BRCA, LUAD")
p.add_argument("--driver_list", choices=["cohort", "cancer_type", "pancancer"], default="cohort")
p.add_argument("--exclude_genes", nargs="*", default=[])
p.add_argument("--out_dir", required=True)
p.add_argument("--block_bins", type=int, default=500, help="non-driver bins per bootstrap block (contiguous, within chromosome)")
p.add_argument("--driver_gap", type=int, default=1, help="driver bins at most this many bins apart = one bootstrap cluster")
p.add_argument("--n_boot", type=int, default=1000)
p.add_argument("--n_null", type=int, default=200)
p.add_argument("--n_strata", type=int, default=20, help="quantile strata of log total coding sites for gene matching")
p.add_argument("--n_placebo", type=int, default=20)
p.add_argument("--n_boot_placebo", type=int, default=200)
p.add_argument("--min_driver_levels", type=int, default=3)
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
t0 = time.time()
def log(msg): print(f"[{time.time() - t0:.0f}s] {msg}", flush=True)
rng = np.random.default_rng(a.seed)

# ============================ data and feature set (same rules as pilot_class_interaction_test.py) ============================
D = pd.read_csv(a.dataset, sep="\t")
D["chr"] = D["binID"].str.split("_").str[0]
D["idx"] = D["binID"].str.split("_").str[1].astype(int)
chrom_order = lambda v: int(v) if str(v).isdigit() else 100 + ord(str(v)[0])
D = D.sort_values(["chr", "idx"], key=lambda s: s.map(chrom_order) if s.name == "chr" else s).reset_index(drop=True)

needed = ["mut_n_syn", "mut_n_mis", "mut_n_trunc", "sites_syn", "sites_mis", "sites_trunc"]
missing = [c for c in needed if c not in D.columns]
if missing:
    raise SystemExit(f"ERROR: {a.dataset} lacks the mutation columns {missing}: use a feature matrix with COSMIC classes.")

if a.exposure_file:
    X = pd.read_csv(a.exposure_file, sep="\t")
    missing = [c for c in ["binID", "wsites_syn", "wsites_mis", "wsites_trunc"] if c not in X.columns]
    if missing:
        raise SystemExit(f"ERROR: {a.exposure_file} lacks the columns {missing} (expected the output of context_weighted_sites.py)")
    D = D.merge(X[["binID", "wsites_syn", "wsites_mis", "wsites_trunc"]], on="binID", how="left")
    no_value = D["sites_syn"].gt(0) & D["wsites_syn"].isna()
    if no_value.any():
        raise SystemExit(f"ERROR: {int(no_value.sum())} bins with coding sites have no wsites in {a.exposure_file}: "
                         f"different backbone or resolution? e.g. {D.loc[no_value, 'binID'].head(5).tolist()}")
    for c in ("syn", "mis", "trunc"):
        D[f"expo_{c}"] = D[f"wsites_{c}"].fillna(0)
    log(f"exposure: context-weighted sites from {a.exposure_file}")
else:
    for c in ("syn", "mis", "trunc"):
        D[f"expo_{c}"] = D[f"sites_{c}"]
    log("exposure: unweighted sites_* of the dataset")

ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
DRIVER_DEF = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]   # define D: never tested
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_", "expo_"))]
EXONIC = ["cds_bp"] if "cds_bp" in D.columns else []
FEATS = [c for c in D.columns if c not in ALT + DRIVER_DEF + EXPOSURE + EXONIC + ["binID", "chr", "idx"]]
GENE_LEVEL_NA_TO_ZERO = [c for c in FEATS
                         if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                                  "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
                         or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs",
                                          "total_n_paralogs"))]
for c in GENE_LEVEL_NA_TO_ZERO:
    D[c] = D[c].fillna(0)
for c in [c for c in FEATS if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)                                   # NA on a whole chromosome = no gene of that list: cap
too_many_na = [c for c in FEATS if D[c].isna().mean() >= 0.5]
constant = [c for c in FEATS if c not in too_many_na and D[c].nunique() <= 1]
FEATS = [c for c in FEATS if c not in too_many_na + constant]
log(f"{len(D)} bins, {len(FEATS)} features (not tested, they define D: {DRIVER_DEF}); dropped: {len(too_many_na)} with "
    f">= 50% NA {too_many_na}, {len(constant)} constant {constant}")
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

# ============================ driver genes and driver bins ============================
GB = pd.read_csv(a.gene_bin_sites, sep="\t")
missing = [c for c in ["gene", "binID", "cls", "sites"] if c not in GB.columns]
if missing:
    raise SystemExit(f"ERROR: {a.gene_bin_sites} lacks the columns {missing} (expected gene_bin_weighted_sites.tsv)")
GB = GB[GB["sites"] > 0]
gene_sites = GB.groupby("gene")["sites"].sum()                     # total coding sites per gene (all classes, all bins)
gene_bins = GB[["gene", "binID"]].drop_duplicates()
not_in_matrix = set(gene_bins["binID"]) - set(D["binID"])
if not_in_matrix:
    log(f"WARNING: {len(not_in_matrix)} bins of {a.gene_bin_sites} are not in the dataset (ignored)")

cohort_genes = set(pd.read_csv(a.cohort_bed, sep="\t", header=None)[4])
comp = pd.read_csv(a.compendium, sep="\t")
if not {"SYMBOL", "CANCER_TYPE"} <= set(comp.columns):
    raise SystemExit(f"ERROR: {a.compendium} lacks SYMBOL / CANCER_TYPE columns")
if a.driver_list == "cohort":
    drivers = set(cohort_genes)
elif a.driver_list == "cancer_type":
    drivers = cohort_genes | set(comp.loc[comp.CANCER_TYPE == a.cohort, "SYMBOL"])
else:
    drivers = cohort_genes | set(comp.SYMBOL)
not_found = sorted(drivers - set(gene_sites.index))
if not_found:
    log(f"WARNING: {len(not_found)} driver genes without coding sites in {a.gene_bin_sites} (not used): {not_found}")
drivers = drivers & set(gene_sites.index)
any_known_driver = cohort_genes | set(comp.SYMBOL)                 # never drawn into the null gene sets
log(f"driver list '{a.driver_list}': {len(drivers)} genes with coding sites")

# sensitivity: genes removed from every fit
excluded_bins = set()
if a.exclude_genes:
    unknown = [g for g in a.exclude_genes if g not in gene_sites.index]
    if unknown:
        raise SystemExit(f"ERROR: --exclude_genes not found in {a.gene_bin_sites}: {unknown}")
    excluded_bins = set(gene_bins.loc[gene_bins.gene.isin(a.exclude_genes), "binID"])
    drivers = drivers - set(a.exclude_genes)
    n_before = len(D)
    D = D[~D["binID"].isin(excluded_bins)].reset_index(drop=True)
    log(f"--exclude_genes {a.exclude_genes}: {n_before - len(D)} bins removed from every fit")

driver_bins = set(gene_bins.loc[gene_bins.gene.isin(drivers), "binID"])
D["is_driver"] = D["binID"].isin(driver_bins).astype(float)
coding = D["expo_syn"] > 0
log(f"driver bins with coding sites: {int((D.is_driver.astype(bool) & coding).sum())} of {int(coding.sum())} coding bins")

# bootstrap clusters of driver bins (same chromosome, at most --driver_gap bins apart)
D["driver_cluster"] = -1
cid, prev_chr, prev_idx = -1, None, None
for r in np.flatnonzero(D["is_driver"].to_numpy() == 1):
    if D.at[r, "chr"] != prev_chr or D.at[r, "idx"] - prev_idx > a.driver_gap:
        cid += 1
    D.at[r, "driver_cluster"] = cid
    prev_chr, prev_idx = D.at[r, "chr"], D.at[r, "idx"]
log(f"driver bins form {cid + 1} bootstrap clusters (gap <= {a.driver_gap} bins)")

bin_cluster = D.set_index("binID")["driver_cluster"]
units = gene_bins[gene_bins.gene.isin(drivers) & gene_bins.binID.isin(D.binID)].copy()
units["cluster"] = units["binID"].map(bin_cluster)
units = units.groupby("gene").agg(n_bins=("binID", "nunique"), cluster=("cluster", "min")).reset_index()
units["total_sites"] = units["gene"].map(gene_sites)
units.to_csv(f"{a.out_dir}/driver_units.tsv", sep="\t", index=False)

# ============================ random gene-set null: genes matched on total coding sites ============================
candidates = gene_sites[~gene_sites.index.isin(any_known_driver)]
if a.exclude_genes:
    candidates = candidates[~candidates.index.isin(a.exclude_genes)]
log_all = np.log(pd.concat([candidates, gene_sites[list(drivers)]]))
edges = np.unique(np.quantile(log_all, np.linspace(0, 1, a.n_strata + 1)))
stratum_of = lambda s: np.clip(np.searchsorted(edges, np.log(s), side="right") - 1, 0, len(edges) - 2)
cand_stratum = pd.Series(stratum_of(candidates.to_numpy()), index=candidates.index)
driver_strata = pd.Series(stratum_of(gene_sites[list(drivers)].to_numpy())).value_counts()
for k, n_need in driver_strata.items():
    n_have = int((cand_stratum == k).sum())
    if n_have < n_need:
        raise SystemExit(f"ERROR: stratum {k} of total coding sites has {n_need} driver genes but only {n_have} candidate "
                         f"genes: lower --n_strata")
cand_by_stratum = {k: cand_stratum.index[cand_stratum == k].to_numpy() for k in driver_strata.index}
null_sets, null_rows = [], []
for draw in range(a.n_null):
    genes = np.concatenate([rng.choice(cand_by_stratum[k], n, replace=False) for k, n in driver_strata.items()])
    null_sets.append(set(genes))
    null_rows += [dict(draw=draw, gene=g) for g in genes]
pd.DataFrame(null_rows).to_csv(f"{a.out_dir}/null_gene_sets.tsv", sep="\t", index=False)
# null indicator per draw over the rows of D; real driver bins are removed from the null fits (NaN marks them)
null_D = np.empty((a.n_null, len(D)))
for draw, genes in enumerate(null_sets):
    bins = set(gene_bins.loc[gene_bins.gene.isin(genes), "binID"])
    null_D[draw] = D["binID"].isin(bins).to_numpy(float)
null_D[:, D["is_driver"].to_numpy() == 1] = np.nan
n_null_bins = np.nansum(null_D * coding.to_numpy()[None, :], axis=1)
log(f"null: {a.n_null} random gene sets ({len(drivers)} genes each, matched on total coding sites in {len(edges) - 1} "
    f"strata); coding bins per set median {np.median(n_null_bins):.0f} (real {int((D.is_driver.astype(bool) & coding).sum())})")

# ============================ GLM with several covariates, Newton-Raphson ============================
def fit_glm(Xm, offset, family, y, n=None, max_iter=100, tol=1e-10):
    """Coefficients of a GLM (first column of Xm = intercept) with an offset. family 'poisson' or 'binomial' (y successes
    of n trials, logit). Returns an array of NaN if Newton does not converge or the information matrix is singular."""
    b = np.zeros(Xm.shape[1])
    if family == "poisson":
        b[0] = np.log(y.sum() / np.exp(offset).sum())
    else:
        frac = y.sum() / n.sum()
        b[0] = np.log(frac / (1 - frac)) - np.mean(offset)
    for _ in range(max_iter):
        eta = offset + Xm @ b
        if family == "poisson":
            mu = np.exp(eta); grad_w = y - mu; w = mu
        else:
            prob = 1 / (1 + np.exp(-eta)); grad_w = y - n * prob; w = n * prob * (1 - prob)
        try:
            step = np.linalg.solve(Xm.T @ (Xm * w[:, None]), Xm.T @ grad_w)
        except np.linalg.LinAlgError:
            return np.full(Xm.shape[1], np.nan)
        b = b + step
        if not np.all(np.isfinite(b)):
            return np.full(Xm.shape[1], np.nan)
        if np.max(np.abs(step)) < tol:
            return b
    return np.full(Xm.shape[1], np.nan)

def design(z, d):
    return np.column_stack([np.ones_like(z), z, d, z * d])

TERMS = ["z", "D", "zD"]          # coefficient positions 1, 2, 3 of design(); z_driver = z + zD

def rank_normal(x):
    r = pd.Series(x).rank(method="average").to_numpy()
    return norm.ppf((r - 0.5) / len(r))

def make_units(chrom, is_drv, drv_cluster, block_bins):
    """Bootstrap units: contiguous blocks of non-driver rows (within chromosomes) and driver clusters."""
    nd_blocks, d_units = [], []
    for ch in pd.unique(chrom):
        pos = np.flatnonzero((chrom == ch) & (is_drv == 0))
        nd_blocks += [pos[s:s + block_bins] for s in range(0, len(pos), block_bins)]
    for c in np.unique(drv_cluster[is_drv == 1]):
        d_units.append(np.flatnonzero((drv_cluster == c) & (is_drv == 1)))
    return nd_blocks, d_units

def bootstrap(fit_rows, nd_blocks, d_units, n_boot, n_coef):
    out = np.empty((n_boot, n_coef))
    for i in range(n_boot):
        pick_nd = rng.integers(0, len(nd_blocks), len(nd_blocks))
        pick_d = rng.integers(0, len(d_units), len(d_units))
        rows = np.concatenate([nd_blocks[j] for j in pick_nd] + [d_units[j] for j in pick_d])
        out[i] = fit_rows(rows)
    return out

def summarise(est, boot):
    good = boot[np.isfinite(boot)]
    n_bad = len(boot) - len(good)
    if len(good) < 0.9 * len(boot):
        return dict(estimate=est, se_boot=np.nan, ci_low=np.nan, ci_high=np.nan, p=np.nan, n_boot_failed=n_bad)
    se = good.std(ddof=1)
    pval = 2 * norm.sf(abs(est) / se) if se > 0 and np.isfinite(est) else np.nan
    return dict(estimate=est, se_boot=se, ci_low=np.quantile(good, 0.025), ci_high=np.quantile(good, 0.975),
                p=pval, n_boot_failed=n_bad)

def bh(pvals):
    p_arr = np.asarray(pvals, float); m = len(p_arr)
    order = np.argsort(p_arr)
    q = p_arr[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m); out[order] = np.minimum(q, 1)
    return out

# ============================ tests ============================
# effect name -> (family, response column, trials: None = Poisson, offset numerator, offset denominator)
EFFECTS = {"supply": ("poisson", "mut_n_syn", None, "expo_syn", None),
           "selection_mis": ("binomial", "mut_n_mis", "mut_n_syn", "expo_mis", "expo_syn"),
           "selection_trunc": ("binomial", "mut_n_trunc", "mut_n_syn", "expo_trunc", "expo_syn")}
rows_out, null_out, placebo_rows = [], [], []
for feat in FEATS:
    ok_feat = D[feat].notna().to_numpy()
    n_na = int((~ok_feat).sum())
    if n_na:
        log(f"{feat}: {n_na} bins with NA dropped for this feature")
    for eff, (family, ycol, other, num, den) in EFFECTS.items():
        keep = ok_feat & (D["expo_syn"].to_numpy() > 0)
        if family == "binomial":
            keep &= (D[num].to_numpy() > 0) & ((D[ycol] + D[other]).to_numpy() > 0)
        sub = D.loc[keep]
        z = rank_normal(sub[feat].to_numpy())
        d = sub["is_driver"].to_numpy()
        n_drv = int(d.sum())
        sd_z_drv = float(np.std(z[d == 1])) if n_drv > 1 else np.nan
        base = dict(feature=feat, effect=eff, n_bins=int(keep.sum()), n_driver_bins=n_drv, sd_z_driver=sd_z_drv,
                    n_dropped_na=n_na)
        if len(np.unique(z[d == 1])) < a.min_driver_levels or len(np.unique(z[d == 0])) < 2:
            log(f"WARNING: {feat} {eff}: {len(np.unique(z[d == 1]))} distinct values in the {n_drv} driver bins: "
                f"zD not identifiable, not tested")
            rows_out.append(dict(**base, term="not_testable"))
            continue
        y = sub[ycol].to_numpy(float)
        if family == "poisson":
            ntot = None
            off = np.log(sub[num].to_numpy(float))
        else:
            ntot = y + sub[other].to_numpy(float)
            off = np.log(sub[num].to_numpy(float) / sub[den].to_numpy(float))
        Xm = design(z, d)
        fit_rows = lambda r: fit_glm(Xm[r], off[r], family, y[r], None if ntot is None else ntot[r])
        est = fit_rows(np.arange(len(sub)))
        nd_blocks, d_units = make_units(sub["chr"].to_numpy(), d, sub["driver_cluster"].to_numpy(), a.block_bins)
        boot = bootstrap(fit_rows, nd_blocks, d_units, a.n_boot, Xm.shape[1])

        # random gene-set null for zD: real driver bins out, random genes' bins as D
        nd = null_D[:, np.flatnonzero(keep)]
        null_est = np.full(a.n_null, np.nan)
        for draw in range(a.n_null):
            r = np.flatnonzero(np.isfinite(nd[draw]))
            dn = nd[draw, r]
            if len(np.unique(z[r][dn == 1])) < a.min_driver_levels:
                continue
            null_est[draw] = fit_glm(design(z[r], dn), off[r], family, y[r], None if ntot is None else ntot[r])[3]
            null_out.append(dict(feature=feat, effect=eff, draw=draw, estimate=null_est[draw]))
        good_null = null_est[np.isfinite(null_est)]

        for k, term in enumerate(TERMS, start=1):
            row = dict(**base, term=term, **summarise(est[k], boot[:, k]))
            if term == "zD":
                row["mde80"] = (norm.ppf(0.975) + norm.ppf(0.8)) * row["se_boot"]
                row["n_null_ok"] = len(good_null)
                if len(good_null) >= 0.9 * a.n_null and np.isfinite(est[k]):
                    row["null_mean"] = good_null.mean(); row["null_sd"] = good_null.std(ddof=1)
                    row["p_null"] = (1 + np.sum(np.abs(good_null) >= abs(est[k]))) / (1 + len(good_null))
                else:
                    log(f"WARNING: {feat} {eff}: only {len(good_null)} of {a.n_null} null fits usable: p_null not computed")
            rows_out.append(row)
        rows_out.append(dict(**base, term="z_driver", **summarise(est[1] + est[3], boot[:, 1] + boot[:, 3])))

    # ---- placebo: random half of synonymous vs the other half (same sites -> offset 0), zD term ----
    keep = ok_feat & (D["mut_n_syn"].to_numpy() > 0)
    sub = D.loc[keep]
    z = rank_normal(sub[feat].to_numpy())
    d = sub["is_driver"].to_numpy()
    if len(np.unique(z[d == 1])) >= a.min_driver_levels:
        syn = sub["mut_n_syn"].to_numpy().astype(int)
        off = np.zeros(len(sub)); ntot = syn.astype(float); Xm = design(z, d)
        nd_blocks, d_units = make_units(sub["chr"].to_numpy(), d, sub["driver_cluster"].to_numpy(), a.block_bins)
        for split in range(a.n_placebo):
            half_a = rng.binomial(syn, 0.5).astype(float)
            fr = lambda r: fit_glm(Xm[r], off[r], "binomial", half_a[r], ntot[r])
            est = fr(np.arange(len(sub)))
            s = summarise(est[3], bootstrap(fr, nd_blocks, d_units, a.n_boot_placebo, Xm.shape[1])[:, 3])
            placebo_rows.append(dict(feature=feat, split=split, estimate=s["estimate"], se_boot=s["se_boot"], p=s["p"]))
    log(f"done {feat}")

E = pd.DataFrame(rows_out)
for col in ["group", "mut_syn", "mut_mis", "mut_trunc"]:
    if col in L.columns:
        E["group" if col == "group" else f"prior_{col}"] = E["feature"].map(L[col])
E["q"] = np.nan; E["q_null"] = np.nan
for (eff, term), g in E[E["term"] != "not_testable"].groupby(["effect", "term"]):
    g_p = g[g["p"].notna()]
    if len(g_p):
        E.loc[g_p.index, "q"] = bh(g_p["p"].to_numpy())
    if "p_null" in g.columns:
        g_n = g[g["p_null"].notna()]
        if len(g_n):
            E.loc[g_n.index, "q_null"] = bh(g_n["p_null"].to_numpy())
untrusted = E[(E["term"] != "not_testable") & E["p"].isna()]
if len(untrusted):
    log(f"WARNING: {len(untrusted)} rows with < 90% converged bootstrap replicates (SE/CI/p = NaN): "
        + ", ".join(f"{f} {e} {t}" for f, e, t in zip(untrusted.feature, untrusted.effect, untrusted.term)))
E.to_csv(f"{a.out_dir}/effects.tsv", sep="\t", index=False)
pd.DataFrame(null_out).to_csv(f"{a.out_dir}/null_draws.tsv", sep="\t", index=False)

P = pd.DataFrame(placebo_rows, columns=["feature", "split", "estimate", "se_boot", "p"])
P.to_csv(f"{a.out_dir}/placebo.tsv", sep="\t", index=False)
PS = P.groupby("feature").agg(n_placebo=("p", "count"), fpr_05=("p", lambda s: float((s.dropna() < 0.05).mean())),
                              mean_zD=("estimate", "mean"), sd_zD=("estimate", "std")).reset_index()
PS.to_csv(f"{a.out_dir}/placebo_summary.tsv", sep="\t", index=False)
log(f"placebo zD: false-positive rate at p < 0.05 = {(P['p'].dropna() < 0.05).mean():.3f} (expected 0.05)")
pos = E[(E.term == "D") & (E.effect != "supply")].groupby("effect")["estimate"].median()
log("positive control, median driver excess (term D) over features: " + ", ".join(f"{k} {v:+.3f}" for k, v in pos.items()))
log(f"written {a.out_dir}/effects.tsv, null_draws.tsv, null_gene_sets.tsv, driver_units.tsv, placebo*.tsv")
