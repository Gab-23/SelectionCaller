"""
PILOT - spike-in power calibration of the per-feature selection test (pilot_class_interaction_test.py).
Gab, 25/11/2026.

Question: for each predictor, could the selection test detect an effect of a given size in data like ours? A null result
means "no selection" only if the test had power for the effect sizes we care about.

Design (agreed 25/11/2026; inference changed the same day, see NOTE below):
  1) REALISED EFFECT by Poisson-resampling injection, one predictor and one class (mis / trunc) at a time:
        C'_i ~ Poisson( C_i * exp((gamma - gamma_obs) * z_i) / c )
    C_i = observed count of the class in bin i, z_i = rank-normal score of the predictor, gamma_obs = the predictor's own
    fitted selection slope on the observed data (REMOVED, so the injected effect replaces it), c = rescaling so that the
    total count of the class is unchanged in expectation. Synonymous counts and exposures stay as observed. Each replicate
    is fitted as in pilot_class_interaction_test.py (binomial C | (C + syn), offset log(expo_C / expo_syn), z re-ranked on
    the bins kept); realised effect m(gamma) = mean slope over --n_rep replicates. m(gamma) < gamma means the effect is
    diluted by the data (e.g. bins with C_i = 0 cannot receive mutations, frequent at 20 kb).
  2) UNCERTAINTY from the real run: SE = se_boot of the predictor and class in --effects (moving-block bootstrap of the
    real data, so it includes overdispersion and autocorrelation).
  3) POWER, two detection criteria, in the injected direction:
      - p < 0.05:  power = Phi(m/SE - 1.96) for gamma > 0, Phi(-m/SE - 1.96) for gamma < 0 (normal approximation, the
                   same one behind the p-values of the main test);
      - BH q < 0.05 in the discovery setting: --n_draw estimates drawn from N(m, SE^2) -> p -> BH together with the
                   p-values of the OTHER predictors for the same class from the real run (model features with a trusted
                   p, the real BH family). Approximation: the other predictors keep their observed p.
    At gamma = 0 the reported rate is 5% (p) by construction: it is NOT a calibration check; calibration is the placebo
    of the main test.
  4) Minimum detectable effect (MDE): smallest |gamma| of the grid with power >= --power_target (default 0.8), per
    predictor, class, sign and criterion; NaN = not reached within the grid (--gammas).

NOTE (25/11/2026): first version tested each resampled replicate with its own block bootstrap. On toy data the bootstrap SE
inside the replicates was ~sqrt(2) times the real SE (resampling adds a second round of Poisson noise), so p-values were too
large: false-positive rate 0.00 at gamma = 0 and MDE ~1.5x too large (power conservative). Replaced by the scheme above
(agreed with Gab): resampling only for the realised effect, uncertainty from the real run. Much faster (no bootstrap).

Limits (read before interpreting):
  * Normal approximation for the estimate with the real-run SE (same assumption as the p-values of the main test).
  * The SE is the one of the observed data; with a large injected effect the true SE changes slightly (ignored).
  * The injected effect is log-linear in the rank-normal score, the same shape the test fits: power for effects of other
    shapes (e.g. concentrated in a few bins, like OG hotspots) can be much lower.
  * Removing gamma_obs removes only the linear part of the predictor's own effect.
  * Predictor/class pairs without a trusted SE in --effects (p = NaN there) are skipped and logged.

Preprocessing: IDENTICAL to pilot_class_interaction_test.py (copied block below, keep the two in sync): gene-level NA -> 0,
distance NA -> 2 Mb, drop >= 50% NA and constant features, Hi-C eigenvectors oriented by GC per chromosome,
--exposure_file for context-weighted sites. Features constant on bins with coding sites are skipped (not testable).

Inputs : --dataset        feature matrix (as for pilot_class_interaction_test.py)
         --effects        effects.tsv of pilot_class_interaction_test.py on the SAME dataset / exposure (for BH and as
                          a check that the observed slope matches)
         --exposure_file  optional bin_weighted_sites.tsv (use it if the real run used it)
         --features       optional comma-separated subset (default: all features + positive controls)
Outputs (--out_dir):
  spikein_replicates.tsv  feature, role, cls, gamma, rep, estimate, n_bins   (slopes of the resampled replicates)
  spikein_power.tsv       feature, role, cls, gamma, gamma_obs, se_real, n_rep, n_valid, mean_estimate (= m),
                          sd_estimate, realised_fraction (m / gamma), power_p05, power_q05
  spikein_mde.tsv         feature, role, cls, sign, se_real, mde_p05, mde_q05
Usage (cluster, pilot env):
  python pilot_spikein_power.py --dataset <feature_matrix> --effects analyses/v4/LUAD/1000000bp/class_interaction_wsites_pooled/effects.tsv \
      --exposure_file analyses/v4/LUAD/1000000bp/context_weighted_sites_pooled/bin_weighted_sites.tsv \
      --out_dir analyses/v4/LUAD/1000000bp/spikein_power --n_rep 200 --n_jobs 8
"""
import argparse, os, time
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.stats import norm

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--effects", required=True, help="effects.tsv of pilot_class_interaction_test.py on the same data")
p.add_argument("--out_dir", required=True)
p.add_argument("--exposure_file", default=None)
p.add_argument("--features", default=None, help="comma-separated subset of features/controls (default: all)")
p.add_argument("--gammas", default="0.01,0.02,0.035,0.05,0.075,0.1",
               help="positive injected effects; their negatives and 0 are added")
p.add_argument("--n_rep", type=int, default=200, help="spike-in replicates per feature, class and gamma (realised effect)")
p.add_argument("--n_draw", type=int, default=5000, help="normal draws per feature, class and gamma for the BH power")
p.add_argument("--power_target", type=float, default=0.8)
p.add_argument("--n_jobs", type=int, default=1, help="parallel processes (one feature per task)")
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
t0 = time.time()
def log(msg): print(f"[{time.time() - t0:.0f}s] {msg}", flush=True)

pos_gammas = sorted({float(g) for g in a.gammas.split(",")})
if any(g <= 0 for g in pos_gammas):
    raise SystemExit(f"ERROR: --gammas must be positive values (their negatives are added automatically), got {a.gammas}")
GAMMAS = [0.0] + pos_gammas + [-g for g in pos_gammas]

# ======== data and feature set: SAME rules as pilot_class_interaction_test.py (copied, keep in sync) ========
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
CONTROLS = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_", "expo_"))]
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
    D[c] = D[c].fillna(2e6)
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

TO_RUN = FEATS + CONTROLS
if a.features:
    asked = [f.strip() for f in a.features.split(",") if f.strip()]
    unknown = [f for f in asked if f not in TO_RUN]
    if unknown:
        raise SystemExit(f"ERROR: --features not among the testable features/controls of this dataset: {unknown}")
    TO_RUN = asked

# ======== real-run p-values (BH family) ========
E = pd.read_csv(a.effects, sep="\t")
missing = [c for c in ["feature", "role", "effect", "estimate", "p"] if c not in E.columns]
if missing:
    raise SystemExit(f"ERROR: {a.effects} lacks the columns {missing} (expected effects.tsv of pilot_class_interaction_test.py)")
if "se_boot" not in E.columns:
    raise SystemExit(f"ERROR: {a.effects} lacks the column se_boot (expected effects.tsv of pilot_class_interaction_test.py)")
REAL_P, REAL_EST, REAL_SE = {}, {}, {}
for cls in ("mis", "trunc"):
    fam = E[(E.role == "feature") & (E.effect == f"selection_{cls}") & E["p"].notna()]
    REAL_P[cls] = fam.set_index("feature")["p"]
    sel = E[E.effect == f"selection_{cls}"].set_index("feature")
    REAL_EST[cls] = sel["estimate"]
    REAL_SE[cls] = sel["se_boot"].where(sel["p"].notna())      # untrusted bootstrap (p NaN) -> no SE
    log(f"BH family for {cls}: {len(fam)} features with a trusted p in {a.effects}")

# ======== test functions: same as pilot_class_interaction_test.py ========
def fit_glm_slope(z, offset, y, n, max_iter=50, tol=1e-10):
    """Binomial (logit) GLM with intercept, one covariate z and an offset; returns the slope (NaN if no convergence)."""
    X = np.column_stack([np.ones_like(z), z])
    frac = y.sum() / n.sum()
    b = np.array([np.log(frac / (1 - frac)) - np.mean(offset), 0.0])
    for _ in range(max_iter):
        prob = 1 / (1 + np.exp(-(offset + X @ b)))
        try:
            step = np.linalg.solve(X.T @ (X * (n * prob * (1 - prob))[:, None]), X.T @ (y - n * prob))
        except np.linalg.LinAlgError:
            return np.nan
        b = b + step
        if np.max(np.abs(step)) < tol:
            return b[1]
    return np.nan

def rank_normal(x):
    r = pd.Series(x).rank(method="average").to_numpy()
    return norm.ppf((r - 0.5) / len(r))

def bh_threshold(p_others, alpha=0.05):
    """Largest p that a NEW test would need to be called at BH q < alpha, when added to the fixed set p_others.
    BH with m tests calls every p <= alpha * R / m, R = largest r with #{p <= alpha * r / m} >= r. With the new p
    counted in that set, it is called iff for some r: p <= alpha * r / m and (#others <= alpha * r / m) + 1 >= r.
    The threshold is the largest such alpha * r / m (0 if none)."""
    others = np.asarray(p_others, float)
    m = len(others) + 1
    thr = 0.0
    for r in range(1, m + 1):
        cut = alpha * r / m
        if (others <= cut).sum() + 1 >= r:
            thr = cut
    return thr

def fit_once(y, syn, off, x_raw):
    """Slope of the selection test (no bootstrap) on the bins with y + syn > 0, z re-ranked on those bins."""
    keep = (y + syn) > 0
    z = rank_normal(x_raw[keep])
    return fit_glm_slope(z, off[keep], y[keep], (y + syn)[keep]), int(keep.sum())

# ======== one feature: all classes, gammas, replicates ========
def run_feature(job):
    k, feat = job
    rng = np.random.default_rng([a.seed, k])
    role = "positive_control" if feat in CONTROLS else "feature"
    rows, msgs = [], []
    ok_feat = D[feat].notna().to_numpy()
    has_sites = ok_feat & (D["sites_syn"].to_numpy() > 0)
    if D.loc[has_sites, feat].nunique() <= 1:
        return rows, [f"WARNING: {feat} is constant on the bins with coding sites: not testable, skipped"]
    for cls in ("mis", "trunc"):
        ycol, scol = f"mut_n_{cls}", f"expo_{cls}"
        base = ok_feat & (D["expo_syn"].to_numpy() > 0) & (D[scol].to_numpy() > 0)
        sub = D.loc[base]
        c_obs = sub[ycol].to_numpy(float); syn = sub["mut_n_syn"].to_numpy(float)
        off = np.log(sub[scol].to_numpy(float) / sub["expo_syn"].to_numpy(float))
        x_raw = sub[feat].to_numpy(float)
        z_base = rank_normal(x_raw)
        se_real = REAL_SE[cls].get(feat, np.nan)
        if not np.isfinite(se_real) or se_real <= 0:
            msgs.append(f"WARNING: {feat} {cls}: no trusted bootstrap SE in {a.effects} (p NaN there): skipped")
            continue
        # the predictor's own slope on the observed data (bins with C + syn > 0, as in the real run)
        has_n = (c_obs + syn) > 0
        gamma_obs = fit_glm_slope(rank_normal(x_raw[has_n]), off[has_n], c_obs[has_n], (c_obs + syn)[has_n])
        if not np.isfinite(gamma_obs):
            msgs.append(f"WARNING: {feat} {cls}: own slope did not converge, skipped")
            continue
        real = REAL_EST[cls].get(feat, np.nan)
        if np.isfinite(real) and abs(real - gamma_obs) > 1e-6:
            msgs.append(f"WARNING: {feat} {cls}: own slope here {gamma_obs:+.5f} differs from {a.effects} ({real:+.5f}): "
                        f"same dataset / exposure / feature preprocessing?")
        for gamma in GAMMAS:
            mult = c_obs * np.exp((gamma - gamma_obs) * z_base)
            lam = mult * c_obs.sum() / mult.sum()          # total of the class unchanged in expectation
            for rep in range(a.n_rep):
                y = rng.poisson(lam).astype(float)
                est, nb = fit_once(y, syn, off, x_raw)
                rows.append(dict(feature=feat, role=role, cls=cls, gamma=gamma, gamma_obs=gamma_obs, se_real=se_real,
                                 rep=rep, estimate=est, n_bins=nb))
    msgs.append(f"done {feat}")
    return rows, msgs

log(f"{len(TO_RUN)} features x 2 classes x {len(GAMMAS)} gammas x {a.n_rep} replicates; {a.n_jobs} processes")
all_rows = []
jobs = list(enumerate(TO_RUN))
if a.n_jobs > 1:
    with Pool(a.n_jobs) as pool:
        for rows, msgs in pool.imap_unordered(run_feature, jobs):
            all_rows += rows
            for m in msgs: log(m)
else:
    for job in jobs:
        rows, msgs = run_feature(job)
        all_rows += rows
        for m in msgs: log(m)
if not all_rows:
    raise SystemExit("ERROR: no feature could be run (see warnings above)")

R = pd.DataFrame(all_rows)
R.to_csv(f"{a.out_dir}/spikein_replicates.tsv", sep="\t", index=False)
n_nan = int(R["estimate"].isna().sum())
if n_nan:
    log(f"WARNING: {n_nan} of {len(R)} replicate fits did not converge: excluded from the realised effect (see n_valid)")

# ======== power and MDE ========
Z_CRIT = norm.isf(0.025)                                   # 1.96: two-sided p < 0.05
BH_THR = {cls: bh_threshold(REAL_P[cls].to_numpy()) for cls in ("mis", "trunc")}
draw_rng = np.random.default_rng(a.seed + 1)

def power_p(m, se, gamma):
    """P(p < 0.05 and estimate in the injected direction) for estimate ~ N(m, se^2); gamma = 0: either direction."""
    up, down = norm.sf(Z_CRIT - m / se), norm.cdf(-Z_CRIT - m / se)
    return up if gamma > 0 else down if gamma < 0 else up + down

def power_q(m, se, gamma, feat, cls):
    """Same with BH q < 0.05 among the real family (the predictor's own real p replaced by the draw)."""
    thr = BH_THR[cls] if feat not in REAL_P[cls].index else bh_threshold(REAL_P[cls].drop(feat).to_numpy())
    est = draw_rng.normal(m, se, a.n_draw)
    pval = 2 * norm.sf(np.abs(est) / se)
    sig = pval <= thr
    if gamma > 0: sig &= est > 0
    if gamma < 0: sig &= est < 0
    return sig.mean()

pw_rows = []
for (feat, role, cls, gamma), g in R.groupby(["feature", "role", "cls", "gamma"]):
    v = g[g["estimate"].notna()]
    se = g["se_real"].iloc[0]
    m = v["estimate"].mean() if len(v) else np.nan
    pw_rows.append(dict(feature=feat, role=role, cls=cls, gamma=gamma, gamma_obs=g["gamma_obs"].iloc[0], se_real=se,
                        n_rep=len(g), n_valid=len(v), mean_estimate=m, sd_estimate=v["estimate"].std(),
                        realised_fraction=m / gamma if gamma != 0 else np.nan,
                        power_p05=power_p(m, se, gamma) if len(v) else np.nan,
                        power_q05=power_q(m, se, gamma, feat, cls) if len(v) else np.nan))
PW = pd.DataFrame(pw_rows).sort_values(["feature", "cls", "gamma"])
PW.to_csv(f"{a.out_dir}/spikein_power.tsv", sep="\t", index=False)

mde_rows = []
for (feat, role, cls), g in PW.groupby(["feature", "role", "cls"]):
    for sign, s in (("positive", 1), ("negative", -1)):
        side = g[np.sign(g.gamma) == s].assign(abs_gamma=lambda d: d.gamma.abs()).sort_values("abs_gamma")
        mde = {}
        for col in ("power_p05", "power_q05"):
            hit = side[side[col] >= a.power_target]
            mde[col] = hit["abs_gamma"].iloc[0] if len(hit) else np.nan
        mde_rows.append(dict(feature=feat, role=role, cls=cls, sign=sign, se_real=g["se_real"].iloc[0],
                             mde_p05=mde["power_p05"], mde_q05=mde["power_q05"]))
MDE = pd.DataFrame(mde_rows)
MDE.to_csv(f"{a.out_dir}/spikein_mde.tsv", sep="\t", index=False)

g0 = PW[PW.gamma == 0]
log(f"check, gamma = 0 (own effect removed): mean realised slope {g0.mean_estimate.mean():+.4f} "
    f"(should be ~0; max |.| {g0.mean_estimate.abs().max():.4f})")
nz = PW[PW.gamma != 0]
log(f"realised fraction m / gamma: median {nz.realised_fraction.median():.2f}, min {nz.realised_fraction.min():.2f} "
    f"(< 1 = effect diluted by the data)")
for cls in ("mis", "trunc"):
    m = MDE[(MDE.cls == cls) & (MDE.role == "feature")]
    log(f"{cls}: features reaching {a.power_target:.0%} power within the grid (p < 0.05): "
        f"{int(m.mde_p05.notna().sum())} of {len(m)} feature x sign; median MDE {m.mde_p05.median():.3f}")
log(f"written {a.out_dir}/spikein_replicates.tsv, spikein_power.tsv, spikein_mde.tsv")
