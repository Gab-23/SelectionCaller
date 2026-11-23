#!/usr/bin/env python
"""
PILOT - identity of the predictors from two chained models (background on a neutral class, selection on top), mutations.

Idea (discussed with Gab, 23/11/2026): background (occurrence) acts on the mutation SUPPLY of a bin, which is shared by
every consequence class; selection acts on the CONSEQUENCE. Synonymous mutations are (nearly) neutral, so:
  1. background model: XGBoost on SYNONYMOUS counts only (offset log(sites_syn)) -> features important here = background;
  2. selection model:  XGBoost on the CONTRAST class only (missense or truncating), offset = log(sites of the class) +
     the background log-rate predicted by model 1 for that bin -> it only sees what the background does not explain;
     features important here = selection.
  label: O = important only in 1, S = only in 2, both = in 1 and 2, none = neither.
This is the same structure Gab plans for the final method (background model -> offset of the selection model), used
here to assign the identity of the predictors.
NOTE (23/11/2026): first written as ONE stacked model (two rows per bin + class indicator, selection = feature x class
interaction). Rejected after the toy test: with a real interaction present XGBoost splits on the class indicator and
every feature, even pure noise and every shadow, got a "significant" interaction (64 of 64 noise features called S).

Details
  bins: cds_bp >= --min_exonic_bp and sites > 0 for both classes (as the pilot).
  CV: interleaved blocks with buffer on the bins (INTERLEAVED_KFOLD_IDXS), the same folds for both models.
  The background prediction used as offset is OUT-OF-FOLD (each bin gets the log-rate from the background model that did
  not see it), so the selection model never sees a background fitted on its own synonymous counts. Limit: the offsets of
  the training bins of a fold come from background models that saw that fold's test bins (synonymous counts only).
  Importance: mean |SHAP| (log-rate scale) on the held-out bins of each fold, mean over --seeds.
  Shadows: one per feature (the feature circularly shifted within each chromosome: same distribution and autocorrelation,
  no alignment with the target), inside BOTH models. A feature is important in a model only if its importance is above
  the maximum importance of the shadows in that model (Boruta rule): this is the floor of generic overfitting.
  Background under-fit (null of the selection model): if model 1 captures only part of a strong background feature, the
  rest appears in model 2 as false selection. Parametric bootstrap: a background model fitted on ALL synonymous counts
  generates synonymous and contrast counts with NO selection (contrast = sites x background rate x class mean); the whole
  procedure (both models, folds, seeds) is rerun on --n_null simulated datasets. p_sel_j = (1 + #null with importance >=
  observed) / (1 + n_null), uncorrected. Selection is called only if the feature beats its shadows AND p_sel < --alpha.
  direction_j (selection model) = Spearman(feature j, SHAP_j) on held-out bins: > 0 = higher values raise the contrast
  class above the background (excess), < 0 = depletion (selection in both directions).
Assumption (cannot be tested here): synonymous mutations are neutral and share the supply of the contrast class once the
exposure (possible sites) is accounted for. Supply differences the exposure misses (trinucleotide context, signatures)
appear as "selection". Correlated features share credit in XGBoost: read the table by blocks of correlated features.

Inputs
  --dataset   feature matrix (data/feature_matrix/<name>), one resolution (start at 1 Mb)
  --labels    optional: prior labels (config/predictor_identity_labels.tsv), only for the comparison table
Outputs (in --out_dir)
  chained_identity.tsv     one row per (contrast, feature, incl. shadows): imp_bg, imp_sel (mean and sd over seeds),
                           shadow thresholds, null mean of imp_sel, p_sel, z_sel, direction, label, prior label
  chained_null.tsv         imp_sel per (contrast, feature, null repeat)
  model_fit.tsv            per contrast and seed: held-out deviance explained of the background model (synonymous, vs
                           constant rate), of the background offset alone on the contrast class, and of the selection
                           model on top of it
  label_vs_prior_<contrast>.tsv  crosstab of the new labels vs the prior labels (only with --labels)
Usage (cluster, pilot env, from the root of the repository)
  python pilot_interaction_identity.py --dataset data/feature_matrix/<1 Mb dataset> --gene_na_to_zero \
      --labels config/predictor_identity_labels.tsv --block_bins 2 --max_buffer_bins 2 --out_dir <dir>
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr

sys.path.insert(0, os.getcwd())
try:
    from workflow.scripts.helpers.INTERLEAVED_KFOLD_IDXS import INTERLEAVED_KFOLD_IDXS
except ImportError:
    sys.exit("workflow/scripts/helpers/INTERLEAVED_KFOLD_IDXS.py is needed: run from the root of the repository.")

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--labels", default=None, help="optional prior labels, only for the comparison table")
p.add_argument("--contrasts", nargs="+", default=["mut_mis", "mut_trunc"], help="classes modelled on top of mut_syn")
p.add_argument("--gene_na_to_zero", action="store_true", help="as in pilot_predictor_identity.py")
p.add_argument("--min_exonic_bp", type=float, default=500)
p.add_argument("--n_splits", type=int, default=5)
p.add_argument("--block_bins", type=int, default=2, help="2 at 1 Mb, 100 at 20 kb (as the identity runs)")
p.add_argument("--max_buffer_bins", type=int, default=2)
p.add_argument("--n_trees", type=int, default=300)
p.add_argument("--n_jobs", type=int, default=4)
p.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
p.add_argument("--n_null", type=int, default=20, help="parametric-bootstrap reruns (null of the selection model)")
p.add_argument("--alpha", type=float, default=0.05, help="uncorrected empirical p threshold of the selection model")
p.add_argument("--random_seed", type=int, default=1, help="shadows and simulated null counts")
p.add_argument("--out_dir", default="chained_identity")
a = p.parse_args()
if 1 / (1 + a.n_null) >= a.alpha:
    sys.exit(f"--n_null {a.n_null}: the smallest possible p_sel is 1/(1+n_null) = {1 / (1 + a.n_null):.3f} >= --alpha "
             f"{a.alpha}: no feature could ever be called S; use n_null >= {int(np.ceil(1 / a.alpha))}")
os.makedirs(a.out_dir, exist_ok=True)
t0 = time.time()
def log(msg): print(f"[{time.time() - t0:.0f}s] {msg}", flush=True)

TARGETS = {"mut_syn": ("mut_n_syn", "sites_syn"), "mut_mis": ("mut_n_mis", "sites_mis"),
           "mut_trunc": ("mut_n_trunc", "sites_trunc")}
bad = [c for c in a.contrasts if c not in TARGETS or c == "mut_syn"]
if bad:
    sys.exit(f"--contrasts: only mut_mis / mut_trunc (modelled on top of mut_syn); got {bad}")
XGB_PARAMS = dict(max_depth=5, eta=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=10, tree_method="hist")

# ============================ data and feature set (same rules as pilot_predictor_identity.py) ============================
D = pd.read_csv(a.dataset, sep="\t")
D["chr"] = D["binID"].str.split("_").str[0]
D["idx"] = D["binID"].str.split("_").str[1].astype(int)
chrom_order = lambda v: int(v) if str(v).isdigit() else 100 + ord(str(v)[0])
D = D.sort_values(["chr", "idx"], key=lambda s: s.map(chrom_order) if s.name == "chr" else s).reset_index(drop=True)
ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
VALIDATION_ONLY = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_"))]
FEATS = [c for c in D.columns if c not in ALT + VALIDATION_ONLY + EXPOSURE + ["cds_bp", "binID", "chr", "idx"]]
if a.gene_na_to_zero:
    gene_level = [c for c in FEATS
                  if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                           "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
                  or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs", "total_n_paralogs"))]
    for c in gene_level:
        D[c] = D[c].fillna(0)
    log(f"--gene_na_to_zero: NA set to 0 in {len(gene_level)} gene-level features")
for c in [c for c in FEATS if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)
too_many_na = [c for c in FEATS if D[c].isna().mean() >= 0.5]
constant = [c for c in FEATS if c not in too_many_na and D[c].nunique() <= 1]
FEATS = [c for c in FEATS if c not in too_many_na + constant]
log(f"{len(D)} bins, {len(FEATS)} features; dropped {len(too_many_na)} with >= 50% NA {too_many_na}, "
    f"{len(constant)} constant {constant}")
eig = [c for c in FEATS if c.startswith("eigen_")]
if eig and "GC_content" in FEATS:
    for ch, rows in D.groupby("chr").groups.items():
        for c in eig:
            r = D.loc[rows, c].corr(D.loc[rows, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0:
                D.loc[rows, c] = -D.loc[rows, c]
L = None
if a.labels:
    L = pd.read_csv(a.labels, sep="\t", dtype=str).set_index("feature")
    missing = [f for f in FEATS if f not in L.index]
    if missing:
        sys.exit(f"{a.labels}: {len(missing)} features have no row: {missing}")


def deviance_poisson(y, mu):
    mu = np.clip(mu, 1e-12, None)
    return 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0) - (y - mu))


def corr_len(y, max_lag=3000):
    """Lag (bins) at which the autocorrelation of y drops below 0.5 (buffer of the CV, as in the pilot)."""
    y = y - y.mean(); v = np.dot(y, y)
    if v == 0: return 1
    for lag in range(1, min(max_lag, len(y) - 1)):
        if np.dot(y[:-lag], y[lag:]) / v < 0.5: return lag
    return max_lag


def fit_poisson(X, y, log_offset, seed):
    params = dict(XGB_PARAMS, objective="count:poisson", nthread=a.n_jobs, seed=seed)
    return xgb.train(params, xgb.DMatrix(X, label=y, base_margin=log_offset), num_boost_round=a.n_trees)


rng = np.random.default_rng(a.random_seed)
result_rows, null_rows, fit_rows = [], [], []
for contrast in a.contrasts:
    syn_count, syn_sites = TARGETS["mut_syn"]
    con_count, con_sites = TARGETS[contrast]
    keep = ((D["cds_bp"] >= a.min_exonic_bp) & (D[syn_sites] > 0) & (D[con_sites] > 0)
            & D[syn_count].notna() & D[con_count].notna())
    d = D[keep].reset_index(drop=True)
    log(f"{contrast} on top of mut_syn: {len(D) - len(d)} of {len(D)} bins removed (cds_bp < {a.min_exonic_bp:g}, no "
        f"sites of one of the two classes or NA count); {len(d)} bins")
    n_bins = len(d)
    chroms = d["chr"].to_numpy()

    # shadow features: one per feature, circularly shifted within each chromosome
    shadow_names, shadow_cols = [], []
    for src in FEATS:
        x = d[src].to_numpy(float).copy()
        for c in np.unique(chroms):
            m = np.where(chroms == c)[0]; n = len(m)
            x[m] = np.roll(x[m], int(rng.integers(int(.1 * n), max(int(.1 * n) + 1, int(.9 * n)))))
        shadow_names.append(f"shadow_{src}"); shadow_cols.append(x)
    X = np.column_stack([d[FEATS].to_numpy(float)] + shadow_cols)
    names = FEATS + shadow_names
    is_shadow = np.array([n.startswith("shadow_") for n in names])
    k = len(names)

    y_syn_obs = d[syn_count].to_numpy(float); s_syn = d[syn_sites].to_numpy(float)
    y_con_obs = d[con_count].to_numpy(float); s_con = d[con_sites].to_numpy(float)
    rate = y_con_obs / s_con
    ell = int(np.median([corr_len(rate[chroms == c]) for c in np.unique(chroms)]))
    buffer = min(ell, a.max_buffer_bins)
    folds = list(INTERLEAVED_KFOLD_IDXS(d["binID"].to_numpy(), a.n_splits, a.block_bins, buffer,
                                        MIN_TRAIN_FRACTION=0.05).values())
    tested = np.zeros(n_bins, bool)
    for tr, te in folds:
        tested[te] = True
    if not tested.all():
        sys.exit(f"{contrast}: {int((~tested).sum())} bins are in no test fold: the importances would be wrong")
    log(f"  buffer {buffer} bins, {len(folds)} folds, {len(FEATS)} features + {len(shadow_names)} shadows, seeds {a.seeds}")

    def chained(y_syn, y_con, seed, want_fit=False):
        """Both models, out of fold. Returns imp_bg, imp_sel (mean |SHAP| over held-out bins), the held-out SHAP of the
        selection model (for the direction) and, with want_fit, the held-out deviances."""
        # 1. background model on synonymous counts; out-of-fold log-rate per site (without the offset) for every bin
        imp_bg, bg_lograte = np.zeros(k), np.full(n_bins, np.nan)
        pred_syn = np.full(n_bins, np.nan)
        for tr, te in folds:
            syn_mean = np.log(max(y_syn[tr].sum(), 1e-12) / s_syn[tr].sum())
            booster = fit_poisson(X[tr], y_syn[tr], np.log(s_syn[tr]) + syn_mean, seed)
            # margin = syn_mean only (no sites): the prediction is the log rate per site of the held-out bins
            dte = xgb.DMatrix(X[te], base_margin=np.full(len(te), syn_mean))
            bg_lograte[te] = booster.predict(dte, output_margin=True)
            imp_bg += np.abs(booster.predict(dte, pred_contribs=True)[:, :k]).sum(axis=0)
            pred_syn[te] = s_syn[te] * np.exp(bg_lograte[te])
        # 2. selection model on the contrast class, offset = log(sites) + background log-rate (+ class scale)
        imp_sel, shap_sel = np.zeros(k), np.zeros((n_bins, k))
        pred_con, pred_off = np.full(n_bins, np.nan), np.full(n_bins, np.nan)
        expected = s_con * np.exp(bg_lograte)
        for tr, te in folds:
            scale = np.log(max(y_con[tr].sum(), 1e-12) / expected[tr].sum())      # class mean on top of the background
            offset = np.log(expected) + scale
            booster = fit_poisson(X[tr], y_con[tr], offset[tr], seed)
            dte = xgb.DMatrix(X[te], base_margin=offset[te])
            phi = booster.predict(dte, pred_contribs=True)[:, :k]
            imp_sel += np.abs(phi).sum(axis=0)
            shap_sel[te] = phi
            pred_con[te] = booster.predict(dte)
            pred_off[te] = np.exp(offset[te])
        out = (imp_bg / n_bins, imp_sel / n_bins, shap_sel)
        if want_fit:
            const_syn = s_syn * y_syn.sum() / s_syn.sum()
            const_con = s_con * y_con.sum() / s_con.sum()
            fit = dict(dev_expl_bg_on_syn=1 - deviance_poisson(y_syn, pred_syn) / deviance_poisson(y_syn, const_syn),
                       dev_expl_bg_offset_on_contrast=1 - deviance_poisson(y_con, pred_off) / deviance_poisson(y_con, const_con),
                       dev_expl_sel_on_top=1 - deviance_poisson(y_con, pred_con) / deviance_poisson(y_con, pred_off))
            out = out + (fit,)
        return out

    # observed data, every seed
    imp_bg_s, imp_sel_s, shap_sum = [], [], np.zeros((n_bins, k))
    for seed in a.seeds:
        imp_bg, imp_sel, shap_sel, fit = chained(y_syn_obs, y_con_obs, seed, want_fit=True)
        imp_bg_s.append(imp_bg); imp_sel_s.append(imp_sel); shap_sum += shap_sel
        fit_rows.append(dict(contrast=contrast, seed=seed, n_bins=n_bins, **fit))
        log(f"  seed {seed}: deviance explained: background on synonymous {fit['dev_expl_bg_on_syn']:.4f}, background "
            f"offset on {contrast} {fit['dev_expl_bg_offset_on_contrast']:.4f}, selection model on top "
            f"{fit['dev_expl_sel_on_top']:.4f}")
    imp_bg_s, imp_sel_s = np.array(imp_bg_s), np.array(imp_sel_s)
    imp_bg, imp_sel = imp_bg_s.mean(axis=0), imp_sel_s.mean(axis=0)
    shap_mean = shap_sum / len(a.seeds)

    # null of the selection model: parametric bootstrap from a background model fitted on ALL synonymous counts
    syn_mean_all = np.log(y_syn_obs.sum() / s_syn.sum())
    gen = fit_poisson(X, y_syn_obs, np.log(s_syn) + syn_mean_all, a.seeds[0])
    gen_rate = np.exp(gen.predict(xgb.DMatrix(X, base_margin=np.full(n_bins, syn_mean_all)), output_margin=True))
    mu_syn_null = s_syn * gen_rate
    mu_con_null = s_con * gen_rate * y_con_obs.sum() / (s_con * gen_rate).sum()
    log(f"  null generator (background only, no selection): expected totals syn {mu_syn_null.sum():.0f} (observed "
        f"{y_syn_obs.sum():.0f}), {contrast} {mu_con_null.sum():.0f} (observed {y_con_obs.sum():.0f})")
    imp_sel_null = []
    for r in range(a.n_null):
        y_syn_sim = rng.poisson(mu_syn_null).astype(float)
        y_con_sim = rng.poisson(mu_con_null).astype(float)
        imp_sel_r = np.mean([chained(y_syn_sim, y_con_sim, seed)[1] for seed in a.seeds], axis=0)
        imp_sel_null.append(imp_sel_r)
        for j, name in enumerate(names):
            null_rows.append(dict(contrast=contrast, feature=name, null_rep=r, imp_sel=imp_sel_r[j]))
        log(f"  null {r + 1}/{a.n_null} done")
    imp_sel_null = np.array(imp_sel_null)

    # thresholds and labels
    thr_bg = imp_bg[is_shadow].max()
    thr_sel = imp_sel[is_shadow].max()
    p_sel = (1 + (imp_sel_null >= imp_sel).sum(axis=0)) / (1 + a.n_null)
    z_sel = (imp_sel - imp_sel_null.mean(axis=0)) / imp_sel_null.std(axis=0, ddof=1)
    for j, name in enumerate(names):
        bg_sig = bool(imp_bg[j] > thr_bg)
        sel_sig = bool(imp_sel[j] > thr_sel and p_sel[j] < a.alpha)
        label = {(True, False): "O", (False, True): "S", (True, True): "both", (False, False): "none"}[(bg_sig, sel_sig)]
        xj = X[:, j]; ok = ~np.isnan(xj)
        direction = spearmanr(xj[ok], shap_mean[ok, j])[0] if ok.sum() > 2 else np.nan
        prior = (L.loc[name, contrast].strip().upper() if (L is not None and name in L.index) else "")
        result_rows.append(dict(contrast=contrast, feature=name, shadow=bool(is_shadow[j]),
                                imp_bg=imp_bg[j], imp_bg_sd=imp_bg_s[:, j].std(ddof=1), thr_bg=thr_bg, bg_sig=bg_sig,
                                imp_sel=imp_sel[j], imp_sel_sd=imp_sel_s[:, j].std(ddof=1), thr_sel=thr_sel,
                                imp_sel_null_mean=imp_sel_null[:, j].mean(), p_sel=p_sel[j], z_sel=z_sel[j],
                                sel_sig=sel_sig, direction=direction, label=label, prior_label=prior))
    R = pd.DataFrame([r for r in result_rows if r["contrast"] == contrast])
    real, shad = R[~R["shadow"]], R[R["shadow"]]
    log(f"  {contrast}: labels {real['label'].value_counts().to_dict()}; shadows with p_sel < {a.alpha}: "
        f"{int((shad['p_sel'] < a.alpha).sum())} of {len(shad)} (shadows can never pass the max-of-shadows rule)")
    if L is not None:
        T = pd.crosstab(real["prior_label"], real["label"])
        T.to_csv(f"{a.out_dir}/label_vs_prior_{contrast}.tsv", sep="\t")
        log(f"  {contrast}: new label (columns) vs prior label (rows):\n{T.to_string()}")

pd.DataFrame(result_rows).to_csv(f"{a.out_dir}/chained_identity.tsv", sep="\t", index=False)
pd.DataFrame(null_rows).to_csv(f"{a.out_dir}/chained_null.tsv", sep="\t", index=False)
pd.DataFrame(fit_rows).to_csv(f"{a.out_dir}/model_fit.tsv", sep="\t", index=False)
log(f"written chained_identity.tsv, chained_null.tsv, model_fit.tsv in {a.out_dir}")
