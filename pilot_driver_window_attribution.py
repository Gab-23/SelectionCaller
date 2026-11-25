"""
BURNER - feature attribution on driver bins vs flanking windows of increasing size (Gab, 25/11/2026).

Question: which features carry SELECTION? Selection acts on the functional elements (here: bins overlapping the
cohort's IntOGen drivers), supply acts everywhere. If we train a count model on driver bins + matched non-driver bins and
then let the window around each bin grow (target and exposure summed over the centre bin +- w bins), the events of the
neighbours dilute the driver excess. Features whose importance DROPS with w are candidate selection features; supply
features should keep (or gain) importance.

Design (decided by Gab, 25/11/2026):
  * units      : driver bins = bins with dist.to.closest.OG == 0 or dist.to.closest.TSG == 0 (cohort-specific IntOGen
                 lists of the dataset) and coding sites > 0; controls = --n_controls non-driver bins per driver bin,
                 matched ONLY on exposure (total coding sites, log scale, nearest neighbour without replacement).
                 Controls are > --max_window bins from every driver bin (their windows never contain driver events) and
                 > 2 * --max_window bins from each other (no shared events between control windows).
  * windows    : w = 0..--max_window bins on each side (20 kb bins: 0 to +-200 kb). Bins missing from the matrix
                 (blacklist) are skipped; windows are truncated at chromosome ends.
  * features   : option 1 = the features of the CENTRE bin at every w (only target and exposure grow), so a change of
                 importance comes from the target, not from smoothing the features.
  * model      : XGBoost Poisson, offset log(exposure of the class in the window), start = training mean rate, same
                 hyperparameters as pilot_omics_background.py (depth 5, eta 0.05, 300 rounds, subsample/colsample 0.8,
                 min_child_weight 10).
  * folds      : over DRIVER CLUSTERS (driver bins closer than 2 * --max_window + 1 bins are one cluster, so windows of
                 the same region never fall in training and test); each control follows the fold of its driver bin.
                 Importance = mean |SHAP| (log-rate scale) on the held-out units; deviance explained out-of-fold.
  * classes    : mis, trunc and syn. syn is the CONTROL class: selection on synonymous changes is ~absent, so a feature
                 whose importance drops in syn too drops for a reason other than selection (dilution of supply signal,
                 exposure shift).
  * placebo    : --n_placebo random N(0,1) features per unit (same values at every w): importance floor per window.
  * excluded   : alteration columns, mut_n_*, sites_*, cds_bp, and dist.to.closest.OG / TSG (they define the driver bins:
                 circular). Same preprocessing as pilot_class_interaction_test.py (gene-level NA -> 0, distance NA ->
                 2 Mb, drop >= 50% NA and constant features, Hi-C eigenvectors oriented by GC per chromosome).

Limits (known): a long driver has several driver bins, and the window of one of them contains the others (less
dilution for long drivers); controls are matched on exposure only, so driver and control units differ in every other
feature (intended: that contrast is what the model can learn).

Inputs : --dataset feature matrix (TSV: binID, features, mut_n_syn/mis/trunc, sites_syn/mis/trunc, dist.to.closest.OG/TSG)
Outputs (--out_dir):
  units.tsv        unit binID, role (driver / control), cluster, fold, matched driver bin, total sites
  importance.tsv   class, window, fold, feature, mean_abs_shap, share (of the total mean |SHAP| of that model)
  performance.tsv  class, window, fold, n_units, n_driver_units, dev_expl_oof; plus the raw rate ratio driver / control
                   (events / exposure) per class and window (fold = 'all'): the dilution of the driver excess itself
  summary.tsv      class, feature, median importance at w = 0 and w = max over folds, ratio max / 0, Spearman of the
                   median importance vs w, placebo max at w = 0 and w = max
Usage (env py_ML / pilot env, from the repo root):
  python pilot_driver_window_attribution.py --dataset <20 kb feature matrix> --out_dir analyses/v4/BRCA/20000bp/driver_window_attribution
"""
import argparse, os, time
import numpy as np, pandas as pd
import xgboost as xgb
from scipy.stats import spearmanr

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--max_window", type=int, default=10, help="largest window, in bins on each side of the centre bin")
p.add_argument("--n_controls", type=int, default=3, help="matched non-driver bins per driver bin")
p.add_argument("--n_folds", type=int, default=5)
p.add_argument("--n_placebo", type=int, default=5)
p.add_argument("--classes", nargs="+", default=["mis", "trunc", "syn"])
p.add_argument("--n_jobs", type=int, default=4)
p.add_argument("--seed", type=int, default=1)
# NOTE (25/11/2026): toy test showed that with the centre-bin features (option 1) a SUPPLY feature with short-range
# variation loses importance as fast as a selection feature (the neighbours' events follow the neighbours' own supply,
# which the centre bin does not see). window_mean = features averaged over the same window as the target (weighted by
# the coding sites of each bin): supply stays aligned with the window target, only the centre-concentrated selection
# dilutes. Kept as an option so the two can be compared; default = Gab's choice (centre).
p.add_argument("--feature_mode", choices=["centre", "window_mean"], default="centre",
               help="centre = features of the centre bin at every w (option 1); window_mean = coding-site-weighted mean of "
                    "each feature over the window (bins without coding sites get weight 0; NA ignored)")
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

needed = ["mut_n_syn", "mut_n_mis", "mut_n_trunc", "sites_syn", "sites_mis", "sites_trunc",
          "dist.to.closest.OG", "dist.to.closest.TSG"]
missing = [c for c in needed if c not in D.columns]
if missing:
    raise SystemExit(f"ERROR: {a.dataset} lacks the columns {missing}: use a feature matrix with COSMIC classes and IntOGen distances.")

ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
DRIVER_DEF = ["dist.to.closest.OG", "dist.to.closest.TSG"]          # define the driver bins: never features
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_"))]
EXONIC = ["cds_bp"] if "cds_bp" in D.columns else []
FEATS = [c for c in D.columns if c not in ALT + DRIVER_DEF + EXPOSURE + EXONIC + ["binID", "chr", "idx"]]
GENE_LEVEL_NA_TO_ZERO = [c for c in FEATS
                         if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                                  "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
                         or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs",
                                          "total_n_paralogs"))]
for c in GENE_LEVEL_NA_TO_ZERO:
    D[c] = D[c].fillna(0)
for c in [c for c in FEATS + DRIVER_DEF if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)                                   # NA on a whole chromosome = no gene of that list: cap
too_many_na = [c for c in FEATS if D[c].isna().mean() >= 0.5]
constant = [c for c in FEATS if c not in too_many_na and D[c].nunique() <= 1]
FEATS = [c for c in FEATS if c not in too_many_na + constant]
log(f"{len(D)} bins, {len(FEATS)} features; dropped: {len(too_many_na)} with >= 50% NA {too_many_na}, "
    f"{len(constant)} constant {constant}")
eig = [c for c in FEATS if c.startswith("eigen_")]
if eig and "GC_content" in FEATS:
    for ch, rows in D.groupby("chr").groups.items():
        for c in eig:
            r = D.loc[rows, c].corr(D.loc[rows, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0:
                D.loc[rows, c] = -D.loc[rows, c]

# ============================ driver bins, clusters, folds ============================
D["sites_total"] = D["sites_syn"] + D["sites_mis"] + D["sites_trunc"]
coding = D["sites_total"] > 0
is_driver_body = (D["dist.to.closest.OG"] == 0) | (D["dist.to.closest.TSG"] == 0)
driver_rows = np.flatnonzero(is_driver_body & coding)
log(f"driver bins: {int(is_driver_body.sum())} overlap an OG/TSG gene body, {len(driver_rows)} of them with coding sites "
    f"(the others are dropped: no coding events possible)")
if len(driver_rows) < 5 * a.n_folds:
    raise SystemExit(f"ERROR: only {len(driver_rows)} driver bins with coding sites: check that dist.to.closest.OG/TSG are 0 "
                     f"for bins overlapping a driver in {a.dataset}")

# clusters: driver bins of the same chromosome closer than 2 * max_window + 1 bins (their windows can overlap)
gap = 2 * a.max_window + 1
cluster = np.full(len(D), -1)
cid = -1; prev_chr, prev_idx = None, None
for r in driver_rows:
    if D.at[r, "chr"] != prev_chr or D.at[r, "idx"] - prev_idx > gap:
        cid += 1
    cluster[r] = cid
    prev_chr, prev_idx = D.at[r, "chr"], D.at[r, "idx"]
n_clusters = cid + 1
cluster_fold = rng.permutation(np.arange(n_clusters) % a.n_folds)   # shuffled round-robin: balanced number of clusters
fold = np.full(len(D), -1)
fold[driver_rows] = cluster_fold[cluster[driver_rows]]
log(f"{n_clusters} driver clusters; driver bins per fold: {np.bincount(fold[driver_rows], minlength=a.n_folds).tolist()}")

# ============================ matched controls ============================
# pool: coding non-driver bins farther than max_window bins from every driver bin (same chromosome)
far_from_driver = np.ones(len(D), bool)
for ch, rows in D.groupby("chr").groups.items():
    rows = np.asarray(rows)
    drv_idx = D.loc[rows, "idx"].to_numpy()[is_driver_body.to_numpy()[rows]]
    if len(drv_idx) == 0:
        continue
    idx = D.loc[rows, "idx"].to_numpy()
    pos = np.searchsorted(drv_idx, idx)                       # nearest driver index on both sides
    left = np.abs(idx - drv_idx[np.clip(pos - 1, 0, len(drv_idx) - 1)])
    right = np.abs(drv_idx[np.clip(pos, 0, len(drv_idx) - 1)] - idx)
    far_from_driver[rows] = np.minimum(left, right) > a.max_window
pool = np.flatnonzero(coding.to_numpy() & ~is_driver_body.to_numpy() & far_from_driver)
log(f"control pool: {len(pool)} coding non-driver bins > {a.max_window} bins from any driver bin")

log_sites = np.log(D["sites_total"].to_numpy(float).clip(min=1))
chr_arr, idx_arr = D["chr"].to_numpy(), D["idx"].to_numpy()
available = np.ones(len(pool), bool)
pool_logs = log_sites[pool]
matched_to = np.full(len(D), -1)
for r in rng.permutation(driver_rows):                       # random order: no driver gets systematically the best match
    for _ in range(a.n_controls):
        cand = np.flatnonzero(available)
        if len(cand) == 0:
            raise SystemExit("ERROR: control pool exhausted: lower --n_controls or --max_window")
        best = cand[np.argmin(np.abs(pool_logs[cand] - log_sites[r]))]
        c_row = pool[best]
        matched_to[c_row] = r; fold[c_row] = fold[r]
        # exclusion zone: no other control within 2 * max_window bins (control windows never share bins)
        close = (chr_arr[pool] == chr_arr[c_row]) & (np.abs(idx_arr[pool] - idx_arr[c_row]) <= 2 * a.max_window)
        available[close] = False
control_rows = np.flatnonzero(matched_to >= 0)
match_logratio = log_sites[control_rows] - log_sites[matched_to[control_rows]]
log(f"{len(control_rows)} controls; matching quality |log(sites control / sites driver)|: median "
    f"{np.median(np.abs(match_logratio)):.3f}, 95th pct {np.quantile(np.abs(match_logratio), 0.95):.3f}")

units = np.concatenate([driver_rows, control_rows])
is_driver_unit = np.r_[np.ones(len(driver_rows), bool), np.zeros(len(control_rows), bool)]
unit_fold = fold[units]
pd.DataFrame({"binID": D["binID"].to_numpy()[units], "role": np.where(is_driver_unit, "driver", "control"),
              # controls carry the cluster of their matched driver bin
              "cluster": np.where(is_driver_unit, cluster[units], cluster[np.maximum(matched_to[units], 0)]),
              "fold": unit_fold,
              "matched_driver": np.where(is_driver_unit, "", D["binID"].to_numpy()[np.maximum(matched_to[units], 0)]),
              "sites_total": D["sites_total"].to_numpy()[units]}).to_csv(f"{a.out_dir}/units.tsv", sep="\t", index=False)

# ============================ window sums (target and exposure) ============================
def window_sums(col, w):
    """Sum of column col over the bins centre-w .. centre+w of the same chromosome, for every unit. Bins absent from the
    matrix (blacklist) contribute nothing."""
    out = np.zeros(len(units))
    for ch in pd.unique(chr_arr[units]):
        rows = np.flatnonzero(chr_arr == ch)
        idx = idx_arr[rows]
        full = np.zeros(idx.max() + 2 * w + 2)                 # dense vector over bin numbers (shifted by w)
        full[idx + w] = D[col].to_numpy(float)[rows]
        csum = np.r_[0, np.cumsum(full)]
        sel = np.flatnonzero(chr_arr[units] == ch)
        centre = idx_arr[units[sel]] + w
        out[sel] = csum[centre + w + 1] - csum[centre - w]
    return out

X_placebo = rng.standard_normal((len(units), a.n_placebo))
names = FEATS + [f"placebo_{i + 1}" for i in range(a.n_placebo)]

def window_features(w):
    """Feature matrix of the units at window w: centre-bin features, or (window_mean) the coding-site-weighted mean over
    the window; NA values get weight 0 (if a whole window is NA the feature stays NA, XGBoost handles it)."""
    if a.feature_mode == "centre" or w == 0:
        return np.column_stack([D[FEATS].to_numpy(float)[units], X_placebo])
    weight = D["sites_total"].to_numpy(float)
    cols = []
    for f in FEATS:
        v = D[f].to_numpy(float)
        ok = np.isfinite(v)
        D["_wv"] = np.where(ok, v * weight, 0.0); D["_w"] = np.where(ok, weight, 0.0)
        num, den = window_sums("_wv", w), window_sums("_w", w)
        cols.append(np.where(den > 0, num / np.where(den > 0, den, 1), np.nan))
    return np.column_stack(cols + [X_placebo])
log(f"feature mode: {a.feature_mode}")

params = dict(objective="count:poisson", max_depth=5, eta=0.05, subsample=0.8, colsample_bytree=0.8,
              min_child_weight=10, tree_method="hist", nthread=a.n_jobs, seed=a.seed)

def poisson_dev(y, mu):
    mu = np.maximum(mu, 1e-12)
    term = np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)
    return 2 * np.sum(term - (y - mu))

imp_rows, perf_rows = [], []
for cls in a.classes:
    for w in range(a.max_window + 1):
        y_all = window_sums(f"mut_n_{cls}", w)
        e_all = window_sums(f"sites_{cls}", w)
        X = window_features(w)
        keep = e_all > 0
        if (~keep).any():
            log(f"{cls} w={w}: {int((~keep).sum())} units without {cls} sites in the window dropped")
        rate_drv = y_all[keep & is_driver_unit].sum() / e_all[keep & is_driver_unit].sum()
        rate_ctl = y_all[keep & ~is_driver_unit].sum() / e_all[keep & ~is_driver_unit].sum()
        perf_rows.append(dict(**{"class": cls}, window=w, fold="all", n_units=int(keep.sum()),
                              n_driver_units=int((keep & is_driver_unit).sum()), dev_expl_oof=np.nan,
                              rate_ratio_driver_vs_control=rate_drv / rate_ctl))
        for k in range(a.n_folds):
            tr = np.flatnonzero(keep & (unit_fold != k)); te = np.flatnonzero(keep & (unit_fold == k))
            off = np.log(e_all)
            c_tr = np.log(max(y_all[tr].sum(), 1e-12) / e_all[tr].sum())     # training mean rate (no test units)
            booster = xgb.train(params, xgb.DMatrix(X[tr], label=y_all[tr], base_margin=off[tr] + c_tr), num_boost_round=300)
            dte = xgb.DMatrix(X[te], base_margin=off[te] + c_tr)
            mu = booster.predict(dte)
            null = e_all[te] * np.exp(c_tr)
            dev_expl = 1 - poisson_dev(y_all[te], mu) / poisson_dev(y_all[te], null)
            shap = booster.predict(dte, pred_contribs=True)[:, :-1]           # last column = bias
            mean_abs = np.abs(shap).mean(axis=0)
            total = mean_abs.sum()
            for f, v in zip(names, mean_abs):
                imp_rows.append(dict(**{"class": cls}, window=w, fold=k, feature=f, mean_abs_shap=v,
                                     share=v / total if total > 0 else np.nan))
            perf_rows.append(dict(**{"class": cls}, window=w, fold=k, n_units=len(te),
                                  n_driver_units=int(is_driver_unit[te].sum()), dev_expl_oof=dev_expl,
                                  rate_ratio_driver_vs_control=np.nan))
        log(f"{cls} w={w}: rate ratio driver/control {rate_drv / rate_ctl:.2f}, OOF dev. explained "
            f"{np.mean([r['dev_expl_oof'] for r in perf_rows[-a.n_folds:]]):.3f} (mean over folds)")

IMP = pd.DataFrame(imp_rows)
PERF = pd.DataFrame(perf_rows)
IMP.to_csv(f"{a.out_dir}/importance.tsv", sep="\t", index=False)
PERF.to_csv(f"{a.out_dir}/performance.tsv", sep="\t", index=False)

# ============================ summary: change of importance with the window ============================
med = IMP.groupby(["class", "feature", "window"])["mean_abs_shap"].median().reset_index()
sum_rows = []
for (cls, f), g in med.groupby(["class", "feature"]):
    g = g.sort_values("window")
    i0, imax = g["mean_abs_shap"].iloc[0], g["mean_abs_shap"].iloc[-1]
    rho = spearmanr(g["window"], g["mean_abs_shap"]).correlation if g["mean_abs_shap"].nunique() > 1 else np.nan
    sum_rows.append(dict(**{"class": cls}, feature=f, imp_w0=i0, imp_wmax=imax,
                         ratio_wmax_w0=imax / i0 if i0 > 0 else np.nan, spearman_vs_window=rho))
S = pd.DataFrame(sum_rows)
plc = med[med["feature"].str.startswith("placebo_")]
for cls in S["class"].unique():
    pc = plc[plc["class"] == cls]
    S.loc[S["class"] == cls, "placebo_max_w0"] = pc.loc[pc["window"] == 0, "mean_abs_shap"].max()
    S.loc[S["class"] == cls, "placebo_max_wmax"] = pc.loc[pc["window"] == a.max_window, "mean_abs_shap"].max()
S = S.sort_values(["class", "imp_w0"], ascending=[True, False])
S.to_csv(f"{a.out_dir}/summary.tsv", sep="\t", index=False)
log(f"done: {a.out_dir}/units.tsv, importance.tsv, performance.tsv, summary.tsv")
