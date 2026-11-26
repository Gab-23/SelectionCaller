#!/usr/bin/env python
"""
BURNER (26/11/2026) - potential of a background trained ONLY on bins without functional elements, tested on the elements.
ATAC logFC (tumour vs normal), BRCA, one resolution. Gab's idea: the non-functional genome represents the "random"
part of an alteration; a model of it, applied to functional elements, should leave selection in the residuals.

Bins
  neutral : no protein-coding gene (annotated backbone), no lncRNA gene (num.overlapping.genes_lncRNA == 0),
            0 bp of ChromHMM enhancer states (10-18) and 0 bp of promoter states (1-4, 22, 23).
  element : every other bin.
Arms (same XGBoost as pilot_omics_background.py, Gaussian; same interleaved folds with buffer for every arm, so every
bin gets an out-of-fold prediction in every arm)
  A : trained on the training bins of the fold (all bins)          -> current approach
  B : trained on the NEUTRAL training bins of the fold only        -> Gab's idea
  each with feature set FULL (pilot features) and REDUCED (FULL minus the features that define the split or are
  undefined without genes: enhancer/promoter states, gene overlaps, gene distances, gene-level scores). In arm B the
  FULL set contains features that are ~constant in training (extrapolation by construction); REDUCED avoids that.
Readouts
  1. calibration.tsv : per arm, on neutral and on element bins separately: mean residual (obs - pred), SD, R2,
                       calibration slope (obs ~ pred), Spearman.
  2. adversarial.tsv : out-of-fold AUC of an XGBoost classifier neutral vs element (FULL and REDUCED features): how
                       different the element bins are from the training bins of arm B (0.5 = indistinguishable).
     support.tsv     : per feature, share of element bins outside the 1st-99th percentile of the neutral bins.
  3. contrast.tsv    : element bins of protein-coding genes; driver = holds a gene of the IntOGen pan-cancer compendium,
                       non-driver = holds only other protein-coding genes. Each driver bin is matched to its --k_match
                       nearest non-driver bins (rank-normal GC_content, replication.timing.smoothed, eigen_1,
                       num.overlapping.genes_protein_coding). Per arm: signed residual difference (driver - matched;
                       meaningful per role, Act and LoF may have opposite signs and cancel in "all"), absolute deviation
                       from the local matched baseline and its AUC (direction-free; see contrast_stats for why raw
                       |residual| is not used). 95% CI from a bootstrap over
                       runs of consecutive driver bins; B - A differences with the same resamples (paired).
                       Groups: all drivers, Act (oncogene role), LoF (tumour suppressor role) by the compendium ROLE.
  predictions.tsv.gz : binID, neutral, driver, role, logFC, out-of-fold prediction of every arm.
Limits: no neutral readout INSIDE elements for ATAC (unlike synonymous mutations), so a residual shift of all element bins
cannot be split into selection vs poor transfer: only the driver vs matched contrast cancels a shared transfer bias.
Tumour vs normal logFC is dominated by cell composition (epithelial tumour vs adipose/stromal normal breast).
Usage (repository root, pilot conda env):
  python pilot_atac_neutral_training.py --dataset <BRCA feature matrix> --out_dir analyses/v4/BRCA/20000bp/atac_neutral_training
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import xgboost as xgb
from scipy.stats import norm, spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, os.getcwd())
try:
    from workflow.scripts.helpers.INTERLEAVED_KFOLD_IDXS import INTERLEAVED_KFOLD_IDXS   # same folds as the pilot
except ImportError:
    sys.exit("Needs workflow/scripts/helpers/INTERLEAVED_KFOLD_IDXS.py: run from the repository root.")

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--backbone_annotated", default="resources/backbone/hg38_binned_20000_0_v45.filtered.sorted.annotated.bed")
p.add_argument("--compendium", default="resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv")
p.add_argument("--out_dir", required=True)
p.add_argument("--n_splits", type=int, default=5)
p.add_argument("--block_bins", type=int, default=100)
p.add_argument("--buffer_bins", type=int, default=1, help="pilot v4 used 1 bin for ATAC (correlation length l = 1)")
p.add_argument("--k_match", type=int, default=5)
p.add_argument("--n_boot", type=int, default=1000)
p.add_argument("--n_jobs", type=int, default=8)
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
rng = np.random.default_rng(a.seed)
t0 = time.time()
LOG = []
def log(msg):
    print(msg, flush=True); LOG.append(str(msg))

# ---------------- data, same preprocessing as pilot_omics_background.py (genome_nopos run of v4) ----------------
D = pd.read_csv(a.dataset, sep="\t", dtype={"binID": str})
if "logFC" not in D.columns:
    sys.exit(f"{a.dataset}: no logFC column (ATAC is available for BRCA only).")
D["chr"] = D["binID"].str.split("_").str[0]
D["idx"] = D["binID"].str.split("_").str[1].astype(int)
chrom_order = lambda s: s.map(lambda v: int(v) if str(v).isdigit() else v) if s.name == "chr" else s
D = D.sort_values(["chr", "idx"], key=chrom_order).reset_index(drop=True)
n_na_target = int(D["logFC"].isna().sum())
if n_na_target:
    log(f"Dropped {n_na_target} bins with NA logFC")
    D = D[D["logFC"].notna()].reset_index(drop=True)

ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
VALIDATION_ONLY = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_"))]
POSITION = ["dist_to_centromere", "dist_to_telomere"]          # excluded in the v4 genome_nopos run
FEATS = [c for c in D.columns if c not in ALT + VALIDATION_ONLY + EXPOSURE + POSITION + ["cds_bp", "binID", "chr", "idx"]]
for c in [c for c in FEATS if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)                                    # NA on whole chromosomes = no element of that list: cap
dropped_na = [c for c in FEATS if D[c].isna().mean() >= 0.5]
FEATS = [c for c in FEATS if c not in dropped_na]
log(f"{len(D)} bins, {len(FEATS)} features; dropped (>= 50% NA): {dropped_na}")
for c in [c for c in FEATS if c.startswith("eigen_")]:          # Hi-C eigenvector sign oriented by GC per chromosome
    for ch, rows in D.groupby("chr").groups.items():
        r = D.loc[rows, c].corr(D.loc[rows, "GC_content"], method="spearman")
        if pd.notna(r) and r < 0:
            D.loc[rows, c] = -D.loc[rows, c]

# ---------------- neutral vs element bins ----------------
ENH_STATES = ["10_TxEnh5'", "11_TxEnh3'", "12_TxEnhW", "13_EnhA1", "14_EnhA2", "15_EnhAF", "16_EnhW1", "17_EnhW2", "18_EnhAc"]
PROM_STATES = ["1_TssA", "2_PromU", "3_PromD1", "4_PromD2", "22_PromP", "23_PromBiv"]
needed = ENH_STATES + PROM_STATES + ["num.overlapping.genes_lncRNA", "num.overlapping.genes_protein_coding",
                                     "GC_content", "replication.timing.smoothed", "eigen_1"]
missing = [c for c in needed if c not in D.columns]
if missing:
    sys.exit(f"{a.dataset}: missing columns needed to define neutral bins / matching: {missing}")
B = pd.read_csv(a.backbone_annotated, sep="\t", header=None, names=["chrom", "start", "end", "binID", "genes"], dtype=str)
genes_of_bin = dict(zip(B["binID"], B["genes"].str.split(",")))
has_pc_gene = D["binID"].map(lambda b: b in genes_of_bin).to_numpy()
n_disagree = int((has_pc_gene != (D["num.overlapping.genes_protein_coding"].fillna(0) > 0).to_numpy()).sum())
log(f"Bins with a protein-coding gene: annotated backbone {int(has_pc_gene.sum())}; "
    f"disagreement with num.overlapping.genes_protein_coding > 0: {n_disagree} bins")
enh_bp = D[ENH_STATES].fillna(0).sum(axis=1).to_numpy()
prom_bp = D[PROM_STATES].fillna(0).sum(axis=1).to_numpy()
has_lnc = (D["num.overlapping.genes_lncRNA"].fillna(0) > 0).to_numpy()
neutral = ~has_pc_gene & ~has_lnc & (enh_bp == 0) & (prom_bp == 0)
log(f"Neutral bins: {int(neutral.sum())} ({100 * neutral.mean():.1f}%); element bins: {int((~neutral).sum())} "
    f"(pc gene {int(has_pc_gene.sum())}, lncRNA {int(has_lnc.sum())}, enhancer > 0 bp {int((enh_bp > 0).sum())}, "
    f"promoter > 0 bp {int((prom_bp > 0).sum())}; categories overlap)")

# ---------------- driver bins (IntOGen pan-cancer compendium) ----------------
C = pd.read_csv(a.compendium, sep="\t")
# one role per gene: the most frequent ROLE over the gene's rows (cohorts); ties -> "ambiguous"
role_of_gene = {}
for gene, roles in C.groupby("SYMBOL")["ROLE"]:
    counts = roles.fillna("ambiguous").value_counts()
    role_of_gene[gene] = counts.index[0] if (len(counts) == 1 or counts.iloc[0] > counts.iloc[1]) else "ambiguous"
log(f"Compendium: {len(role_of_gene)} genes; roles {pd.Series(role_of_gene).value_counts().to_dict()}")
driver_genes_of_bin = [[g for g in genes_of_bin.get(b, []) if g in role_of_gene] for b in D["binID"]]
driver = np.array([len(g) > 0 for g in driver_genes_of_bin])
def bin_role(genes):
    roles = {role_of_gene[g] for g in genes}
    return roles.pop() if len(roles) == 1 else ("mixed" if roles else "")
role = np.array([bin_role(g) for g in driver_genes_of_bin])
genes_found = {g for gs in driver_genes_of_bin for g in gs}
log(f"Driver bins: {int(driver.sum())} (genes found in bins: {len(genes_found)} of {len(role_of_gene)}); "
    f"neutral driver bins (must be 0): {int((driver & neutral).sum())}; by role {pd.Series(role[driver]).value_counts().to_dict()}")

# ---------------- feature sets ----------------
GENE_LEVEL = [c for c in FEATS if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                                        "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
              or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs", "total_n_paralogs"))]
SPLIT_DEFINING = [c for c in FEATS if c in ENH_STATES + PROM_STATES or c.startswith("num.overlapping.genes_")
                  or c in ("dist.to.closest.gene_protein_coding", "dist.to.closest.gene_lncRNA")]
REDUCED = [c for c in FEATS if c not in SPLIT_DEFINING + GENE_LEVEL]
log(f"REDUCED = FULL ({len(FEATS)}) minus split-defining {SPLIT_DEFINING} and gene-level {GENE_LEVEL} -> {len(REDUCED)} features")
FEATURE_SETS = {"FULL": FEATS, "REDUCED": REDUCED}

# ---------------- out-of-fold fits ----------------
y = D["logFC"].to_numpy(float)
folds = INTERLEAVED_KFOLD_IDXS(D["binID"].to_numpy(), a.n_splits, a.block_bins, a.buffer_bins, MIN_TRAIN_FRACTION=0.05)
def xgb_regressor():
    return xgb.XGBRegressor(n_estimators=300, max_depth=5, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                            min_child_weight=10, tree_method="hist", n_jobs=a.n_jobs, random_state=a.seed)
PRED = {}
for fs_name, cols in FEATURE_SETS.items():
    X = D[cols].to_numpy(float)
    for arm in ("A", "B"):
        pred = np.full(len(D), np.nan)
        for k, (tr, te) in folds.items():
            train = tr if arm == "A" else tr[neutral[tr]]          # arm B: neutral training bins of the fold only
            model = xgb_regressor().fit(X[train], y[train])
            pred[te] = model.predict(X[te])
        PRED[f"{arm}_{fs_name}"] = pred
        log(f"[{time.time() - t0:.0f}s] arm {arm} {fs_name}: done")

# ---------------- 1. calibration ----------------
rows = []
for name, pred in PRED.items():
    for set_name, mask in (("neutral", neutral), ("element", ~neutral)):
        obs, prd = y[mask], pred[mask]
        res = obs - prd
        slope = np.polyfit(prd, obs, 1)[0]                     # 1 = calibrated spread; < 1 over-dispersed predictions
        rows.append(dict(arm=name, bins=set_name, n=int(mask.sum()), mean_residual=res.mean(), sd_residual=res.std(),
                         R2=1 - (res ** 2).sum() / ((obs - obs.mean()) ** 2).sum(), calibration_slope=slope,
                         spearman=spearmanr(obs, prd).correlation))
CAL = pd.DataFrame(rows)
CAL.to_csv(f"{a.out_dir}/calibration.tsv", sep="\t", index=False)
log("\n1. CALIBRATION (residual = obs - pred)\n" + CAL.round(4).to_string(index=False))

# ---------------- 2. how different are element bins from neutral bins ----------------
rows = []
for fs_name, cols in FEATURE_SETS.items():
    X = D[cols].to_numpy(float)
    prob = np.full(len(D), np.nan)
    for k, (tr, te) in folds.items():
        clf = xgb.XGBClassifier(n_estimators=300, max_depth=5, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                                min_child_weight=10, tree_method="hist", n_jobs=a.n_jobs, random_state=a.seed)
        clf.fit(X[tr], (~neutral[tr]).astype(int))
        prob[te] = clf.predict_proba(X[te])[:, 1]
    rows.append(dict(features=fs_name, n_features=len(cols), auc_neutral_vs_element=roc_auc_score(~neutral, prob)))
ADV = pd.DataFrame(rows)
ADV.to_csv(f"{a.out_dir}/adversarial.tsv", sep="\t", index=False)
rows = []
for c in FEATS:
    ref = D.loc[neutral, c].dropna()
    el = D.loc[~neutral, c]
    lo, hi = ref.quantile([0.01, 0.99]) if len(ref) else (np.nan, np.nan)
    rows.append(dict(feature=c, in_reduced=c in REDUCED, na_frac_neutral=D.loc[neutral, c].isna().mean(),
                     na_frac_element=el.isna().mean(), n_unique_neutral=ref.nunique(),
                     frac_element_outside_neutral_1_99=((el < lo) | (el > hi)).mean()))
SUP = pd.DataFrame(rows).sort_values("frac_element_outside_neutral_1_99", ascending=False)
SUP.to_csv(f"{a.out_dir}/support.tsv", sep="\t", index=False)
log("\n2. ADVERSARIAL AUC (neutral vs element)\n" + ADV.round(3).to_string(index=False))
log("Features most outside the neutral range (top 15):\n" + SUP.head(15).round(3).to_string(index=False))

# ---------------- 3. driver vs matched non-driver contrast (element bins of protein-coding genes) ----------------
def rank_normal(x):
    r = pd.Series(x).rank(method="average").to_numpy()
    return norm.ppf((r - 0.5) / len(r))
MATCH_COLS = ["GC_content", "replication.timing.smoothed", "eigen_1", "num.overlapping.genes_protein_coding"]
pool = np.where(~neutral & has_pc_gene)[0]                     # drivers are always here (they hold a pc gene)
Z = np.column_stack([rank_normal(D[c].fillna(D[c].median()).to_numpy()[pool]) for c in MATCH_COLS])
is_drv = driver[pool]
drv_rows, ctl_rows = pool[is_drv], pool[~is_drv]
nn = NearestNeighbors(n_neighbors=a.k_match).fit(Z[~is_drv])
_, nbr = nn.kneighbors(Z[is_drv])
matched = ctl_rows[nbr]                                        # (n_driver_bins, k) row positions of matched non-drivers
balance = pd.DataFrame({c: [D[c].to_numpy()[drv_rows].mean(), D[c].to_numpy()[matched].mean(), D[c].to_numpy()[ctl_rows].mean()]
                        for c in MATCH_COLS}, index=["driver", "matched", "all non-driver"])
log("\n3. CONTRAST: matching balance (means)\n" + balance.round(3).to_string())
# bootstrap units: runs of consecutive driver bins on the same chromosome (one gene usually spans several bins)
chr_d, idx_d = D["chr"].to_numpy()[drv_rows], D["idx"].to_numpy()[drv_rows]
run_id = np.cumsum(np.r_[True, (chr_d[1:] != chr_d[:-1]) | (idx_d[1:] != idx_d[:-1] + 1)])
runs = [np.where(run_id == r)[0] for r in np.unique(run_id)]
role_d = role[drv_rows]
def contrast_stats(res_d, res_m):
    """res_d: residuals of driver bins (n); res_m: residuals of their matches (n, k).
    signed: driver - mean of its matches (a residual shift shared by all element bins cancels).
    abs / auc: deviations from the LOCAL matched baseline, so a shared shift also cancels here (raw |residual| does not:
    a shift of all element bins inflates |res| of drivers on one side and deflates it on the other - seen on toy data).
    Driver: |res_d - mean(matches)|; match j: |res_mj - mean(other k-1 matches)| (leave-one-out, so the match is not
    compared with itself). The two baselines have slightly different noise (1/k vs 1/(k-1) of the residual variance):
    with k = 5 this makes the driver deviations ~2% smaller in SD, i.e. a small conservative bias."""
    k = res_m.shape[1]
    base_d = res_m.mean(axis=1)
    signed = (res_d - base_d).mean()
    dev_d = np.abs(res_d - base_d)
    loo_mean = (res_m.sum(axis=1, keepdims=True) - res_m) / (k - 1)
    dev_m = np.abs(res_m - loo_mean)
    absolute = (dev_d - dev_m.mean(axis=1)).mean()
    lab = np.r_[np.ones(len(dev_d)), np.zeros(dev_m.size)]
    auc = roc_auc_score(lab, np.r_[dev_d, dev_m.ravel()])
    return np.array([signed, absolute, auc])
STAT_NAMES = ["signed_diff", "abs_dev_diff", "auc_abs_dev"]
GROUPS = {"all": np.ones(len(drv_rows), bool), "Act": role_d == "Act", "LoF": role_d == "LoF"}
boot_picks = [rng.integers(0, len(runs), len(runs)) for _ in range(a.n_boot)]   # same resamples for every arm (paired)
rows, BOOT = [], {}
for name, pred in PRED.items():
    res = y - pred
    for g, gmask in GROUPS.items():
        sel = np.where(gmask)[0]
        if len(sel) < 10:
            log(f"Group {g}: only {len(sel)} driver bins, skipped"); continue
        est = contrast_stats(res[drv_rows[sel]], res[matched[sel]])
        boots = []
        for pick in boot_picks:
            rows_b = np.concatenate([runs[r] for r in pick])
            rows_b = rows_b[gmask[rows_b]]
            boots.append(contrast_stats(res[drv_rows[rows_b]], res[matched[rows_b]]) if len(rows_b) >= 5 else [np.nan] * 3)
        boots = np.array(boots); BOOT[(name, g)] = boots
        for s, sname in enumerate(STAT_NAMES):
            rows.append(dict(arm=name, group=g, n_driver_bins=len(sel), stat=sname, estimate=est[s],
                             ci_low=np.nanpercentile(boots[:, s], 2.5), ci_high=np.nanpercentile(boots[:, s], 97.5)))
for fs_name in FEATURE_SETS:                                   # paired B - A difference, same resamples
    for g in GROUPS:
        if (f"A_{fs_name}", g) not in BOOT: continue
        diff_boot = BOOT[(f"B_{fs_name}", g)] - BOOT[(f"A_{fs_name}", g)]
        for s, sname in enumerate(STAT_NAMES):
            est_a = [r["estimate"] for r in rows if r["arm"] == f"A_{fs_name}" and r["group"] == g and r["stat"] == sname][0]
            est_b = [r["estimate"] for r in rows if r["arm"] == f"B_{fs_name}" and r["group"] == g and r["stat"] == sname][0]
            rows.append(dict(arm=f"B-A_{fs_name}", group=g, n_driver_bins=int(GROUPS[g].sum()), stat=sname,
                             estimate=est_b - est_a, ci_low=np.nanpercentile(diff_boot[:, s], 2.5),
                             ci_high=np.nanpercentile(diff_boot[:, s], 97.5)))
CON = pd.DataFrame(rows)
CON.to_csv(f"{a.out_dir}/contrast.tsv", sep="\t", index=False)
log(f"Driver bins {len(drv_rows)} in {len(runs)} runs; matched to {a.k_match} non-driver bins each "
    f"(pool {len(ctl_rows)}); bootstrap {a.n_boot} over runs")
log(CON.round(4).to_string(index=False))

OUT = pd.DataFrame({"binID": D["binID"], "neutral": neutral, "driver": driver, "role": role, "logFC": y})
for name, pred in PRED.items():
    OUT[f"pred_{name}"] = pred
OUT.to_csv(f"{a.out_dir}/predictions.tsv.gz", sep="\t", index=False)
with open(f"{a.out_dir}/log.txt", "w") as f:
    f.write("\n".join(LOG) + "\n")
log(f"[{time.time() - t0:.0f}s] Output in {a.out_dir}/")
