#!/usr/bin/env python
"""
BURNER (27/11/2026) - step 2 of the gene-structure map idea: REDUCED REPRESENTATION of the gene x (cell x class) count
matrix written by gene_structure_map.py, and per-gene counts REBUILT from it (Gab's choices: GLM-PCA, K by held-out
deviance, one joint matrix for syn / mis / trunc).

MODEL (Poisson GLM-PCA with offset, Townes et al. 2019 Genome Biol, written here with numpy/scipy)
  rows g = genes, columns j = (map cell, class): 31 exonic cells x {syn, mis, trunc} + splice x {trunc} = 94 columns
  y[g,j] ~ Poisson(mu[g,j]),  log mu[g,j] = log sites[g,j] + a_j + b_g + sum_k U[g,k] V[j,k]
    sites[g,j]  possible substitutions of that class in that cell of gene g (exposure, fixed offset)
    a_j         column intercept (rate of the cell x class over all genes)
    b_g         gene intercept (overall rate of the gene)  -> K = 0 is "gene rate x cell-class rate", no shared pattern
    V[:,k]      pattern k over the 94 columns (shared by all genes); U[g,k] loading of gene g on pattern k
  Entries with sites = 0 are not modelled (no exposure). Penalty lambda/2 (||U||^2 + ||V||^2) (glmpca default 1).
  NOTE (27/11/2026): gene intercepts are SHRUNK towards 0 with penalty --gene_penalty/2 * sum_g b_g^2 (Gab's choice,
  option a), i.e. b_g ~ N(0, 1/gene_penalty): the overall level sits in the column intercepts a_j, each gene deviates from
  it only as far as its own counts support. Reason: the first cluster run (job 9512042) failed - with a free b_g, genes
  whose events all fall in the held-out fold get b_g -> -inf (no finite estimate on 89-98% zero entries), and genes with
  NO training entries started at b_g = log(1e-8/1e-300) = +672 and never moved -> held-out deviance ~1e14 for every K and
  no convergence for K >= 2. Consequence: per-gene totals are no longer reproduced exactly (sparse genes are pulled
  towards the average rate); the log reports how far.
  Fitted by L-BFGS on the penalized Poisson likelihood. After fitting, U V^T is rotated by SVD so that patterns are
  orthonormal and ordered by strength (the fit itself is unchanged by the rotation).

CHOICE OF K: entries with sites > 0 split at random into --n_folds folds; for each K in --k_grid the model is fitted
  without one fold and the Poisson deviance of the held-out entries is summed over folds. Chosen K = lowest mean held-out
  deviance (the 1-SE choice is also printed, not used).

PLACEBO (the representation must not invent differences between classes): synonymous counts are split at random into
  two halves A and B (binomial thinning, p = 0.5, sites halved for each), the joint matrix becomes synA / synB / mis /
  trunc, and the model is refitted with the chosen K. Since A and B are the same class by construction, their patterns
  and their rebuilt per-gene counts should agree: reported as correlation of the pattern weights on synA vs synB columns
  (same cell) and the distribution of log(rebuilt synA / rebuilt synB) per gene compared with the observed one.

OUTPUTS (--out_dir)
  k_selection.tsv        K, fold, held-out deviance, held-out entries
  patterns.tsv           column (cell, cls), intercept a_j, pattern weights V (orthonormal, chosen K)
  gene_loadings.tsv      gene, b_g, loadings on each pattern
  rebuilt_counts.tsv     gene, cell, cls, sites, events, rebuilt (fitted mean of the final model)
  gene_rebuilt_totals.tsv  gene, per class: sites, observed events, rebuilt events
  placebo_syn_split.tsv  per gene: observed and rebuilt synA, synB
  log.txt

Usage (repository root on the cluster, pilot env):
  python gene_structure_representation.py --map analyses/v4/BRCA/20000bp/gene_structure_map/gene_map_long.tsv \
    --out_dir analyses/v4/BRCA/20000bp/gene_structure_representation
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy.optimize import minimize

p = argparse.ArgumentParser()
p.add_argument("--map", required=True, help="gene_map_long.tsv of gene_structure_map.py")
p.add_argument("--k_grid", default="0,1,2,3,4,5,6,8,10,12,15")
p.add_argument("--n_folds", type=int, default=5)
p.add_argument("--penalty", type=float, default=1.0, help="L2 penalty on U and V (glmpca default 1)")
p.add_argument("--gene_penalty", type=float, default=1.0,
               help="L2 shrinkage of the gene intercepts towards 0 (b_g ~ N(0, 1/gene_penalty)); 0 = free intercepts")
p.add_argument("--max_iter", type=int, default=3000)
p.add_argument("--seed", type=int, default=1)
p.add_argument("--out_dir", default="gene_structure_representation")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
rng = np.random.default_rng(a.seed)
K_GRID = [int(x) for x in a.k_grid.split(",")]

LOG = []
def log(m):
    print(m, file=sys.stderr, flush=True)
    LOG.append(str(m))


# ---------------- matrix ----------------
def to_matrix(long_df, col_order=None):
    """Long table (gene, col, sites, events) -> dense arrays Y, S (genes x columns); S = 0 where no sites."""
    genes = np.array(sorted(long_df["gene"].unique()))
    cols = np.array(col_order) if col_order is not None else np.array(sorted(long_df["col"].unique()))
    gi = pd.Series(np.arange(len(genes)), index=genes)
    ci = pd.Series(np.arange(len(cols)), index=cols)
    Y = np.zeros((len(genes), len(cols)))
    S = np.zeros((len(genes), len(cols)))
    r, c = gi[long_df["gene"]].to_numpy(), ci[long_df["col"]].to_numpy()
    Y[r, c] = long_df["events"].to_numpy(float)
    S[r, c] = long_df["sites"].to_numpy(float)
    return genes, cols, Y, S


# ---------------- Poisson GLM-PCA with offset ----------------
def fit_glmpca(Y, S, mask, K, penalty, max_iter, seed, gene_penalty):
    """Fit log mu = log S + a_j + b_g + U V^T on the entries in mask (bool, genes x columns). Returns dict of params."""
    G, J = Y.shape
    logS = np.where(S > 0, np.log(np.where(S > 0, S, 1.0)), 0.0)
    m = mask.astype(float)
    # NOTE (27/11/2026): starting values = column rates over all genes, gene intercepts at 0 (their prior mean). The old
    # start (alternating closed-form updates) gave b_g = +672 to genes without training entries (see header).
    a_ = np.log(np.maximum((Y * m).sum(0), 0.5) / np.maximum((S * m).sum(0), 1.0))
    b_ = np.zeros(G)
    r = np.random.default_rng(seed)
    U0 = r.normal(0, 0.01, (G, K))
    V0 = r.normal(0, 0.01, (J, K))

    def unpack(x):
        a1 = x[:J]
        b1 = x[J:J + G]
        U = x[J + G:J + G + G * K].reshape(G, K)
        V = x[J + G + G * K:].reshape(J, K)
        return a1, b1, U, V

    def objective(x):
        a1, b1, U, V = unpack(x)
        eta = logS + a1[None, :] + b1[:, None] + U @ V.T
        eta = np.minimum(eta, 30.0)                     # guard against overflow in early iterations
        mu = np.exp(eta)
        nll = (np.sum(m * (mu - Y * eta)) + 0.5 * penalty * (np.sum(U ** 2) + np.sum(V ** 2))
               + 0.5 * gene_penalty * np.sum(b1 ** 2))
        R = m * (mu - Y)                                # d nll / d eta
        grad = np.concatenate([R.sum(0), R.sum(1) + gene_penalty * b1,
                               (R @ V + penalty * U).ravel(), (R.T @ U + penalty * V).ravel()])
        return nll, grad

    x0 = np.concatenate([a_, b_, U0.ravel(), V0.ravel()])
    res = minimize(objective, x0, jac=True, method="L-BFGS-B", options={"maxiter": max_iter, "maxfun": 2 * max_iter})
    a1, b1, U, V = unpack(res.x)
    mu = np.exp(np.minimum(logS + a1[None, :] + b1[:, None] + U @ V.T, 30.0)) * (S > 0)
    return {"a": a1, "b": b1, "U": U, "V": V, "mu": mu, "converged": res.success, "n_iter": res.nit, "message": str(res.message)}


def poisson_deviance(y, mu):
    """Sum of unit Poisson deviances 2 [y log(y/mu) - (y - mu)] (y = 0 -> 2 mu)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(y > 0, y * np.log(y / mu), 0.0)
    return float(np.sum(2 * (term - (y - mu))))


def rotate(U, V):
    """Rotate U V^T by SVD: patterns orthonormal (columns of V), loadings carry the scale, ordered by strength."""
    if U.shape[1] == 0:
        return U, V
    P, d, Qt = np.linalg.svd(U @ V.T, full_matrices=False)
    k = U.shape[1]
    return P[:, :k] * d[:k], Qt[:k].T


# ---------------- data ----------------
L = pd.read_csv(a.map, sep="\t")
L["col"] = L["cell"] + "|" + L["cls"]
CELL_ORDER = list(dict.fromkeys(sorted(L["cell"].unique(), key=lambda c: (
    ["first", "internal", "pen55", "last", "splice"].index(c.split("_")[0]), c))))
COLS = [f"{c}|{k}" for c in CELL_ORDER for k in ("syn", "mis", "trunc") if f"{c}|{k}" in set(L["col"])]
genes, cols, Y, S = to_matrix(L, COLS)
obs = S > 0
log(f"matrix: {len(genes)} genes x {len(cols)} columns; entries with sites {int(obs.sum())}; events {int(Y.sum())}")
for k in ("syn", "mis", "trunc"):
    jj = [i for i, c in enumerate(cols) if c.endswith("|" + k)]
    log(f"  {k:5s}: columns {len(jj)}, events {int(Y[:, jj].sum())}, entries with events "
        f"{int((Y[:, jj] > 0).sum())} of {int(obs[:, jj].sum())}")

# ---------------- K by held-out deviance ----------------
idx = np.argwhere(obs)
fold = rng.integers(0, a.n_folds, len(idx))
rows = []
for K in K_GRID:
    for f in range(a.n_folds):
        held = np.zeros_like(obs)
        held[idx[fold == f, 0], idx[fold == f, 1]] = True
        fit = fit_glmpca(Y, S, obs & ~held, K, a.penalty, a.max_iter, a.seed + f, a.gene_penalty)
        dev = poisson_deviance(Y[held], fit["mu"][held])
        rows.append({"K": K, "fold": f, "heldout_deviance": dev, "heldout_entries": int(held.sum()),
                     "converged": fit["converged"], "n_iter": fit["n_iter"]})
    d = pd.DataFrame(rows)
    d = d[d["K"] == K]
    log(f"K = {K:2d}: held-out deviance {d['heldout_deviance'].sum():.1f} (fold sd {d['heldout_deviance'].std():.1f}); "
        f"converged in {int(d['converged'].sum())}/{a.n_folds} folds")
KS = pd.DataFrame(rows)
KS.to_csv(f"{a.out_dir}/k_selection.tsv", sep="\t", index=False)
summ = KS.groupby("K")["heldout_deviance"].agg(["mean", "std"])
summ["se"] = summ["std"] / np.sqrt(a.n_folds)
K_best = int(summ["mean"].idxmin())
K_1se = int(summ.index[summ["mean"] <= summ.loc[K_best, "mean"] + summ.loc[K_best, "se"]].min())
log(f"chosen K (lowest mean held-out deviance): {K_best}; 1-SE rule would give {K_1se} (not used)")
if K_best == max(K_GRID):
    log("WARNING: the chosen K is the largest of the grid: the held-out optimum may lie beyond it, extend --k_grid")
if not KS["converged"].all():
    log(f"WARNING: {int((~KS['converged']).sum())} fits did not converge within --max_iter; held-out deviances of those "
        f"fits are approximate")

# ---------------- final fit ----------------
fit = fit_glmpca(Y, S, obs, K_best, a.penalty, a.max_iter, a.seed, a.gene_penalty)
log(f"final fit K = {K_best}: converged {fit['converged']} ({fit['n_iter']} iterations, {fit['message']}); "
    f"in-sample deviance {poisson_deviance(Y[obs], fit['mu'][obs]):.1f}")
U, V = rotate(fit["U"], fit["V"])
pat = pd.DataFrame({"column": cols, "cell": [c.split("|")[0] for c in cols], "cls": [c.split("|")[1] for c in cols],
                    "intercept": fit["a"], **{f"pattern_{k + 1}": V[:, k] for k in range(K_best)}})
pat.to_csv(f"{a.out_dir}/patterns.tsv", sep="\t", index=False)
pd.DataFrame({"gene": genes, "gene_intercept": fit["b"], **{f"loading_{k + 1}": U[:, k] for k in range(K_best)}}) \
    .to_csv(f"{a.out_dir}/gene_loadings.tsv", sep="\t", index=False)
gi_, ci_ = np.nonzero(obs)
RB = pd.DataFrame({"gene": genes[gi_], "cell": pat["cell"].to_numpy()[ci_], "cls": pat["cls"].to_numpy()[ci_],
                   "sites": S[gi_, ci_].astype(np.int64), "events": Y[gi_, ci_].astype(np.int64),
                   "rebuilt": fit["mu"][gi_, ci_]})
RB.to_csv(f"{a.out_dir}/rebuilt_counts.tsv", sep="\t", index=False)
TOT = RB.groupby(["gene", "cls"])[["sites", "events", "rebuilt"]].sum().unstack("cls")
TOT.columns = [f"{v}_{c}" for v, c in TOT.columns]
TOT.to_csv(f"{a.out_dir}/gene_rebuilt_totals.tsv", sep="\t")
# NOTE (27/11/2026): with shrunk gene intercepts the gene totals are not reproduced exactly: report how far they move
g_obs = Y.sum(1)
g_reb = fit["mu"].sum(1)
log(f"gene intercepts: sd {np.std(fit['b']):.3f} (gene_penalty {a.gene_penalty}); gene totals rebuilt / observed over all "
    f"genes {g_reb.sum() / g_obs.sum():.4f}; genes with 0 events: {int((g_obs == 0).sum())}, their rebuilt total median "
    f"{np.median(g_reb[g_obs == 0]) if (g_obs == 0).any() else float('nan'):.3f}; genes with >= 10 events: median "
    f"rebuilt/observed {np.median(g_reb[g_obs >= 10] / g_obs[g_obs >= 10]) if (g_obs >= 10).any() else float('nan'):.3f}")
for k in ("syn", "mis", "trunc"):
    log(f"  {k:5s}: observed {int(TOT[f'events_{k}'].sum())}, rebuilt {TOT[f'rebuilt_{k}'].sum():.1f}; "
        f"per-gene correlation observed vs rebuilt {np.corrcoef(TOT[f'events_{k}'].fillna(0), TOT[f'rebuilt_{k}'].fillna(0))[0, 1]:.3f}")

# ---------------- placebo: synonymous split in two halves ----------------
syn = L[L["cls"] == "syn"]
yA = rng.binomial(syn["events"].to_numpy(np.int64), 0.5)
half = []
for name, y in (("synA", yA), ("synB", syn["events"].to_numpy(np.int64) - yA)):
    h = syn.copy()
    h["cls"] = name
    h["events"] = y
    h["sites"] = h["sites"] * 0.5                      # each half has half of the exposure
    half.append(h)
LP = pd.concat(half + [L[L["cls"] != "syn"]], ignore_index=True)
LP["col"] = LP["cell"] + "|" + LP["cls"]
COLS_P = [f"{c}|{k}" for c in CELL_ORDER for k in ("synA", "synB", "mis", "trunc") if f"{c}|{k}" in set(LP["col"])]
gp, cp, YP, SP = to_matrix(LP, COLS_P)
fp = fit_glmpca(YP, SP, SP > 0, K_best, a.penalty, a.max_iter, a.seed, a.gene_penalty)
log(f"placebo fit K = {K_best}: converged {fp['converged']} ({fp['n_iter']} iterations)")
UP, VP = rotate(fp["U"], fp["V"])
colsP = list(cp)
jA = [colsP.index(f"{c}|synA") for c in CELL_ORDER if f"{c}|synA" in colsP]
jB = [colsP.index(f"{c}|synB") for c in CELL_ORDER if f"{c}|synB" in colsP]
for k in range(K_best):
    log(f"  placebo pattern {k + 1}: correlation of weights synA vs synB over cells {np.corrcoef(VP[jA, k], VP[jB, k])[0, 1]:.3f}")
log(f"  placebo column intercepts synA vs synB: correlation {np.corrcoef(fp['a'][jA], fp['a'][jB])[0, 1]:.3f}, "
    f"mean difference {np.mean(fp['a'][jA] - fp['a'][jB]):+.4f}")
obsA, obsB = YP[:, jA].sum(1), YP[:, jB].sum(1)
rebA, rebB = fp["mu"][:, jA].sum(1), fp["mu"][:, jB].sum(1)
PL = pd.DataFrame({"gene": gp, "obs_synA": obsA, "obs_synB": obsB, "rebuilt_synA": rebA, "rebuilt_synB": rebB})
PL.to_csv(f"{a.out_dir}/placebo_syn_split.tsv", sep="\t", index=False)
keep = (obsA + obsB) > 0
lr_obs = np.log((obsA[keep] + 0.5) / (obsB[keep] + 0.5))
lr_reb = np.log(rebA[keep] / rebB[keep])
log(f"  placebo per gene (genes with >= 1 synonymous event, n = {int(keep.sum())}): log(A/B) observed median "
    f"{np.median(lr_obs):+.3f} IQR {np.subtract(*np.percentile(lr_obs, [75, 25])):.3f}; rebuilt median "
    f"{np.median(lr_reb):+.3f} IQR {np.subtract(*np.percentile(lr_reb, [75, 25])):.3f}")

open(f"{a.out_dir}/log.txt", "w").write("\n".join(LOG) + "\n")
print(f"Output in {a.out_dir}/", file=sys.stderr)
