#!/usr/bin/env python
"""
BURNER (26/11/2026) - mutation version of the neutral-training test (ATAC version: pilot_atac_neutral_training.py).
Gab's idea: train the background only on units without functional elements, test on the elements. With WES every
mutated bin holds genes, so here "functional" = bins of the IntOGen pan-cancer compendium genes (gene body overlap,
annotated backbone) and "neutral" = every other coding bin. Synonymous mutations INSIDE the held-out genes are the
neutral readout that ATAC lacks: they tell whether the background transfers into the elements.

Arms (same pilot_omics_background.py run, mutation classes only, same folds):
  A : standard pilot (compendium bins in training, predicted out-of-fold)
  B : pilot with --exclude_training_bins <compendium bins> (never trained on, predicted out-of-fold)
  then gene_level_drivers_classes.py on each arm's oof_predictions.tsv (unchanged).

Steps
  make_bins : writes the list of compendium bins (one binID per line) for --exclude_training_bins, and logs how many
              coding bins and mutations they hold.
  compare   : reads both arms and writes, in --out_dir:
    bin_level.tsv     per arm, class and bin set (compendium / other), from oof_predictions.tsv: sum observed /
                      sum predicted. Synonymous in compendium bins = TRANSFER readout (should be ~1 in B).
    gene_level.tsv    per arm, score and gene set (compendium / other), from gene_scores_<score>.tsv: sum O / sum E
                      (classes only) and share of genes with p < 0.05 and p < 0.01 (p = 10^-residual). Synonymous in
                      compendium genes = transfer readout; mis / trunc in OTHER genes = calibration (A and B should
                      agree: both predict those genes out-of-fold from models that never saw them).
    power.tsv         AP and AUROC of the cohort drivers (is_driver of the gene scores) per score and arm, universe =
                      all genes and = compendium genes only; B - A with a paired bootstrap over genes (drivers and
                      non-drivers resampled separately, same resamples for both arms).
Limits: the held-out set is defined by a driver list and AP is evaluated on drivers (cohort IntOGen drivers are mostly
inside the compendium), so B's gain measures the ABSORPTION effect in a favourable setting, not real-world performance.
Removing ~5% of coding bins keeps training almost unchanged, so transfer is easier than in the real version (all
genes removed, needs WGS): a positive result is necessary, not sufficient.
Usage (repository root, pilot env):
  python pilot_mut_neutral_training.py make_bins --dataset <DS> --out <dir>/compendium_bins.txt
  python pilot_mut_neutral_training.py compare --arm_a <dir A> --arm_b <dir B> --excluded_bins <dir>/compendium_bins.txt \
      --out_dir <dir>/comparison
  where <dir X> holds pilot/oof_predictions.tsv and drivers_gene_level/gene_scores_*.tsv
"""
import argparse, os, sys
import numpy as np, pandas as pd

p = argparse.ArgumentParser()
sub = p.add_subparsers(dest="step", required=True)
pb = sub.add_parser("make_bins")
pb.add_argument("--dataset", required=True)
pb.add_argument("--backbone_annotated", default="resources/backbone/hg38_binned_20000_0_v45.filtered.sorted.annotated.bed")
pb.add_argument("--compendium", default="resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv")
pb.add_argument("--min_exonic_bp", type=float, default=500, help="same as the pilot: training bins have >= this CDS bp")
pb.add_argument("--out", required=True)
pc = sub.add_parser("compare")
pc.add_argument("--arm_a", required=True); pc.add_argument("--arm_b", required=True)
pc.add_argument("--excluded_bins", required=True)
pc.add_argument("--compendium", default="resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv")
pc.add_argument("--scores", nargs="+", default=["syn", "mis", "trunc", "nonsyn", "nonsyn_fisher"])
pc.add_argument("--n_boot", type=int, default=1000); pc.add_argument("--seed", type=int, default=1)
pc.add_argument("--out_dir", required=True)
a = p.parse_args()

def need(path):
    if not os.path.exists(path):
        sys.exit(f"Missing input: {path}")
    return path

# =============================== make_bins ===============================
if a.step == "make_bins":
    comp_genes = set(pd.read_csv(need(a.compendium), sep="\t")["SYMBOL"])
    B = pd.read_csv(need(a.backbone_annotated), sep="\t", header=None, names=["chrom", "start", "end", "binID", "genes"], dtype=str)
    B["comp"] = B["genes"].str.split(",").map(lambda gs: [g for g in gs if g in comp_genes])
    sel = B[B["comp"].map(len) > 0]
    found = {g for gs in sel["comp"] for g in gs}
    D = pd.read_csv(need(a.dataset), sep="\t", dtype={"binID": str},
                    usecols=["binID", "cds_bp", "mut_n_syn", "mut_n_mis", "mut_n_trunc"]).set_index("binID")
    in_ds = sel["binID"][sel["binID"].isin(D.index)]
    coding = D.loc[in_ds][D.loc[in_ds, "cds_bp"] >= a.min_exonic_bp]
    train_all = D[D["cds_bp"] >= a.min_exonic_bp]
    print(f"Compendium: {len(comp_genes)} genes, {len(found)} found in the annotated backbone; missing (symbol not in GENCODE "
          f"protein-coding genes): {sorted(comp_genes - found)}")
    print(f"Compendium bins: {len(sel)} in the backbone, {len(in_ds)} in the dataset, {len(coding)} with >= {a.min_exonic_bp:g} "
          f"CDS bp = {100 * len(coding) / len(train_all):.1f}% of the {len(train_all)} mutation training bins")
    for c in ("syn", "mis", "trunc"):
        n_c, n_all = coding[f"mut_n_{c}"].sum(), train_all[f"mut_n_{c}"].sum()
        print(f"  {c}: {n_c:.0f} of {n_all:.0f} events of the training bins ({100 * n_c / n_all:.1f}%) are in compendium bins")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    in_ds.to_csv(a.out, index=False, header=False)
    print(f"Written {a.out} ({len(in_ds)} binIDs)")
    sys.exit(0)

# =============================== compare ===============================
os.makedirs(a.out_dir, exist_ok=True)
rng = np.random.default_rng(a.seed)
LOG = []
def log(msg):
    print(msg, flush=True); LOG.append(str(msg))
excluded = set(pd.read_csv(need(a.excluded_bins), header=None, dtype=str)[0])
comp_genes = set(pd.read_csv(need(a.compendium), sep="\t")["SYMBOL"])
ARMS = {"A": a.arm_a, "B": a.arm_b}

# ---- bin level: observed / predicted per class in compendium vs other bins ----
rows = []
for arm, d in ARMS.items():
    O = pd.read_csv(need(f"{d}/pilot/oof_predictions.tsv"), sep="\t", dtype={"binID": str})
    O = O[O["omic"].isin(["mut_syn", "mut_mis", "mut_trunc"])]
    for om, g in O.groupby("omic"):
        in_comp = g["binID"].isin(excluded)
        if arm == "B" and g.loc[in_comp, "in_training"].astype(str).eq("True").any():
            sys.exit(f"arm B {om}: some compendium bins have in_training = True: the exclusion did not work")
        for set_name, m in (("compendium", in_comp), ("other", ~in_comp)):
            rows.append(dict(arm=arm, cls=om.replace("mut_", ""), bins=set_name, n_bins=int(m.sum()),
                             n_in_training=int(g.loc[m, "in_training"].astype(str).eq("True").sum()),
                             observed=g.loc[m, "obs"].sum(), predicted=g.loc[m, "pred"].sum(),
                             obs_over_pred=g.loc[m, "obs"].sum() / g.loc[m, "pred"].sum()))
BL = pd.DataFrame(rows)
BL.to_csv(f"{a.out_dir}/bin_level.tsv", sep="\t", index=False)
log("BIN LEVEL (observed / predicted; synonymous in compendium bins = transfer readout)\n" + BL.round(3).to_string(index=False))

# ---- gene level: obs/exp and share of small p-values ----
def load_scores(d, score):
    G = pd.read_csv(need(f"{d}/drivers_gene_level/gene_scores_{score}.tsv"), sep="\t", index_col=0)
    G["is_driver"] = G["is_driver"].astype(str).eq("True")
    return G
rows = []
SC = {(arm, s): load_scores(d, s) for arm, d in ARMS.items() for s in a.scores}
for (arm, s), G in SC.items():
    in_comp = G.index.isin(comp_genes)
    for set_name, m in (("compendium", in_comp), ("other", ~in_comp)):
        Gm = G[m & G["residual"].notna()]
        pval = 10 ** (-Gm["residual"])
        row = dict(arm=arm, score=s, genes=set_name, n_genes=len(Gm), n_drivers=int(Gm["is_driver"].sum()),
                   frac_p05=(pval < 0.05).mean(), frac_p01=(pval < 0.01).mean())
        if "E" in Gm.columns:                                  # per-class scores; the Fisher score has no expectation
            row.update(O_sum=Gm["O"].sum(), E_sum=Gm["E"].sum(), O_over_E=Gm["O"].sum() / Gm["E"].sum())
        rows.append(row)
GL = pd.DataFrame(rows)
GL.to_csv(f"{a.out_dir}/gene_level.tsv", sep="\t", index=False)
log("\nGENE LEVEL (syn in compendium genes = transfer; mis/trunc in other genes = calibration, A and B should agree)\n"
    + GL.round(4).to_string(index=False))

# ---- power: AP / AUROC of cohort drivers, paired bootstrap B - A ----
def ap(s, l):
    """Average precision; tied scores form one step (as compare_with_dndscv.py); NaN scores at the bottom."""
    s = np.where(np.isfinite(s), s, -np.inf)
    o = np.argsort(-s, kind="mergesort"); ss, lo = s[o], l[o]
    last = np.r_[np.where(np.diff(ss) != 0)[0], len(ss) - 1]
    tp = np.cumsum(lo)[last]
    rec, prec = tp / max(l.sum(), 1), tp / (last + 1)
    return float(np.sum(np.diff(np.r_[0, rec]) * prec))
def auroc(s, l):
    s = np.where(np.isfinite(s), s, -np.inf)
    r = pd.Series(s).rank().to_numpy(); n1 = l.sum(); n0 = len(l) - n1
    return (r[l].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
rows = []
for s in a.scores:
    GA, GB = SC[("A", s)], SC[("B", s)]
    genes = GA.index.intersection(GB.index)
    if (GA.loc[genes, "is_driver"] != GB.loc[genes, "is_driver"]).any():
        sys.exit(f"{s}: driver labels differ between arms")
    for universe in ("all", "compendium"):
        u = genes if universe == "all" else genes[genes.isin(comp_genes)]
        lab = GA.loc[u, "is_driver"].to_numpy()
        if lab.sum() < 3:
            log(f"{s} {universe}: only {lab.sum()} drivers, skipped"); continue
        sa, sb = GA.loc[u, "residual"].to_numpy(float), GB.loc[u, "residual"].to_numpy(float)
        pos, neg = np.where(lab)[0], np.where(~lab)[0]
        est = {m: (f(sa, lab), f(sb, lab)) for m, f in (("AP", ap), ("AUROC", auroc))}
        boot = {m: [] for m in est}
        for _ in range(a.n_boot):
            idx = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
            for m, f in (("AP", ap), ("AUROC", auroc)):
                boot[m].append((f(sa[idx], lab[idx]), f(sb[idx], lab[idx])))
        for m in est:
            bt = np.array(boot[m]); diff = bt[:, 1] - bt[:, 0]
            rows.append(dict(score=s, universe=universe, n_genes=len(u), n_drivers=int(lab.sum()), metric=m,
                             A=est[m][0], A_lo=np.percentile(bt[:, 0], 2.5), A_hi=np.percentile(bt[:, 0], 97.5),
                             B=est[m][1], B_lo=np.percentile(bt[:, 1], 2.5), B_hi=np.percentile(bt[:, 1], 97.5),
                             B_minus_A=est[m][1] - est[m][0], diff_lo=np.percentile(diff, 2.5), diff_hi=np.percentile(diff, 97.5)))
PW = pd.DataFrame(rows)
PW.to_csv(f"{a.out_dir}/power.tsv", sep="\t", index=False)
log("\nPOWER (cohort drivers; B - A paired bootstrap over genes)\n" + PW.round(4).to_string(index=False))
with open(f"{a.out_dir}/log.txt", "w") as f:
    f.write("\n".join(LOG) + "\n")
log(f"Output in {a.out_dir}/")
