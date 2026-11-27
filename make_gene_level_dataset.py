#!/usr/bin/env python
"""
BURNER (27/11/2026) - GENE-level input for pilot_class_interaction_test.py (Gab: run the per-feature supply vs selection
test with genes as units instead of bins; same test script, unchanged).

One row per gene (canonical transcript as in gene_structure_map.py / gene_level_drivers_classes.py):
  mut_n_syn / mis / trunc   events of the gene (gene_map_long.tsv of gene_structure_map.py, summed over cells;
                            positions inside the backbone bins)
  sites_syn / mis / trunc   possible substitutions of the gene (same file)
  features                  coding-sites-weighted mean of the BIN values over the bins holding the gene's coding sites,
                            weight = the gene's syn + mis + trunc sites in the bin (gene_bin_weighted_sites.tsv of
                            context_weighted_sites.py); same rule Gab approved for the weighted gene features (25/11).
                            Before averaging, the test's own NA rules are applied to the bin values (gene-level features
                            NA -> 0, distances NA -> 2 Mb); other NAs (ChromHMM, Hi-C, Repli-seq) are skipped in the mean
                            (a gene with all its bins NA stays NA).
  binID                     "<chrom>_<rank>": rank of the gene along its chromosome (order of its first bin, then name),
                            so that the test's moving-block bootstrap resamples blocks of CONSECUTIVE genes.
Context-weighted exposure for --exposure_file: gene_wsites.tsv (binID, wsites_syn / mis / trunc) = sum of the gene's
wsites over its bins (gene_bin_weighted_sites.tsv).
CHECK: per-gene sites in gene_bin_weighted_sites.tsv vs gene_map_long.tsv (both canonical transcript, inside the backbone;
the weighted file skips sites whose trinucleotide context contains N, so tiny differences are possible and are logged).

Usage: python make_gene_level_dataset.py <bin dataset> <gene_map_long.tsv> <gene_bin_weighted_sites.tsv> <out_dir>
Outputs in <out_dir>: gene_dataset.tsv, gene_wsites.tsv, gene_order.tsv (gene, chrom, first bin index, binID), log.txt
"""
import os, sys
import numpy as np, pandas as pd

ds_path, map_path, gbw_path, out_dir = sys.argv[1:5]
os.makedirs(out_dir, exist_ok=True)
LOG = []
def log(m):
    print(m, flush=True)
    LOG.append(str(m))

# ---------------- bin dataset: features only, with the test's NA rules (same lists as pilot_class_interaction_test.py) ----------------
D = pd.read_csv(ds_path, sep="\t", dtype={"binID": str})
ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith(("length_coding", "length_noncoding", "mut_n_"))]
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_", "expo_"))]
FEATS = [c for c in D.columns if c not in ALT + EXPOSURE + ["binID", "cds_bp"]]     # includes dist.to.closest.OG/TSG (controls)
GENE_LEVEL_NA_TO_ZERO = [c for c in FEATS
                         if c in ("ess_score", "haploinsufficiency_score", "complex_density", "n_complex_partners_trans",
                                  "n_PPIs_trans", "n_ohnolog_bin", "n_paralogs_bin")
                         or c.startswith(("GTEx_", "mean_n_ohnologs", "total_n_ohnologs", "mean_n_paralogs", "total_n_paralogs"))]
for c in GENE_LEVEL_NA_TO_ZERO:
    D[c] = D[c].fillna(0)
for c in [c for c in FEATS if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)
non_numeric = [c for c in FEATS if not pd.api.types.is_numeric_dtype(D[c])]
if non_numeric:
    sys.exit(f"ERROR: non-numeric feature columns in {ds_path}: {non_numeric}")
F = D.set_index("binID")[FEATS]
log(f"bin dataset: {len(D)} bins, {len(FEATS)} feature columns (incl. positive controls)")

# ---------------- gene x bin weights and context-weighted exposure ----------------
W = pd.read_csv(gbw_path, sep="\t", dtype={"binID": str})
w_gb = W.groupby(["gene", "binID"])["sites"].sum().rename("w").reset_index()        # syn + mis + trunc sites of the gene in the bin
missing_bins = set(w_gb["binID"]) - set(F.index)
if missing_bins:
    sys.exit(f"ERROR: {len(missing_bins)} bins of {gbw_path} are not in {ds_path}: e.g. {sorted(missing_bins)[:3]}")

# weighted mean per gene, skipping NA bins feature by feature
vals = F.loc[w_gb["binID"]].to_numpy(float)
w = w_gb["w"].to_numpy(float)[:, None]
has = ~np.isnan(vals)
num = pd.DataFrame(np.where(has, vals, 0.0) * w, columns=FEATS).groupby(w_gb["gene"].to_numpy()).sum()
den = pd.DataFrame(has * w, columns=FEATS).groupby(w_gb["gene"].to_numpy()).sum()
GF = num / den.where(den > 0)
log(f"genes with weights: {len(GF)}; features with NA genes: " +
    ", ".join(f"{c} {int(GF[c].isna().sum())}" for c in FEATS if GF[c].isna().any()))

WS = W.groupby(["gene", "cls"])["wsites"].sum().unstack("cls").reindex(columns=["syn", "mis", "trunc"]).fillna(0)
SW = W.groupby(["gene", "cls"])["sites"].sum().unstack("cls").reindex(columns=["syn", "mis", "trunc"]).fillna(0)

# ---------------- events and sites per gene from the map ----------------
M = pd.read_csv(map_path, sep="\t")
EV = M.groupby(["gene", "cls"])["events"].sum().unstack("cls").reindex(columns=["syn", "mis", "trunc"]).fillna(0)
ST = M.groupby(["gene", "cls"])["sites"].sum().unstack("cls").reindex(columns=["syn", "mis", "trunc"]).fillna(0)
common = ST.index.intersection(SW.index)
log(f"genes: map {len(ST)}, weighted-sites file {len(SW)}, common {len(common)}; only in map {len(ST.index.difference(SW.index))}, "
    f"only in weighted file {len(SW.index.difference(ST.index))}")
for c in ("syn", "mis", "trunc"):
    diff = (ST.loc[common, c] - SW.loc[common, c]).abs()
    log(f"CHECK sites {c:5s}: genes with different sites (map vs weighted file) {int((diff > 0).sum())}; total |difference| "
        f"{int(diff.sum())} of {int(ST.loc[common, c].sum())}")
genes = common.intersection(GF.index)

# ---------------- order along chromosomes -> binID for the block bootstrap ----------------
first_bin = w_gb.assign(chrom=w_gb["binID"].str.split("_").str[0], idx=w_gb["binID"].str.split("_").str[1].astype(int)) \
               .groupby("gene").agg(chrom=("chrom", "first"), idx=("idx", "min"))
multi_chrom = w_gb.assign(chrom=w_gb["binID"].str.split("_").str[0]).groupby("gene")["chrom"].nunique()
if (multi_chrom > 1).any():
    log(f"genes with bins on more than one chromosome (PAR genes): {int((multi_chrom > 1).sum())}: placed on their first chromosome")
O = first_bin.loc[genes].rename_axis("gene").reset_index().sort_values(["chrom", "idx", "gene"])
O["rank"] = O.groupby("chrom").cumcount()
O["binID"] = O["chrom"] + "_" + O["rank"].astype(str)
O.to_csv(f"{out_dir}/gene_order.tsv", sep="\t", index=False)

G = O[["gene", "binID"]].set_index("gene")
for c in ("syn", "mis", "trunc"):
    G[f"mut_n_{c}"] = EV.reindex(G.index)[c].fillna(0).astype(int)
    G[f"sites_{c}"] = ST.loc[G.index, c].astype(int)
G = G.join(GF.loc[G.index])
G.reset_index(drop=True).to_csv(f"{out_dir}/gene_dataset.tsv", sep="\t", index=False)
X = pd.DataFrame({"binID": G["binID"].to_numpy(), **{f"wsites_{c}": WS.reindex(G.index)[c].fillna(0).to_numpy() for c in ("syn", "mis", "trunc")}})
X.to_csv(f"{out_dir}/gene_wsites.tsv", sep="\t", index=False)
log(f"gene_dataset.tsv: {len(G)} genes; events syn {int(G['mut_n_syn'].sum())}, mis {int(G['mut_n_mis'].sum())}, "
    f"trunc {int(G['mut_n_trunc'].sum())}; genes per chromosome median {int(O.groupby('chrom').size().median())}")
open(f"{out_dir}/log.txt", "w").write("\n".join(LOG) + "\n")
