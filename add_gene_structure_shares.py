#!/usr/bin/env python
"""
BURNER (27/11/2026) - copy of a feature matrix with 5 extra columns: share of the bin's coding sites in each exon category
of the gene-structure map (Gab: test whether map categories are selection predictors with the class interaction test).
Same pattern as analyses/v4/deleteriousness/add_am_feature.py. The original dataset is never modified.

  gs_share_<cat> = (sites_syn + sites_mis + sites_trunc of category cat in the bin) / (all coding sites of the bin)
  cat in first / internal / pen55 / last / splice (splice = essential splice sites, truncating only)
  from bin_category_sites.tsv of gene_structure_map.py (counted as the pipeline's sites_*). The 5 shares sum to 1.
Bins WITHOUT coding sites get the median share of the bins with coding sites: they have no exposure, so the class
interaction test drops them anyway (0 sites); without this the columns are NaN in > 50% of 20 kb bins and the test drops
the features before testing them (as happened with am_damaging_share).

Also writes a labels file = the input labels + one row per new feature (group 13_gene_structure). The prior labels in it
are ANNOTATION ONLY in pilot_class_interaction_test.py; they are set to "S" as for am_damaging_share.

Usage: python add_gene_structure_shares.py <dataset> <bin_category_sites.tsv> <labels.tsv> <out_dataset.tsv> <out_labels.tsv>
"""
import sys
import pandas as pd

dataset, cat_path, labels_path, out_path, out_labels = sys.argv[1:6]
CATS = ["first", "internal", "pen55", "last", "splice"]
D = pd.read_csv(dataset, sep="\t", dtype={"binID": str})
C = pd.read_csv(cat_path, sep="\t", dtype={"binID": str}).set_index("binID")
new_cols = [f"gs_share_{c}" for c in CATS]
clash = [c for c in new_cols if c in D.columns]
if clash:
    sys.exit(f"ERROR: {dataset} already has {clash}")
missing_bins = C.index.difference(D["binID"])
if len(missing_bins):
    sys.exit(f"ERROR: {len(missing_bins)} bins of {cat_path} are not in {dataset} (different backbone?): e.g. {list(missing_bins[:3])}")

cat_sites = pd.DataFrame({c: C[[k for k in C.columns if k.startswith(f"sites_{c}_")]].sum(axis=1) for c in CATS})
total = cat_sites.sum(axis=1)
shares = cat_sites.div(total, axis=0)
shares.columns = new_cols
D = D.merge(shares, left_on="binID", right_index=True, how="left")
# per-bin totals must match the dataset's coding sites (same counting rule)
tot_ds = D.set_index("binID")[["sites_syn", "sites_mis", "sites_trunc"]].sum(axis=1)
n_diff = int((tot_ds.reindex(total.index) != total).sum())
print(f"check: bins whose total coding sites differ from the dataset sites_*: {n_diff} of {len(total)}")
if n_diff:
    sys.exit("ERROR: category sites do not add up to the dataset sites_*: run gene_structure_map.py with the same backbone/GTF")
no_expo = D[new_cols[0]].isna()
for c in new_cols:
    D.loc[no_expo, c] = shares[c].median()
print(f"{out_path}: {len(D)} bins; with coding sites {int((~no_expo).sum())}; {int(no_expo.sum())} bins without coding sites "
      f"set to the median shares (not tested: no exposure)")
print("mean share over bins with coding sites: " + ", ".join(f"{c} {shares[c].mean():.3f}" for c in new_cols))
D.to_csv(out_path, sep="\t", index=False)

L = pd.read_csv(labels_path, sep="\t")
add = pd.DataFrame([{**{col: "S" for col in L.columns if col not in ("feature", "group")},
                     "feature": c, "group": "13_gene_structure"} for c in new_cols])[L.columns]
pd.concat([L, add], ignore_index=True).to_csv(out_labels, sep="\t", index=False)
print(f"{out_labels}: {len(L)} + {len(add)} features")
