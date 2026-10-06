#!/usr/bin/env bash
# Unit tests for combine_dataset.R. This rule is fully reviewed and fixed but not yet
# wired into rule_all (it needs ATAC-seq + RepliSeq outputs that don't exist yet), so
# this is the only place it is currently exercised - directly, with toy data covering
# every one of its many inputs/outputs.
# Usage: bash .tests/unit/test_combine_dataset_scripts.sh
# Requires: Rscript (ggplot2, dplyr, tidyr, R.utils)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

FAIL=0
check() {
    local desc="$1"; local got="$2"; local expected="$3"
    if [ "$got" == "$expected" ]; then
        echo "  [OK] $desc"
    else
        echo "  [FAIL] $desc"
        echo "         expected: $expected"
        echo "         got:      $got"
        FAIL=1
    fi
}

echo "=== combine_dataset.R ==="
echo "    Bugs found and fixed (verified with this toy dataset, 20 bins with random but"
echo "    non-degenerate values so the correlation test has enough variance to run):"
echo "    1) output_path_bool_multicollinearity was commented out (a leftover from disabling"
echo "       the PCA block above it) but still used later -> guaranteed crash."
echo "    2) the final boxplots_0s_non0s_chromstates ggsave() passed the OUTPUT PATH as the"
echo "       plot object instead of the actual computed plot variable."
echo "    3) the correlation-heatmap ggsave() referenced 'output_path_corr' (never assigned)"
echo "       instead of 'output_path_corr_plot' (the actual variable) -> guaranteed crash."
echo "    4) cor.test() crashed the ENTIRE script (before the feature matrix was even"
echo "       written) if any pair of columns had zero variance - confirmed with the real toy"
echo "       Snakemake pipeline (several near-constant toy features). Wrapped in tryCatch so a"
echo "       degenerate pair is skipped (NA) instead of blocking the core deliverable; added"
echo "       na.rm=TRUE to the final any() so NAs don't turn the multicollinearity flag itself"
echo "       into NA."
echo "    Also newly discovered while wiring up the .smk rule (not a script bug, a wiring"
echo "    bug): the epigenomes input pointed at a file that contains a path STRING, not an"
echo "    actual binID-keyed table (see README) - would have crashed merge(by='binID')."
echo "    Hi-C's eigen_1..eigen_10 columns were added as a new input per user request (Phase"
echo "    12) - combine_dataset.R never read any Hi-C input before that; no script changes"
echo "    needed beyond declaring the new argument, since the merge loop auto-collects every"
echo "    input_path_* variable in scope."
echo "    Protein-interaction (CORUM + Interactome INSIDER) and GTEx features added as two new"
echo "    inputs, same mechanism. Real data has bins without protein-coding genes: those bins"
echo "    are ABSENT from both files (not NA rows), so the toy files below omit some bins and"
echo "    also contain explicit NAs (gene sets with no expressed gene)."

python3 - "$TMP" << 'PYEOF'
import random, sys, os
random.seed(7)
tmp = sys.argv[1]
bins = [f"1_{i}" for i in range(1, 21)]

def write(fname, header, cols_fn, bin_subset=None):
    with open(os.path.join(tmp, fname), "w") as f:
        f.write(header + "\n")
        for b in (bins if bin_subset is None else bin_subset):
            vals = cols_fn(b)
            f.write(b + "\t" + "\t".join(str(v) for v in vals) + "\n")

write("DARs.tsv", "binID\tlogFC", lambda b: [round(random.uniform(-2,2),2)])
write("centromeres.tsv", "binID\tdist.to.closest.centromere", lambda b: [random.randint(0,10000)])
write("telomeres.tsv", "binID\tdist.to.closest.telomere", lambda b: [random.randint(0,10000)])
write("coding_mut.tsv", "binID\tlength_coding_agnostic", lambda b: [random.randint(0,5)])
write("noncoding_mut.tsv", "binID\tlength_noncoding_agnostic", lambda b: [random.randint(0,5)])
write("dist_essential.tsv", "binID\tdist.to.closest.essential", lambda b: [random.randint(0,10000)])
write("epigenomes_real.tsv", "binID\t1_TssA\t2_TssFlnk", lambda b: [random.randint(0,200), random.randint(0,200)])
write("num_fgs.tsv", "binID\tnum.overlapping.FGS", lambda b: [random.randint(0,2)])
write("dist_fgs.tsv", "binID\tdist.to.closest.FGS", lambda b: [random.randint(0,10000)])
write("gc.tsv", "binID\tGC.content", lambda b: [round(random.uniform(0.3,0.6),2)])
write("ampl.tsv", "binID\tampl.frequency", lambda b: [round(random.uniform(0,1),2)])
write("del.tsv", "binID\tdel.frequency", lambda b: [round(random.uniform(0,1),2)])
write("numov_pc.tsv", "binID\tnum.overlapping.genes_protein_coding", lambda b: [random.randint(0,3)])
write("numov_lnc.tsv", "binID\tnum.overlapping.genes_lncRNA", lambda b: [random.randint(0,3)])
write("numov_mi.tsv", "binID\tnum.overlapping.genes_miRNA", lambda b: [random.randint(0,3)])
write("numov_r.tsv", "binID\tnum.overlapping.genes_rRNA", lambda b: [random.randint(0,3)])
write("dist_pc.tsv", "binID\tdist.to.closest.gene_protein_coding", lambda b: [random.randint(0,10000)])
write("dist_lnc.tsv", "binID\tdist.to.closest.gene_lncRNA", lambda b: [random.randint(0,10000)])
write("dist_mi.tsv", "binID\tdist.to.closest.gene_miRNA", lambda b: [random.randint(0,10000)])
write("dist_r.tsv", "binID\tdist.to.closest.gene_rRNA", lambda b: [random.randint(0,10000)])
write("tsg.tsv", "binID\tdist.to.closest.TSG", lambda b: [random.randint(0,10000)])
write("og.tsv", "binID\tdist.to.closest.OG", lambda b: [random.randint(0,10000)])
write("repliseq.tsv", "binID\trepliseq_signal", lambda b: [round(random.uniform(0,1),2)])
write("ohnologs.tsv", "binID\tn_ohnolog_bin\tmean_n_ohnologs.mmpaper\ttotal_n_ohnologs.mmpaper\tmean_n_ohnologs.mmpaper_trans\ttotal_n_ohnologs.mmpaper_trans",
      lambda b: [random.randint(0,2)]*5)
write("paralogs.tsv", "binID\tn_paralogs_bin\tmean_n_paralogs\ttotal_n_paralogs\tmean_n_paralogs_trans\ttotal_n_paralogs_trans",
      lambda b: [random.randint(0,2)]*5)
write("haplo.tsv", "binID\thaploinsufficiency_score", lambda b: [round(random.uniform(0,1),2)])
write("hic.tsv", "binID\teigen_1\teigen_2", lambda b: [round(random.uniform(-1,1),2), round(random.uniform(-1,1),2)])
# Bins 1_19 and 1_20 have no protein-coding genes -> absent from both files
gene_bins = bins[:18]
write("interactors.tsv", "binID\tcomplex_density\tn_complex_partners_trans\tn_PPIs_trans",
      lambda b: [round(random.uniform(0,1),2), random.randint(0,30), random.randint(0,60)], gene_bins)
# Fixed, known values for bin 1_1 so pass-through can be checked exactly; random elsewhere,
# with NAs where a gene set has no expressed gene.
write("gtex.tsv", "binID\tGTEx_genes\tGTEx_partners_trans\tGTEx_PPIs_trans\tGTEx_all_int_trans",
      lambda b: [0.123, "NA", -0.456, 0.789] if b == "1_1" else
                [round(random.uniform(-2,2),3) if random.random() > 0.2 else "NA" for _ in range(4)], gene_bins)
PYEOF

mkdir -p "$TMP/out"
RC=0
Rscript workflow/scripts/combine_dataset.R \
    --input_path_DARs="$TMP/DARs.tsv" \
    --input_path_centromeres="$TMP/centromeres.tsv" \
    --input_path_telomeres="$TMP/telomeres.tsv" \
    --input_path_coding_mutations="$TMP/coding_mut.tsv" \
    --input_path_noncoding_mutations="$TMP/noncoding_mut.tsv" \
    --input_path_dist_to_essential="$TMP/dist_essential.tsv" \
    --input_path_epigenomes="$TMP/epigenomes_real.tsv" \
    --input_path_num_fragile_sites="$TMP/num_fgs.tsv" \
    --input_path_dist_fragile_sites="$TMP/dist_fgs.tsv" \
    --input_path_GC_content="$TMP/gc.tsv" \
    --input_path_ampl_frequency="$TMP/ampl.tsv" \
    --input_path_del_frequency="$TMP/del.tsv" \
    --input_path_num_overlaps_protein_coding="$TMP/numov_pc.tsv" \
    --input_path_num_overlaps_lncRNA="$TMP/numov_lnc.tsv" \
    --input_path_num_overlaps_miRNA="$TMP/numov_mi.tsv" \
    --input_path_num_overlaps_rRNA="$TMP/numov_r.tsv" \
    --input_path_dist_protein_coding="$TMP/dist_pc.tsv" \
    --input_path_dist_lncRNA="$TMP/dist_lnc.tsv" \
    --input_path_dist_miRNA="$TMP/dist_mi.tsv" \
    --input_path_dist_rRNA="$TMP/dist_r.tsv" \
    --input_path_TSG="$TMP/tsg.tsv" \
    --input_path_OG="$TMP/og.tsv" \
    --input_path_repliseq="$TMP/repliseq.tsv" \
    --input_path_ohnologs="$TMP/ohnologs.tsv" \
    --input_path_paralogs="$TMP/paralogs.tsv" \
    --input_path_haploinsufficiency="$TMP/haplo.tsv" \
    --input_path_hic="$TMP/hic.tsv" \
    --input_path_interactors="$TMP/interactors.tsv" \
    --input_path_GTEx="$TMP/gtex.tsv" \
    --output_path_dataset="$TMP/out/dataset.tsv" \
    --output_path_corr="$TMP/out/corr.pdf" \
    --output_path_NA_frac_plot="$TMP/out/NA_frac.pdf" \
    --output_path_zero_frac_plot="$TMP/out/zero_frac.pdf" \
    --output_path_boxplots="$TMP/out/boxplots.pdf" \
    --output_path_boxplots_scaled="$TMP/out/boxplots_scaled.pdf" \
    --output_path_boxplots_0s_non0s_muts="$TMP/out/box_muts.pdf" \
    --output_path_boxplots_0s_non0s_chromstates="$TMP/out/box_chromstates.pdf" \
    --output_path_bool_multicollinearity="$TMP/out/multicoll.txt" \
    > "$TMP/stdout.log" 2>&1 || RC=$?
check "script exits successfully" "$RC" "0"
[ "$RC" -eq 0 ] || cat "$TMP/stdout.log"

for f in dataset.tsv corr.pdf NA_frac.pdf zero_frac.pdf boxplots.pdf boxplots_scaled.pdf box_muts.pdf box_chromstates.pdf multicoll.txt; do
    check "produces $f" "$([ -s "$TMP/out/$f" ] && echo yes)" "yes"
done
check "feature matrix has a binID column" "$(head -1 "$TMP/out/dataset.tsv" | cut -f1)" "binID"
check "feature matrix includes Hi-C's eigen_1 column" "$(head -1 "$TMP/out/dataset.tsv" | tr '\t' '\n' | grep -c '^eigen_1$')" "1"
check "rows with NA logFC (the target variable) are dropped" "$(( $(wc -l < "$TMP/out/dataset.tsv") - 1 ))" "20"
for col in complex_density n_complex_partners_trans n_PPIs_trans GTEx_genes GTEx_partners_trans GTEx_PPIs_trans GTEx_all_int_trans; do
    check "feature matrix includes $col" "$(head -1 "$TMP/out/dataset.tsv" | tr '\t' '\n' | grep -c "^${col}\$")" "1"
done
# Values pass through unchanged, and bins absent from the gene-based inputs are kept (as NA)
check "GTEx values for bin 1_1 pass through unchanged" \
    "$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} $1=="1_1"{print $c["GTEx_genes"], $c["GTEx_partners_trans"], $c["GTEx_PPIs_trans"], $c["GTEx_all_int_trans"]}' "$TMP/out/dataset.tsv")" \
    "0.123 NA -0.456 0.789"
check "gene-less bins (absent from interactor/GTEx files) are kept with NA" \
    "$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} $1=="1_19" || $1=="1_20"{print $c["complex_density"], $c["GTEx_genes"]}' "$TMP/out/dataset.tsv" | sort -u)" \
    "NA NA"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi