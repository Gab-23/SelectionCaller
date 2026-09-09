#!/usr/bin/env bash
# Unit tests for the GDC_CNA module's R scripts.
# Usage: bash .tests/unit/test_gdc_scripts.sh
# Requires: Rscript (dplyr, ggplot2, gridExtra, R.utils)

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

echo "=== 1) get_samples_GDC.R (cohort+gender filter, drop patients with >1 sample/file) ==="
cat > "$TMP/clinical.tsv" << 'TSV'
project.project_id	cases.case_id	cases.submitter_id	demographic.gender
TCGA-BRCA	cid1	TCGA-P1	female
TCGA-BRCA	cid2	TCGA-P2	female
TCGA-BRCA	cid3	TCGA-P3	male
TCGA-BRCA	cid4	TCGA-P4	female
TCGA-LUAD	cid5	TCGA-P5	female
TSV
cat > "$TMP/sample_sheet.tsv" << 'TSV'
Case.ID	Sample.ID	File.ID	File.Name
TCGA-P1	S1	FILE001	file001.seg.txt
TCGA-P1	S1b	FILE001b	file001b.seg.txt
TCGA-P2	S2	FILE002	file002.seg.txt
TCGA-P3	S3	FILE003	file003.seg.txt
TCGA-P4	S4	FILE004	file004.seg.txt
TSV
Rscript workflow/scripts/get_samples_GDC.R \
    --input_path_clinical="$TMP/clinical.tsv" --input_path_sample_sheet="$TMP/sample_sheet.tsv" \
    --gender_GDC=female --cohort_GDC=TCGA-BRCA \
    --output_path_table="$TMP/table.tsv" --output_path_barcodes="$TMP/barcodes" --output_path_files="$TMP/files"
check "P1 dropped (2 files, ambiguous), P3 dropped (male), P5 dropped (LUAD)" \
    "$(cat "$TMP/barcodes" | tr '\n' ';')" "TCGA-P2;TCGA-P4;"
check "matching file IDs extracted" "$(cat "$TMP/files" | tr '\n' ';')" "FILE002;FILE004;"

echo "=== 2) filter_GDC_CNA_data.R (probe-count and |segment_mean| filter) ==="
printf "chr1\t100\t50100\t60\t0.85\nchr2\t100\t20100\t30\t-0.9\nchr1\t200\t50200\t55\t0.05\nchr2\t200\t20200\t5\t-0.02\n" > "$TMP/agg.bed"
Rscript workflow/scripts/filter_GDC_CNA_data.R \
    --input_path="$TMP/agg.bed" --output_path_filtered="$TMP/filtered.bed" \
    --output_path_plot_probes="$TMP/probes.pdf" --output_path_plot_segment="$TMP/segment.pdf" \
    --min_Num_Probes=10 --min_abs_Segment_Mean=0.1
check "2 of 4 segments pass both thresholds" "$(wc -l < "$TMP/filtered.bed" | tr -d ' ')" "2"
check "plot files created" "$([ -s "$TMP/probes.pdf" ] && [ -s "$TMP/segment.pdf" ] && echo yes)" "yes"

echo "=== 3) split logic (grep -v \"-\" / grep \"-\"): must not be confused by anything but segment_mean ==="
printf "chr1\t100\t50100\t60\t0.85\nchr2\t100\t20100\t30\t-0.9\n" > "$TMP/for_split.bed"
grep -v "-" "$TMP/for_split.bed" > "$TMP/ampl.bed"
grep "-" "$TMP/for_split.bed" > "$TMP/del.bed"
check "amplification (positive segment_mean) correctly isolated" "$(cut -f5 "$TMP/ampl.bed")" "0.85"
check "deletion (negative segment_mean) correctly isolated" "$(cut -f5 "$TMP/del.bed")" "-0.9"

echo "=== 4) process_GDC_CNA_overlap.R (frequency = count / num_patients, plots) ==="
printf "chr1\t0\t200\t1_1\t1\nchr1\t200\t400\t1_2\t0\n" > "$TMP/ampl_overlap.tsv"
printf "chr1\t0\t200\t1_1\t0\nchr1\t200\t400\t1_2\t1\n" > "$TMP/del_overlap.tsv"
Rscript workflow/scripts/process_GDC_CNA_overlap.R \
    --input_path_ampl="$TMP/ampl_overlap.tsv" --input_path_del="$TMP/del_overlap.tsv" \
    --output_path_ampl="$TMP/ampl_out" --output_path_del="$TMP/del_out" \
    --output_path_plot_ampl="$TMP/ampl_plot.pdf" --output_path_plot_del="$TMP/del_plot.pdf" \
    --num_patients=2
check "amplification frequency = 1/2 = 0.5 for bin 1_1" "$(grep '^1_1' "$TMP/ampl_out" | cut -f2)" "0.5"
check "deletion frequency = 1/2 = 0.5 for bin 1_2" "$(grep '^1_2' "$TMP/del_out" | cut -f2)" "0.5"
check "landscape plots created" "$([ -s "$TMP/ampl_plot.pdf" ] && [ -s "$TMP/del_plot.pdf" ] && echo yes)" "yes"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
