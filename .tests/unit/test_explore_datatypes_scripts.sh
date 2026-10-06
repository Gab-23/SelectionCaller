#!/usr/bin/env bash
# Unit tests for explore_datatypes.R (the .smk rule for this was entirely missing/empty
# in the original codebase - built from scratch, see README).
# Usage: bash .tests/unit/test_explore_datatypes_scripts.sh
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

echo "=== explore_datatypes.R (bug: two argument keys did not match their variable names -"
echo "    args\$input_path_ampl for a GC-content input, args\$output_path_mutations_barplot"
echo "    for a variable named output_path_num_mutations_barplot - functionally harmless as"
echo "    long as caller and script agree, but confusing; renamed for clarity since this rule"
echo "    is being wired up for the first time) ==="
echo "    The samples barplot (DepMap/GDC/COSMIC sample counts) moved to explore_cohort_atacseq.R,"
echo "    so this script no longer takes the sample inputs nor produces that plot."

printf "chr1\t0\t200\t1_1\t0.55\nchr1\t200\t400\t1_2\t0.42\n" > "$TMP/gc.bed"
printf "chr1\t99\t500\tENST001\tGENEA\tOG\n" > "$TMP/ogs.bed"
printf "chr2\t99\t500\tENST002\tGENEB\tTSG\n" > "$TMP/tsgs.bed"
printf "chr1\t99\t500\tGENEA\n" > "$TMP/essential.bed"
printf "chr1\t99\t500\tCOSM1\tSubstitution - Missense\nchr1\t600\t900\tCOSM2\tDeletion - In frame\n" > "$TMP/coding.bed"
printf "chr1\t250\t260\tCOSN1\n" > "$TMP/noncoding.bed"

RC=0
Rscript workflow/scripts/explore_datatypes.R \
    --input_path_GC="$TMP/gc.bed" \
    --input_path_OGs="$TMP/ogs.bed" \
    --input_path_TSGs="$TMP/tsgs.bed" \
    --input_path_essential="$TMP/essential.bed" \
    --input_path_coding_mutations="$TMP/coding.bed" \
    --input_path_noncoding_mutations="$TMP/noncoding.bed" \
    --output_path_GC_landscape="$TMP/out_GC_landscape.pdf" \
    --output_path_GC_distribution="$TMP/out_GC_dist.pdf" \
    --output_path_genes_barplot="$TMP/out_genes.pdf" \
    --output_path_num_mutations_barplot="$TMP/out_nmut.pdf" \
    --output_path_coding_mut_descriptions_barplot="$TMP/out_descr.pdf" \
    > "$TMP/stdout.log" 2>&1 || RC=$?
check "script exits successfully" "$RC" "0"
[ "$RC" -eq 0 ] || cat "$TMP/stdout.log"

for f in out_GC_landscape out_GC_dist out_genes out_nmut out_descr; do
    check "produces $f.pdf" "$([ -s "$TMP/$f.pdf" ] && echo yes)" "yes"
done
check "no samples barplot any more (moved to explore_cohort_atacseq.R)" "$(ls "$TMP" | grep -c samples || true)" "0"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi