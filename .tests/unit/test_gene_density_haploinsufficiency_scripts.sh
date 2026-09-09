#!/usr/bin/env bash
# Unit tests for the gene_density and haploinsufficiency modules' R scripts.
# Usage: bash .tests/unit/test_gene_density_haploinsufficiency_scripts.sh

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

echo "=== 1) process_gene_density.R (bug: TSGs was reading the OGs file by mistake, and"
echo "    write.table used row.names=T, misaligning the header against the data) ==="
printf "chr1\t0\t200\t1_1\tGENEA,GENEB\tGENEA\n" > "$TMP/essential.bed"
printf "chr1\t0\t200\t1_1\tGENEA,GENEB\tGENEA_OG\n" > "$TMP/ogs.bed"
printf "chr1\t0\t200\t1_1\tGENEA,GENEB\tGENEB_TSG\n" > "$TMP/tsgs.bed"
Rscript workflow/scripts/process_gene_density.R \
    --input_path_essential="$TMP/essential.bed" --input_path_OGs="$TMP/ogs.bed" --input_path_TSGs="$TMP/tsgs.bed" \
    --output_path_essential="$TMP/out_essential" --output_path_OGs="$TMP/out_ogs" --output_path_TSGs="$TMP/out_tsgs"
check "OGs output has the OG gene, not the TSG one" "$(sed -n 2p "$TMP/out_ogs" | cut -f6)" "GENEA_OG"
check "TSGs output has the TSG gene, NOT the OG one (the actual bug)" \
    "$(sed -n 2p "$TMP/out_tsgs" | cut -f6)" "GENEB_TSG"
check "header and data have the same number of columns (row.names=F)" \
    "$(head -1 "$TMP/out_ogs" | awk -F'\t' '{print NF}')" "$(sed -n 2p "$TMP/out_ogs" | awk -F'\t' '{print NF}')"

echo "=== 2) process_haploinsufficiency.R ==="
echo "    (two bugs: (a) missing library(dplyr) - unqualified select/group_by/summarise"
echo "    crashed; (b) the annotated backbone bed has no header row, but was read with"
echo "    header=T, silently dropping whichever bin happened to be first in the file)"
printf "transcript;gene;pLI\nENST001;GENEA;0.95\nENST002;GENEA;0.90\nENST003;GENEB;0.10\n" > "$TMP/haplo.tsv"
printf "chr1\t0\t200\t1_1\tGENEA,GENEB\nchr1\t200\t400\t1_2\tGENEC\n" > "$TMP/backbone_no_header.tsv"
Rscript workflow/scripts/process_haploinsufficiency.R \
    --input_path_haplo="$TMP/haplo.tsv" --input_path_backbone="$TMP/backbone_no_header.tsv" \
    --output_path="$TMP/haplo_out" --haplo_thr=0.5 --C=1
check "does not crash (dplyr resolved correctly)" "$(wc -l < "$TMP/haplo_out" | tr -d ' ')" "3"
check "the FIRST bin in the file (1_1) is present, not silently dropped by a wrong header=T" \
    "$(grep -c '^1_1' "$TMP/haplo_out")" "1"
check "haploinsufficiency score correct for bin 1_1 (mean pLI of GENEA, above threshold)" \
    "$(grep '^1_1' "$TMP/haplo_out" | cut -f2)" "0.925"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
