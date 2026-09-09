#!/usr/bin/env bash
# Unit tests for the epigenomes/Roadmap module's R scripts.
# Usage: bash .tests/unit/test_epigenomes_scripts.sh
# Requires: Rscript (dplyr, tidyr, magrittr, R.utils, GenomicRanges, rtracklayer)

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

echo "=== 1) process_epigenomes.R / compute_state_length.R (bin x chromatin-state overlap matrix) ==="
printf "chr1\t0\t200\t1_1\nchr1\t200\t400\t1_2\n" > "$TMP/backbone.bed"
printf "chr1\t0\t150\t1_TssA\nchr1\t150\t400\t2_TssFlnk\n" > "$TMP/E001.bed"
mkdir -p "$TMP/roadmap" "$TMP/checks/processing"
cp "$TMP/E001.bed" "$TMP/roadmap/E001_12_coreMarks_mnemonics.bed"
echo "$TMP/roadmap/E001_12_coreMarks_mnemonics.bed" > "$TMP/downloads_list"

Rscript workflow/scripts/process_epigenomes.R \
    --epigenomes_paths="$TMP/downloads_list" \
    --backbone_path="$TMP/backbone.bed" \
    --output_path="$TMP/checks/processing/breast_PrimaryCulture_noliftover_mnemonics_hg00_binned_200_0" \
    --is_test_call=F
tsv_path=$(cat "$TMP/checks/processing/breast_PrimaryCulture_noliftover_mnemonics_hg00_binned_200_0")
check "bin 1_1 correctly split across 2 overlapping states (150bp + 50bp)" \
    "$(grep '^1_1' "$tsv_path" | cut -f2,3)" "$(printf '150\t50')"
check "bin 1_2 fully covered by a single state (200bp), missing state filled with 0" \
    "$(grep '^1_2' "$tsv_path" | cut -f2,3)" "$(printf '0\t200')"

echo "=== 2) aggregate_epigenomes.R (bug: NA-fill must happen BEFORE averaging, not after) ==="
echo "    (a state column missing for one epigenome must be treated as 0bp for that epigenome"
echo "    when computing the mean across epigenomes -- filling NA with 0 only AFTER averaging"
echo "    would silently zero out the whole state whenever a single epigenome lacks it, even"
echo "    if other epigenomes have real signal there.)"

printf "binID\t1_TssA\t2_TssFlnk\t10_EnhA\n1_1\t100\t50\t50\n2_1\t0\t0\t200\n" > "$TMP/E001_out.tsv"
printf "binID\t1_TssA\t2_TssFlnk\n1_1\t80\t120\n2_1\t200\t0\n" > "$TMP/E002_out.tsv"
printf "%s\n%s\n" "$TMP/E001_out.tsv" "$TMP/E002_out.tsv" > "$TMP/agg_filelist.txt"

Rscript workflow/scripts/aggregate_epigenomes.R \
    --input_path="$TMP/agg_filelist.txt" --output_path="$TMP/aggregated.tsv"

check "1_TssA (present in both) -- mean(100,80)=90"    "$(grep '^1_1' "$TMP/aggregated.tsv" | cut -f2)" "90"
check "2_TssFlnk (present in both) -- mean(50,120)=85" "$(grep '^1_1' "$TMP/aggregated.tsv" | cut -f3)" "85"
check "10_EnhA (missing in E002, must count as 0) -- mean(50,0)=25, NOT 0 and NOT 50" \
    "$(grep '^1_1' "$TMP/aggregated.tsv" | cut -f4)" "25"
check "10_EnhA for bin 2_1 -- mean(200,0)=100" "$(grep '^2_1' "$TMP/aggregated.tsv" | cut -f4)" "100"

echo "=== 3) aggregate_epigenomes.R: numeric (not alphabetic) ordering of double-digit states ==="
check "column order is 1_TssA, 2_TssFlnk, 10_EnhA (numeric), not 10_EnhA, 1_TssA, 2_TssFlnk (alphabetic)" \
    "$(head -1 "$TMP/aggregated.tsv")" "$(printf 'binID\t1_TssA\t2_TssFlnk\t10_EnhA')"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
