#!/usr/bin/env bash
# Unit tests for the fragile_sites module's R scripts.
# Usage: bash .tests/unit/test_fragile_sites_scripts.sh

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

echo "=== 1) process_dist_to_fragile_sites.R (bug: colnames<- was missing 'binID', silently"
echo "    shifting every column name by one and losing the binID column entirely) ==="
printf "chr1\t0\t200\t1_1\t0\nchr1\t200\t400\t1_2\t-1\n" > "$TMP/dist.bed"
Rscript workflow/scripts/process_dist_to_fragile_sites.R \
    --input_path_dist_fragile_sites="$TMP/dist.bed" --output_path_dist_fragile_sites="$TMP/dist_out"
check "output has the expected 2-column header" "$(head -1 "$TMP/dist_out")" "$(printf 'binID\tdist.to.closest.FGS')"
check "-1 converted to NA" "$(tail -1 "$TMP/dist_out" | cut -f2)" "NA"
check "regular value passed through unchanged" "$(sed -n 2p "$TMP/dist_out" | cut -f2)" "0"

echo "=== 2) process_num_overlaps_fragile_sites.R ==="
printf "chr1\t0\t200\t1_1\t2\nchr1\t200\t400\t1_2\t0\n" > "$TMP/overlap.bed"
Rscript workflow/scripts/process_num_overlaps_fragile_sites.R \
    --input_path_num_overlaps_fragile_sites="$TMP/overlap.bed" --output_path_num_overlaps_fragile_sites="$TMP/overlap_out"
check "overlap counts passed through correctly" "$(sed -n 2p "$TMP/overlap_out" | cut -f2)" "2"

echo "=== 3) sort_fragile_sites.smk (bug: fragile_sites_li_{assembly}.bed is a manually-"
echo "    curated source file with no producing rule, and was never run through bedtools"
echo "    sort - if it happens to be in natural chromosome order (chr1,chr2,...,chr10) while"
echo "    the backbone (via sort_backbone's plain 'bedtools sort', no -g file) is in"
echo "    lexicographic order (chr1,chr10,chr11,...,chr2,...), bedtools closest/intersect"
echo "    refuses to run at all. Confirmed directly on a real cluster run and reproduced"
echo "    here with the same error message.)"
printf "chr1\t0\t100000\t1_1\nchr1\t100000\t200000\t1_2\nchr10\t0\t100000\t10_1\nchr10\t100000\t200000\t10_2\nchr2\t0\t100000\t2_1\nchr2\t100000\t200000\t2_2\nchr9\t0\t100000\t9_1\n" > "$TMP/backbone_raw.bed"
bedtools sort -i "$TMP/backbone_raw.bed" > "$TMP/backbone.sorted.bed"
printf "chr1\t50000\t50100\nchr2\t50000\t50100\nchr9\t50000\t50100\nchr10\t50000\t50100\n" > "$TMP/fragile_natural_order.bed"

RC_UNFIXED=$(bedtools closest -a "$TMP/backbone.sorted.bed" -b "$TMP/fragile_natural_order.bed" -t first -d > /dev/null 2>&1; echo $?)
check "confirmed: naturally-ordered fragile_sites file breaks bedtools closest against the alphabetically-sorted backbone" "$RC_UNFIXED" "1"

bedtools sort -i "$TMP/fragile_natural_order.bed" > "$TMP/fragile.sorted.bed"
RC_FIXED=$(bedtools closest -a "$TMP/backbone.sorted.bed" -b "$TMP/fragile.sorted.bed" -t first -d > /dev/null 2>&1; echo $?)
check "fix: running fragile_sites through the same plain bedtools sort resolves it" "$RC_FIXED" "0"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
