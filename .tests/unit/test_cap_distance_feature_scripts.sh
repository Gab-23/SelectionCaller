#!/usr/bin/env bash
# Unit tests for cap_distance_feature.R (+ helpers/cap_distance_feature_helper.R): caps one processed
# distance table (binID + one distance column) at --dist_cap bp before combine_dataset.
# Usage: bash .tests/unit/test_cap_distance_feature_scripts.sh
# Requires: Rscript (ggplot2, R.utils)

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

run_cap() {   # $1 input, $2 cap, $3 output prefix
    Rscript workflow/scripts/cap_distance_feature.R \
        --input_path="$1" --dist_cap="$2" \
        --output_path="$3" --output_path_plot="$3.pdf" > "$3.log" 2>&1
}

echo "=== cap_distance_feature.R ==="

# chr1: one locus; chr2: no locus (all NA, as process_dist_to_* writes it); chrX after chr10
printf "binID\tdist.to.closest.OG\n1_1\t0\n1_2\t999999\n1_3\t1000000\n1_4\t5000000\n2_1\tNA\n2_2\tNA\n10_1\t20\nX_1\t3000000\n" > "$TMP/og.tsv"

RC=0; run_cap "$TMP/og.tsv" 1000000 "$TMP/capped" || RC=$?
check "script exits successfully" "$RC" "0"
[ "$RC" -eq 0 ] || cat "$TMP/capped.log"

check "header unchanged" "$(head -1 "$TMP/capped")" "$(head -1 "$TMP/og.tsv")"
check "binIDs and their order unchanged" "$(cut -f1 "$TMP/capped" | tr '\n' ' ')" "$(cut -f1 "$TMP/og.tsv" | tr '\n' ' ')"
check "distances above the cap are set to the cap, the others are unchanged" \
    "$(awk -F'\t' 'NR>1 && $2!="NA"{print $2}' "$TMP/capped" | tr '\n' ' ')" "0 999999 1000000 1000000 20 1000000 "
check "NA (chromosome with no locus) stays NA" "$(awk -F'\t' '$1 ~ /^2_/{print $2}' "$TMP/capped" | sort -u)" "NA"
check "no scientific notation in the output" "$(grep -c 'e+' "$TMP/capped" || true)" "0"
check "produces the diagnostic pdf" "$([ -s "$TMP/capped.pdf" ] && echo yes)" "yes"

RC=0; run_cap "$TMP/og.tsv" 1000000000 "$TMP/noop" || RC=$?
check "a cap above every distance leaves the table identical" "$(cmp -s "$TMP/og.tsv" "$TMP/noop" && echo identical)" "identical"

for bad in 0 -5 abc; do
    RC=0; run_cap "$TMP/og.tsv" "$bad" "$TMP/bad" || RC=$?
    check "rejects dist_cap=$bad" "$([ "$RC" -ne 0 ] && echo rejected)" "rejected"
done

printf "binID\tdist.to.closest.OG\tdist.to.closest.TSG\n1_1\t0\t5\n" > "$TMP/two_cols.tsv"
RC=0; run_cap "$TMP/two_cols.tsv" 1000000 "$TMP/two" || RC=$?
check "rejects a table with more than one distance column" "$([ "$RC" -ne 0 ] && grep -q 'exactly one distance column' "$TMP/two.log" && echo rejected)" "rejected"

check "chromosomes ordered numerically, then X/Y (diagnostic facets)" \
    "$(Rscript -e 'source("workflow/scripts/helpers/cap_distance_feature_helper.R"); cat(levels(split_binID(c("10_1","2_1","X_1","1_1"))$chr))')" "1 2 10 X"

check "contiguous runs split at gaps in the bin index and at chromosome changes (no line across gaps)" \
    "$(Rscript -e 'source("workflow/scripts/helpers/cap_distance_feature_helper.R"); cat(contiguous_run(c("1","1","1","1","2","2"), c(1,2,5,6,7,8)))')" "1 1 2 2 3 3"
check "contiguous runs do not depend on row order" \
    "$(Rscript -e 'source("workflow/scripts/helpers/cap_distance_feature_helper.R"); cat(contiguous_run(c("1","2","1","1"), c(5,1,1,2)))')" "2 3 1 1"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi