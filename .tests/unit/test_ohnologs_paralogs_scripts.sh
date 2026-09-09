#!/usr/bin/env bash
# Unit tests for the ohnologs_and_paralogs module's R scripts.
# Usage: bash .tests/unit/test_ohnologs_paralogs_scripts.sh

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

printf "chr1\t0\t200\t1_1\tGENEA,GENEB\nchr1\t200\t400\t1_2\tGENEC\n" > "$TMP/backbone_annotated.bed"

echo "=== 1) compute_and_process_ohnologs.R (verified against the real OhnologList.RData"
echo "    schema you provided: 'n_ohnolog_mmpaper' already exists as a clean pre-computed"
echo "    integer column - no fix was needed here after all) ==="
Rscript -e '
ohnologs <- list(mmpaper = data.frame(
  Symbol = c("GENEA", "GENEC"),
  ohnologs_mmpaper = c("GENEX::GENEY", "GENEZ"),
  n_ohnolog_mmpaper = c(2L, 1L)
))
save(ohnologs, file = "'"$TMP"'/toy_ohnologs.RData")
'
Rscript workflow/scripts/compute_and_process_ohnologs.R \
    --input_path_ohnologs="$TMP/toy_ohnologs.RData" --input_path_backbone="$TMP/backbone_annotated.bed" \
    --output_path_ohnologs="$TMP/ohnologs_out"
check "mean_n_ohnologs.mmpaper matches the real pre-existing column (bin 1_1, GENEA)" \
    "$(sed -n 2p "$TMP/ohnologs_out" | cut -f3)" "2"
check "total_n_ohnologs.mmpaper matches (bin 1_2, GENEC)" \
    "$(sed -n 3p "$TMP/ohnologs_out" | cut -f4)" "1"

echo "=== 2) compute_and_process_paralogs.R (non-regression: this one already worked, since"
echo "    'n_paralogs' genuinely exists as a column in the raw paralogs data, unlike ohnologs) ==="
Rscript -e '
paralogs <- data.frame(
  symbol = c("GENEA", "GENEC"),
  n_paralogs = c(3, 1),
  paralogs = c("GENEP::GENEQ::GENER", "GENES")
)
save(paralogs, file = "'"$TMP"'/toy_paralogs.RData")
'
Rscript workflow/scripts/compute_and_process_paralogs.R \
    --input_path_paralogs="$TMP/toy_paralogs.RData" --input_path_backbone="$TMP/backbone_annotated.bed" \
    --output_path_paralogs="$TMP/paralogs_out"
check "mean_n_paralogs correct (bin 1_1, GENEA has 3 paralogs)" "$(sed -n 2p "$TMP/paralogs_out" | cut -f3)" "3"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
