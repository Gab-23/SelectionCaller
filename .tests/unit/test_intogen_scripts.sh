#!/usr/bin/env bash
# Unit tests for the IntOGen module's R scripts.
# Usage: bash .tests/unit/test_intogen_scripts.sh

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

echo "=== 1) process_IntOGen.R: Act/LoF/ambiguous -> OG/TSG/ambiguous mapping ==="
echo "    (historical bug: a locale-dependent alphabetical sort() could silently swap the"
echo "    TSG <-> ambiguous labels. Here we test all 3 roles together, in the order they"
echo "    appear in the input file, to make sure the mapping does not depend in any way"
echo "    on input ordering or system locale.)"

cat > "$TMP/IntOGen.tsv" << 'TSV'
SYMBOL	TRANSCRIPT	COHORT	CANCER_TYPE	ROLE	IS_DRIVER
GENEA	ENST001	COHORT1	BRCA	Act	True
GENEB	ENST002	COHORT1	BRCA	LoF	True
GENEC	ENST003	COHORT1	BRCA	ambiguous	True
GENED	ENST004	COHORT1	BRCA	LoF	False
TSV
printf "chr1\t99\t500\tENST001\tGENEA\tprotein_coding\ttranscript\nchr1\t599\t900\tENST002\tGENEB\tprotein_coding\ttranscript\nchr2\t99\t400\tENST003\tGENEC\tprotein_coding\ttranscript\n" > "$TMP/transcripts.bed"

Rscript workflow/scripts/process_IntOGen.R \
    --input_path_IntOGen="$TMP/IntOGen.tsv" \
    --input_path_protein_coding_transcript="$TMP/transcripts.bed" \
    --tumor_type_IntOGen=BRCA --cohort=COHORT1 \
    --output_path="$TMP/out.bed"

check "Act -> OG (GENEA)"              "$(grep GENEA "$TMP/out.bed" | cut -f6)" "OG"
check "LoF -> TSG (GENEB)"             "$(grep GENEB "$TMP/out.bed" | cut -f6)" "TSG"
check "ambiguous -> ambiguous (GENEC)" "$(grep GENEC "$TMP/out.bed" | cut -f6)" "ambiguous"
check "IS_DRIVER=False excluded (GENED absent)" "$(grep -c GENED "$TMP/out.bed" || true)" "0"

echo "=== 2) process_IntOGen.R: unexpected ROLE value -> must fail explicitly, not silently ==="
cat > "$TMP/IntOGen_bad.tsv" << 'TSV'
SYMBOL	TRANSCRIPT	COHORT	CANCER_TYPE	ROLE	IS_DRIVER
GENEA	ENST001	COHORT1	BRCA	SomethingWeird	True
TSV
if Rscript workflow/scripts/process_IntOGen.R \
    --input_path_IntOGen="$TMP/IntOGen_bad.tsv" \
    --input_path_protein_coding_transcript="$TMP/transcripts.bed" \
    --tumor_type_IntOGen=BRCA --cohort=COHORT1 \
    --output_path="$TMP/out_bad.bed" 2>"$TMP/err.log"; then
    echo "  [FAIL] should have failed on an unrecognized ROLE, but it succeeded instead"
    FAIL=1
else
    check "explicit error on unrecognized ROLE" "$(grep -c 'unexpected ROLE value' "$TMP/err.log")" "1"
fi

echo "=== 3) process_dist_to_cancer_genes.R ('.'/-1 -> NA conversion) ==="
printf "chr1\t0\t200\t1_1\t.\t.\t-1\n" > "$TMP/dist_tsg.bed"
printf "chr1\t0\t200\t1_1\t.\t.\t-1\n" > "$TMP/dist_og.bed"
Rscript workflow/scripts/process_dist_to_cancer_genes.R \
    --input_path_TSG="$TMP/dist_tsg.bed" --input_path_OG="$TMP/dist_og.bed" \
    --output_path_TSG="$TMP/out_tsg" --output_path_OG="$TMP/out_og"
check "TSG: '.'/-1 converted to NA" "$(tail -1 "$TMP/out_tsg" | cut -f2)" "NA"
check "OG: '.'/-1 converted to NA"  "$(tail -1 "$TMP/out_og" | cut -f2)" "NA"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
