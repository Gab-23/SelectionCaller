#!/usr/bin/env bash
# Lightweight unit tests for the DepMap (essential genes) module's R scripts.
# Usage: bash .tests/unit/test_depmap_scripts.sh
# Requires: Rscript (base R, no extra package besides R.utils)

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

echo "=== 1) filter_CRISPR_models.R ==="
cat > "$TMP/Model.csv" << 'CSV'
ModelID,OncotreeCode,Sex
M1,BRCA,Female
M2,BRCA,Female
M3,BRCA,Male
M4,LUAD,Female
CSV
Rscript workflow/scripts/filter_CRISPR_models.R \
    --input_path="$TMP/Model.csv" --output_path="$TMP/models_out" \
    --gender_DepMap=Female --oncotree_code_DepMap=BRCA
check "correctly filters by oncotree+gender (only M1,M2)" "$(cat "$TMP/models_out")" "$(printf 'M1\nM2')"

echo "=== 2) get_essential_genes.R (DepMap 'GENE (id)' column-name parsing + gene missing from scores) ==="
cat > "$TMP/probabilities.csv" << 'CSV'
ModelID,"GENEA (1)","GENEB (2)","GENED (4)"
M1,0.9,0.1,0.6
M2,0.8,0.2,0.55
CSV
cat > "$TMP/scores.csv" << 'CSV'
ModelID,"GENEA (1)","GENEB (2)","GENEC (3)"
M1,-0.5,0.05,-0.2
M2,-0.6,0.02,-0.3
CSV
OUT=$(Rscript workflow/scripts/get_essential_genes.R \
    --input_path="$TMP/probabilities.csv" --scores_path="$TMP/scores.csv" \
    --output_path="$TMP/essential_out" \
    --essentiality_cutoff=0.5 --num_assays_cutoff=1 2>&1 >/dev/null)
check "GENEA essential, GENEB not, GENED dropped (missing from scores)" "$(cat "$TMP/essential_out")" "GENEA"
check "explicit warning about the dropped gene" "$(echo "$OUT" | grep -c 'dropping 1 gene')" "1"

echo "=== 3) get_model_average_essentiality_score.R (bug: a single essential gene) ==="
printf "chr1\t99\t500\tGENEA\n" > "$TMP/one_gene.bed"
Rscript workflow/scripts/get_model_average_essentiality_score.R \
    --input_path_scores="$TMP/scores.csv" --input_path_genes="$TMP/one_gene.bed" \
    --output_path="$TMP/avg_one"
check "correct mean with a single gene (no crash on colMeans)" "$(cut -f4,5 "$TMP/avg_one")" "$(printf 'GENEA\t-0.55')"

printf "chr1\t99\t500\tGENEA\nchr2\t599\t900\tGENEB\n" > "$TMP/two_genes.bed"
Rscript workflow/scripts/get_model_average_essentiality_score.R \
    --input_path_scores="$TMP/scores.csv" --input_path_genes="$TMP/two_genes.bed" \
    --output_path="$TMP/avg_two"
check "correct mean with two genes (non-regression)" "$(cut -f4,5 "$TMP/avg_two" | tr '\n' ';')" "GENEA	-0.55;GENEB	0.035;"

echo "=== 4) process_bin_level_essentiality.R (per-bin mean aggregation) ==="
printf "1_1\t-0.5\n1_1\t-0.6\n1_2\t-0.3\n" > "$TMP/bins.tsv"
Rscript workflow/scripts/process_bin_level_essentiality.R --input_path="$TMP/bins.tsv" --output_path="$TMP/bins_out"
check "aggregated per-bin mean (1_1 -> mean(-0.5,-0.6)=-0.55)" "$(grep '^1_1' "$TMP/bins_out" | cut -f2)" "-0.55"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
