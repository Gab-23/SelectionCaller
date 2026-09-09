#!/usr/bin/env bash
# Unit tests for the COSMIC module's R scripts.
# Usage: bash .tests/unit/test_cosmic_scripts.sh

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

echo "=== 1) filter_classification_COSMIC.R ==="
cat > "$TMP/classification.tsv" << 'TSV'
COSMIC_PHENOTYPE_ID	PRIMARY_SITE	PRIMARY_HISTOLOGY
PH001	breast	carcinoma
PH002	lung	carcinoma
TSV
Rscript workflow/scripts/filter_classification_COSMIC.R \
    --input_path="$TMP/classification.tsv" --output_path_phenotype_id="$TMP/class_out" \
    --COSMIC_parameters=breast_carcinoma
check "correctly filters PRIMARY_SITE+PRIMARY_HISTOLOGY (only PH001)" "$(cat "$TMP/class_out")" "PH001"

echo "    NOTE (documented, not fixed): COSMIC_parameters is split by '_' assuming a"
echo "    single-word site and a single-word histology (e.g. 'breast_carcinoma'). Real"
echo "    COSMIC site/histology names can themselves contain underscores (e.g."
echo "    'central_nervous_system'), which this position-based split cannot handle -"
echo "    this only works because the tumor type currently in production (BRCA ->"
echo "    'breast_carcinoma') has neither."
site_parts=$(Rscript -e 'x<-strsplit("central_nervous_system_glioma","_")[[1]]; cat(x[1])')
check "known limitation reproduced: multi-word site is NOT correctly extracted" "$site_parts" "central"

echo "=== 2) filter_samples_COSMIC.R + get_final_samples_COSMIC.R ==="
cat > "$TMP/samples.tsv" << 'TSV'
SAMPLE_NAME	GENDER	COSMIC_SAMPLE_ID	COSMIC_PHENOTYPE_ID
TCGA-S1	f	CS001	PH001
TCGA-S2	m	CS002	PH001
OTHER-S3	f	CS003	PH001
TSV
Rscript workflow/scripts/filter_samples_COSMIC.R \
    --input_path="$TMP/samples.tsv" --output_path="$TMP/samples_f" \
    --gender_COSMIC=f --sample_name_pattern=TCGA
check "filters by gender AND sample-name prefix (only TCGA-S1)" \
    "$(tail -n+2 "$TMP/samples_f" | cut -f1)" "TCGA-S1"

Rscript workflow/scripts/get_final_samples_COSMIC.R \
    --input_path_samples="$TMP/samples_f" --input_path_phenotype_ids="$TMP/class_out" \
    --output_path_samples="$TMP/final_samples"
check "joins on phenotype ID, outputs COSMIC_SAMPLE_ID (CS001)" "$(cat "$TMP/final_samples")" "CS001"

echo "=== 3) get_coding_mutations_description.R (join coordinates back to mutation type) ==="
printf "chr1\t100\t100\tCOSM1001\nchr1\t300\t302\tCOSM1003\n" > "$TMP/coords.sorted.bed"
cat > "$TMP/muts.filtered.bed" << 'TSV'
CHROMOSOME	GENOME_START	GENOME_STOP	GENOMIC_MUTATION_ID	MUTATION_DESCRIPTION
1	100	100	COSM1001	Substitution - Missense
1	300	302	COSM1003	Deletion - In frame
TSV
Rscript workflow/scripts/get_coding_mutations_description.R \
    --input_path_coords="$TMP/coords.sorted.bed" --input_path_muts="$TMP/muts.filtered.bed" \
    --output_path="$TMP/w_description.bed"
check "exactly 2 rows (no join duplication)" "$(wc -l < "$TMP/w_description.bed" | tr -d ' ')" "2"
check "COSM1001 correctly joined to its description" \
    "$(grep COSM1001 "$TMP/w_description.bed" | cut -f5)" "Substitution - Missense"

echo "=== 4) compute_coding_mutations_overlap.R (coordinate handling, see README for the"
echo "    full explanation of why GENOME_START/GENOME_STOP are used as-is, unconverted,"
echo "    and why that is correct given how this script builds GRanges directly) ==="
printf "chr1\t0\t200\t1_1\nchr1\t200\t400\t1_2\n" > "$TMP/backbone.bed"
printf "chr1\t100\t100\tCOSM1001\tSubstitution - Missense\nchr1\t300\t302\tCOSM1003\tDeletion - In frame\n" > "$TMP/muts_sorted.bed"
Rscript workflow/scripts/compute_coding_mutations_overlap.R \
    --input_path_mutation="$TMP/muts_sorted.bed" --input_path_backbone="$TMP/backbone.bed" \
    --output_path_mutation="$TMP/overlap_out.tsv"
check "single-base substitution (COSM1001) counts as 1bp in bin 1_1" \
    "$(grep '^1_1' "$TMP/overlap_out.tsv" | cut -f2)" "1"
check "3bp deletion (COSM1003) counts as 3bp in bin 1_2" \
    "$(grep '^1_2' "$TMP/overlap_out.tsv" | cut -f2)" "3"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
