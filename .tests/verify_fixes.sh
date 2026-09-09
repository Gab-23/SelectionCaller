#!/usr/bin/env bash
# Verifies that all fixes from Phase 1 (backbone), Phase 2 (DepMap) and Phase 3 (IntOGen)
# are actually present in the local checkout. Does not run anything, just greps for the
# expected fixed lines. Run from the repository root:
#
#   bash .tests/verify_fixes.sh

set -uo pipefail

PASS=0
FAIL=0

check() {
    local desc="$1"; local file="$2"; local pattern="$3"; local expected_count="$4"

    if [ ! -f "$file" ]; then
        echo "  [MISSING FILE] $desc -- $file does not exist"
        FAIL=$((FAIL+1))
        return
    fi

    local got
    got=$(grep -cE "$pattern" "$file")

    if [ "$got" -ge "$expected_count" ]; then
        echo "  [OK] $desc"
        PASS=$((PASS+1))
    else
        echo "  [MISSING FIX] $desc -- pattern not found (or not often enough) in $file"
        FAIL=$((FAIL+1))
    fi
}

echo "=== Phase 1: backbone/annotations module ==="
check "get_actual_chromsizes_from_backbone.R: namespaced dplyr calls" \
    "workflow/scripts/get_actual_chromsizes_from_backbone.R" \
    "dplyr::group_by.*dplyr::summarise" 1

check "get_TSSs.R: coordinate fix applied to both strands" \
    "workflow/scripts/get_TSSs.R" \
    'strand == "-".*end - 1.*start - 1' 1

check "split_gtf.R: 1-based -> 0-based conversion present" \
    "workflow/scripts/split_gtf.R" \
    'sub_gtf\$start <- sub_gtf\$start - 1' 1

check "get_telomeres_proxy.R: chromosomes derived from input (not hard-coded chr1-22)" \
    "workflow/scripts/get_telomeres_proxy.R" \
    'chrom_sizes_canon\$V1' 1

check "config/sequencers_info: HiSeq_3000 typo fixed" \
    "config/sequencers_info" \
    "HiSeq_3000" 1

check "get_cytobands.smk: explicit grep pattern (not the unquoted bracket expression)" \
    "workflow/rules/centromeres_and_telomeres/get_cytobands.smk" \
    "grep -vE '_\|X\|M\|Y'" 1

check "get_cytobands.smk: output is .txt, not misleadingly .txt.gz" \
    "workflow/rules/centromeres_and_telomeres/get_cytobands.smk" \
    'cytoBand\.txt"' 1

check "get_centromeres_coords.smk: plain grep (matches the .txt, non-gzipped input)" \
    "workflow/rules/centromeres_and_telomeres/get_centromeres_coords.smk" \
    '^\s*grep acen' 1

check "generate_genes_on_bins_obj.smk: _v{gencode_version} present in input/output" \
    "workflow/rules/backbone/generate_genes_on_bins_obj.smk" \
    '_v\{gencode_version\}' 2

echo "=== Phase 2: DepMap module ==="
check "get_model_average_essentiality_score.R: drop = FALSE fix" \
    "workflow/scripts/get_model_average_essentiality_score.R" \
    "drop = FALSE" 2

check "get_essential_genes.smk: log path fixed (get_essential_genes, not filter_CRISPR_probabilities)" \
    "workflow/rules/DepMap/get_essential_genes.smk" \
    "logs/get_essential_genes/" 1

check "find_coords_essential_genes.smk: rule renamed to plural (matches file name)" \
    "workflow/rules/DepMap/find_coords_essential_genes.smk" \
    "^rule find_coords_essential_genes:" 1

echo "=== Phase 3: IntOGen module ==="
check "process_IntOGen.R: explicit name-based role mapping (not sort()-based)" \
    "workflow/scripts/process_IntOGen.R" \
    'role_map <- c\(Act = "OG", LoF = "TSG"' 1

check "process_IntOGen.R: fail-fast on unrecognized ROLE values" \
    "workflow/scripts/process_IntOGen.R" \
    "unexpected ROLE value" 1

echo
echo "-----------------------------------"
echo "Passed: $PASS   Failed: $FAIL"
if [ "$FAIL" -eq 0 ]; then
    echo "ALL FIXES ARE PRESENT"
    exit 0
else
    echo "SOME FIXES ARE MISSING -- see [MISSING FIX] / [MISSING FILE] lines above"
    exit 1
fi
