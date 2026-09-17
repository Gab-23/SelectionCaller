#!/usr/bin/env bash
# Unit tests for the RepliSeq module's R scripts and helpers.
# Usage: bash .tests/unit/test_repliseq_scripts.sh
# Requires: Rscript (dplyr, tidyr, ggplot2, GenomicRanges, R.utils)

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

echo "=== 1) find_contiguous_subsets.R ==="
echo "    (bug: the classic seq(0) gotcha - seq(0) == c(1,0), NOT an empty sequence -"
echo "    combined with an if/else-if chain that could never add both the opening AND"
echo "    closing boundary when there was exactly one break. Confirmed narrow to that"
echo "    exact case (0, 2, 3, 4 breaks all worked in the original code too) - a real"
echo "    genome-wide blacklist typically has several regions per chromosome, so this"
echo "    edge case is unlikely to trigger in practice, but the fix is a zero-risk"
echo "    safety net verified as a non-regression on every other break count.)"
Rscript -e '
source("workflow/scripts/helpers/find_contiguous_subsets.R")
res0 <- find_contiguous_subsets(data.frame(pos = 1:10))
stopifnot(length(res0) == 1)
res1 <- find_contiguous_subsets(data.frame(pos = c(1,2,3,4,5,8,9,10)))
stopifnot(length(res1) == 2)
stopifnot(all(res1[[1]]$pos == 1:5))
stopifnot(all(res1[[2]]$pos == c(8,9,10)))
res2 <- find_contiguous_subsets(data.frame(pos = c(1,2,3,6,7,8,12,13)))
stopifnot(length(res2) == 3)
cat("OK\n")
'
check "no gap -> 1 contiguous subset; 1 gap -> 2; 2 gaps -> 3" "$?" "0"

echo "=== 2) process_repliseq.R (bugs, all confirmed with toy data covering the full"
echo "    ES/LS/G1 pipeline; the epigenomes-style header/chr-prefix bugs were caught by"
echo "    running the real Snakemake rule end-to-end, not just this isolated test) ==="

# NOTE (fix, 2026-09): process_repliseq.R now filters coverage files by an EXACT filename
# suffix built from assembly_code/window_mode/window_size/step/mapq_thr/cov_thr (see the
# script's own comment for why: files from a DIFFERENT parameter combination could otherwise
# slip in). This test previously never passed those five new CLI args and named its toy
# files just "x.filtered.bed" - which used to work when the script only checked for the
# literal substring "filtered", but silently stopped matching anything once the suffix check
# was added, failing with "expected exactly 1 coverage file ... found 0". Confirmed by
# reproducing the exact old (broken) invocation here before applying this fix - it failed
# with that exact message; the fixed invocation below (matching args + filenames) succeeds.
ASSEMBLY_CODE_TOY=hg00
WINDOW_MODE_TOY=binned
STEP_TOY=0
MAPQ_THR_TOY=20
COV_THR_TOY=5
COVERAGE_SUFFIX="${ASSEMBLY_CODE_TOY}_${WINDOW_MODE_TOY}_200_${STEP_TOY}_${MAPQ_THR_TOY}_${COV_THR_TOY}.filtered.bed"

mkdir -p "$TMP/coverage/SRR9001" "$TMP/coverage/SRR9002" "$TMP/coverage/SRR9003"
# coverage files WITHOUT a header (matching filter_coverage_repliseq.R's real col.names=F
# output), with a gap at pos 3 to exercise the find_contiguous_subsets fix too. Filenames now
# match the exact suffix process_repliseq.R expects (see NOTE above), not a placeholder name.
printf "chr1\t0\t200\t1_1\t50\nchr1\t200\t400\t1_2\t60\nchr1\t600\t800\t1_4\t55\nchr1\t800\t1000\t1_5\t45\n" > "$TMP/coverage/SRR9001/$COVERAGE_SUFFIX"
printf "chr1\t0\t200\t1_1\t80\nchr1\t200\t400\t1_2\t90\nchr1\t600\t800\t1_4\t85\nchr1\t800\t1000\t1_5\t70\n" > "$TMP/coverage/SRR9002/$COVERAGE_SUFFIX"
printf "chr1\t0\t200\t1_1\t100\nchr1\t200\t400\t1_2\t110\nchr1\t600\t800\t1_4\t95\nchr1\t800\t1000\t1_5\t90\n" > "$TMP/coverage/SRR9003/$COVERAGE_SUFFIX"
printf "chr1\t0\t200\t1_1\nchr1\t200\t400\t1_2\nchr1\t400\t600\t1_3\nchr1\t600\t800\t1_4\nchr1\t800\t1000\t1_5\n" > "$TMP/backbone_mask.bed"
printf "SRA_id\tsample_type\nSRR9001\tES\nSRR9002\tLS\nSRR9003\tG1\n" > "$TMP/samples.tsv"
printf "SRA_id\tlibrary_size\nSRR9001\t1000000\nSRR9002\t1000000\nSRR9003\t1000000\n" > "$TMP/library_size.tsv"

Rscript workflow/scripts/process_repliseq.R \
    --input_path_dir="$TMP/coverage/" --input_path_samples="$TMP/samples.tsv" \
    --input_path_library_size="$TMP/library_size.tsv" --input_path_backbone_mask="$TMP/backbone_mask.bed" \
    --output_path_data="$TMP/data.tsv" \
    --output_path_plot_distribution="$TMP/dist.pdf" --output_path_plot_landscape="$TMP/landscape.pdf" \
    --output_path_plot_raw_smooth="$TMP/raw_smooth.pdf" \
    --window_size=200 --smoothing_size=200 \
    --assembly_code=$ASSEMBLY_CODE_TOY --window_mode=$WINDOW_MODE_TOY --step=$STEP_TOY \
    --mapq_thr=$MAPQ_THR_TOY --cov_thr=$COV_THR_TOY

check "final data.tsv is NOT empty (bug: 'chr' prefix stripped from binID was never" \
    "$(( $(wc -l < "$TMP/data.tsv") - 1 ))" "4"
echo "     re-added before building the GRanges used for the final mask overlap, so"
echo "     replication_timing_gr's seqnames ('1'..'22') never matched mask_gr's"
echo "     ('chr1'..'chr22') and findOverlaps() always returned zero hits - the plots"
echo "     were and remain unaffected, since they derive chr/pos from binID independently)"
check "bin 1_1 present (header=F fix - it used to be silently eaten as a fake header)" \
    "$(grep -c '^1_1' "$TMP/data.tsv")" "1"

echo "=== 3) process_repliseq.R: multiple replicate files for one sample_type now fails"
echo "    explicitly instead of silently using only the first one ==="
mkdir -p "$TMP/coverage/SRR9004"
printf "chr1\t0\t200\t1_1\t120\nchr1\t200\t400\t1_2\t130\nchr1\t600\t800\t1_4\t115\nchr1\t800\t1000\t1_5\t100\n" > "$TMP/coverage/SRR9004/$COVERAGE_SUFFIX"
printf "SRA_id\tsample_type\nSRR9001\tES\nSRR9002\tLS\nSRR9003\tG1\nSRR9004\tG1\n" > "$TMP/samples_multi.tsv"
printf "SRA_id\tlibrary_size\nSRR9001\t1000000\nSRR9002\t1000000\nSRR9003\t1000000\nSRR9004\t1000000\n" > "$TMP/library_size_multi.tsv"
if Rscript workflow/scripts/process_repliseq.R \
    --input_path_dir="$TMP/coverage/" --input_path_samples="$TMP/samples_multi.tsv" \
    --input_path_library_size="$TMP/library_size_multi.tsv" --input_path_backbone_mask="$TMP/backbone_mask.bed" \
    --output_path_data="$TMP/data_multi.tsv" \
    --output_path_plot_distribution="$TMP/dist2.pdf" --output_path_plot_landscape="$TMP/landscape2.pdf" \
    --output_path_plot_raw_smooth="$TMP/raw_smooth2.pdf" \
    --window_size=200 --smoothing_size=200 \
    --assembly_code=$ASSEMBLY_CODE_TOY --window_mode=$WINDOW_MODE_TOY --step=$STEP_TOY \
    --mapq_thr=$MAPQ_THR_TOY --cov_thr=$COV_THR_TOY 2>"$TMP/err.log"; then
    echo "  [FAIL] should have failed with 2 G1 replicates, but it succeeded instead"
    FAIL=1
else
    check "explicit, actionable error on multiple replicate files" \
        "$(grep -c 'expected exactly 1 coverage file' "$TMP/err.log")" "1"
fi

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
