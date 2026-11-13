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
# NOTE (fix, 2026-10): process_repliseq.R now requires smoothing_size to cover at least 4
# bins of window_size (loess, degree = 1, needs >= 4 points per local fit). The toy calls
# below used smoothing_size = window_size = 200 (1 bin) and were stopped by that guard, so
# they now use 800 (= 4 x 200). Expected values are unchanged: the toy's contiguous stretches
# have 2 bins each, below the script's > 5 bins threshold for smoothing, so the smoothing size
# does not affect the output. Section 4 checks the guard itself.
ASSEMBLY_CODE_TOY=hg00
WINDOW_MODE_TOY=binned
STEP_TOY=0
MAPQ_THR_TOY=20
COV_THR_TOY=5
WINDOW_SIZE_TOY=200
SMOOTHING_SIZE_TOY=800
COVERAGE_SUFFIX="${ASSEMBLY_CODE_TOY}_${WINDOW_MODE_TOY}_${WINDOW_SIZE_TOY}_${STEP_TOY}_${MAPQ_THR_TOY}_${COV_THR_TOY}.filtered.bed"

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

run_repliseq() {   # $1 samples, $2 library sizes, $3 output prefix, $4 smoothing_size
    Rscript workflow/scripts/process_repliseq.R \
        --input_path_dir="$TMP/coverage/" --input_path_samples="$1" \
        --input_path_library_size="$2" --input_path_backbone_mask="$TMP/backbone_mask.bed" \
        --output_path_data="$3_data.tsv" \
        --output_path_plot_distribution="$3_dist.pdf" --output_path_plot_landscape="$3_landscape.pdf" \
        --output_path_plot_raw_smooth="$3_raw_smooth.pdf" --output_path_plot_replicates="$3_replicates.pdf" \
        --window_size=$WINDOW_SIZE_TOY --smoothing_size="$4" \
        --assembly_code=$ASSEMBLY_CODE_TOY --window_mode=$WINDOW_MODE_TOY --step=$STEP_TOY \
        --mapq_thr=$MAPQ_THR_TOY --cov_thr=$COV_THR_TOY
}

RC=0
run_repliseq "$TMP/samples.tsv" "$TMP/library_size.tsv" "$TMP/main" "$SMOOTHING_SIZE_TOY" > "$TMP/main.log" 2>&1 || RC=$?
check "process_repliseq.R exits successfully (smoothing_size = 4 bins)" "$RC" "0"
[ "$RC" -eq 0 ] || cat "$TMP/main.log"
check "final data.tsv is NOT empty (bug: 'chr' prefix stripped from binID was never" \
    "$(( $(wc -l < "$TMP/main_data.tsv") - 1 ))" "4"
echo "     re-added before building the GRanges used for the final mask overlap, so"
echo "     replication_timing_gr's seqnames ('1'..'22') never matched mask_gr's"
echo "     ('chr1'..'chr22') and findOverlaps() always returned zero hits - the plots"
echo "     were and remain unaffected, since they derive chr/pos from binID independently)"
check "bin 1_1 present (header=F fix - it used to be silently eaten as a fake header)" \
    "$(grep -c '^1_1' "$TMP/main_data.tsv")" "1"

echo "=== 3) process_repliseq.R: replicates of a sample_type are SUMMED (coverage and library size) before CPM ==="
mkdir -p "$TMP/coverage/SRR9004" "$TMP/coverage/SRR9005" "$TMP/coverage/SRR9006" "$TMP/coverage/SRRPOOL"
printf "chr1\t0\t200\t1_1\t120\nchr1\t200\t400\t1_2\t130\nchr1\t600\t800\t1_4\t115\nchr1\t800\t1000\t1_5\t100\n" > "$TMP/coverage/SRR9004/$COVERAGE_SUFFIX"
printf "SRA_id\tsample_type\nSRR9001\tES\nSRR9002\tLS\nSRR9003\tG1\nSRR9004\tG1\n" > "$TMP/samples_multi.tsv"
printf "SRA_id\tlibrary_size\nSRR9001\t1000000\nSRR9002\t1000000\nSRR9003\t1000000\nSRR9004\t1000000\n" > "$TMP/library_size_multi.tsv"
# reference: one G1 sample whose coverage is SRR9003 + SRR9004 summed by hand, library size summed (2,000,000)
printf "chr1\t0\t200\t1_1\t220\nchr1\t200\t400\t1_2\t240\nchr1\t600\t800\t1_4\t210\nchr1\t800\t1000\t1_5\t190\n" > "$TMP/coverage/SRRPOOL/$COVERAGE_SUFFIX"
printf "SRA_id\tsample_type\nSRR9001\tES\nSRR9002\tLS\nSRRPOOL\tG1\n" > "$TMP/samples_pool.tsv"
printf "SRA_id\tlibrary_size\nSRR9001\t1000000\nSRR9002\t1000000\nSRRPOOL\t2000000\n" > "$TMP/library_size_pool.tsv"
RC=0; run_repliseq "$TMP/samples_multi.tsv" "$TMP/library_size_multi.tsv" "$TMP/multi" "$SMOOTHING_SIZE_TOY" > "$TMP/multi.log" 2>&1 || RC=$?
check "two G1 replicates: process_repliseq.R exits successfully" "$RC" "0"
[ "$RC" -eq 0 ] || cat "$TMP/multi.log"
RC=0; run_repliseq "$TMP/samples_pool.tsv" "$TMP/library_size_pool.tsv" "$TMP/pool" "$SMOOTHING_SIZE_TOY" > "$TMP/pool.log" 2>&1 || RC=$?
check "two G1 replicates: output identical to the hand-summed coverage (sum before CPM)" \
    "$(cmp -s "$TMP/multi_data.tsv" "$TMP/pool_data.tsv" && echo same || echo different)" "same"
# with the replicate column: replicate 1 = SRR9001/2/3, replicate 2 = SRR9005/6 (ES, LS) + SRR9004 (G1)
printf "chr1\t0\t200\t1_1\t55\nchr1\t200\t400\t1_2\t62\nchr1\t600\t800\t1_4\t50\nchr1\t800\t1000\t1_5\t48\n" > "$TMP/coverage/SRR9005/$COVERAGE_SUFFIX"
printf "chr1\t0\t200\t1_1\t78\nchr1\t200\t400\t1_2\t95\nchr1\t600\t800\t1_4\t80\nchr1\t800\t1000\t1_5\t72\n" > "$TMP/coverage/SRR9006/$COVERAGE_SUFFIX"
printf "SRA_id\tsample_type\tcell_line\treplicate\nSRR9001\tES\ttoy\t1\nSRR9002\tLS\ttoy\t1\nSRR9003\tG1\ttoy\t1\nSRR9005\tES\ttoy\t2\nSRR9006\tLS\ttoy\t2\nSRR9004\tG1\ttoy\t2\n" > "$TMP/samples_reps.tsv"
printf "SRA_id\tlibrary_size\nSRR9001\t1000000\nSRR9002\t1000000\nSRR9003\t1000000\nSRR9004\t1000000\nSRR9005\t1000000\nSRR9006\t1000000\n" > "$TMP/library_size_reps.tsv"
RC=0; run_repliseq "$TMP/samples_reps.tsv" "$TMP/library_size_reps.tsv" "$TMP/reps" "$SMOOTHING_SIZE_TOY" > "$TMP/reps.log" 2>&1 || RC=$?
check "replicate column (2 replicates): process_repliseq.R exits successfully" "$RC" "0"
[ "$RC" -eq 0 ] || cat "$TMP/reps.log"
check "replicate column: per-replicate PDF written" "$([ -s "$TMP/reps_replicates.pdf" ] && echo yes || echo no)" "yes"
check "replicate column: agreement between replicates reported in the log" "$(grep -c 'concordanza tra replicati' "$TMP/reps.log")" "1"
# a sample listed without a coverage file -> explicit error naming it
printf "SRA_id\tsample_type\nSRR9001\tES\nSRR9002\tLS\nSRR9003\tG1\nSRR9999\tG1\n" > "$TMP/samples_missing.tsv"
printf "SRA_id\tlibrary_size\nSRR9001\t1000000\nSRR9002\t1000000\nSRR9003\t1000000\nSRR9999\t1000000\n" > "$TMP/library_size_missing.tsv"
if run_repliseq "$TMP/samples_missing.tsv" "$TMP/library_size_missing.tsv" "$TMP/missing" "$SMOOTHING_SIZE_TOY" 2>"$TMP/err.log" >/dev/null; then
    echo "  [FAIL] a sample without coverage file was accepted"
    FAIL=1
else
    check "sample without coverage file -> explicit error naming the SRA_id" \
        "$(grep -c "expected exactly 1 coverage file for SRA_id 'SRR9999'" "$TMP/err.log")" "1"
fi
rm -rf "$TMP/coverage/SRR9004" "$TMP/coverage/SRR9005" "$TMP/coverage/SRR9006" "$TMP/coverage/SRRPOOL"

echo "=== 4) process_repliseq.R: smoothing_size must cover >= 4 bins of window_size ==="
echo "    (loess, degree = 1, needs >= 4 points per local fit: fewer either fails with"
echo "    'span is too small' or gives degenerate fits, so the script stops up front)"
for smooth in 200 600; do
    if run_repliseq "$TMP/samples.tsv" "$TMP/library_size.tsv" "$TMP/guard_$smooth" "$smooth" > "$TMP/guard_$smooth.log" 2>&1; then
        echo "  [FAIL] smoothing_size = $smooth bp (< 4 bins of 200 bp) was accepted"
        FAIL=1
    else
        check "smoothing_size = $smooth bp (< 4 bins) -> explicit error suggesting >= 800 bp" \
            "$(grep -c 'Use smoothing_size >= 800 bp' "$TMP/guard_$smooth.log")" "1"
    fi
done
check "no data file written when the guard stops the script" "$([ -e "$TMP/guard_200_data.tsv" ] && echo yes || echo no)" "no"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi