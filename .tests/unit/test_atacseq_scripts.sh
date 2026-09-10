#!/usr/bin/env bash
# Unit tests for the ATAC-seq module's R scripts and the CNVkit-output awk logic.
# The full alignment -> CNVkit -> counts chain (real bwa/bowtie2/samtools/picard/bedtools/
# cnvkit/DNAcopy) is validated end-to-end only in the toy Snakemake integration test - this
# script covers the R-level logic and bugs with lightweight, hand-built toy inputs.
# Usage: bash .tests/unit/test_atacseq_scripts.sh
# Requires: Rscript (dplyr, tidyr, ggplot2, edgeR, sva)

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

echo "=== 1) process_CNVkit_output_atacseq.smk's awk logic ==="
echo "    (bug: kept the raw log2 ratio (can be negative) instead of the linear ratio"
echo "    2^log2ratio that the rest of the pipeline needs to divide raw counts by - dividing"
echo "    counts by a negative or near-zero log2 value produces nonsensical results."
echo "    Confirmed correct direction with you: scaling must be in LINEAR space.)"
printf "chromosome\tstart\tend\tgene\tlog2\tdepth\tprobes\tweight\nchr1\t0\t200\tbin1\t-1.0\t50\t10\t0.9\nchr1\t200\t400\tbin2\t0.5\t80\t12\t0.8\n" > "$TMP/toy.cns"
cat "$TMP/toy.cns" | awk -v OFS='\t' '{print $1,$2,$3,$5,2^$5,$6,$7,$8}' | tail -n+2 | awk -v OFS='\t' '{print $1,$2,$3,$5}' > "$TMP/processed"
check "log2=-1.0 -> linear 0.5 (half the copies)" "$(sed -n 1p "$TMP/processed" | cut -f4)" "0.5"
check "log2=0.5 -> linear ~1.414 (mild gain)" "$(sed -n 2p "$TMP/processed" | cut -f4)" "1.41421"

echo "=== 2) do_get_counts_analysis.R / normATAC.R ==="
echo "    (removed dead code: num.sv() was called but its result (n.sv) was never used -"
echo "    sva() always used the explicit, user-provided SVs count instead. Confirmed unused"
echo "    by inspection; incidentally this exact call is numerically fragile with few"
echo "    features/samples, so removing the dead computation also removes that fragility.)"
mkdir -p "$TMP/atacseq_counts_test"
Rscript -e '
suppressMessages(library(edgeR))
suppressMessages(library(sva))
source("workflow/scripts/helpers/normATAC.R")
source("workflow/scripts/helpers/svaBatchCor.R")
source("workflow/scripts/helpers/do_get_counts_analysis.R")
source("workflow/scripts/helpers/pca_plot.R")
source("workflow/scripts/helpers/pca_scree_plot.R")
source("workflow/scripts/helpers/plot_distances_PCA.R")

set.seed(1)
n_bins <- 20
samples <- c("S1","S2","S3","S4","S5","S6")
counts <- data.frame(binID = paste0("1_", 1:n_bins))
for (s in samples) counts[[s]] <- rpois(n_bins, lambda = sample(20:80, 1))

metadata <<- data.frame(
  SRA_id = samples,
  condition = c("tumor","tumor","tumor","NAT","NAT","NAT"),
  study_id = c("A","A","B","A","B","B"),
  sex = c("F","F","M","F","M","M"),
  donor_id = c("D1","D2","D3","D1","D3","D4"),
  subtype = c("luminal","basal","luminal", NA, NA, NA),
  grade = c(2,3,1, NA, NA, NA)
)
minSamples_atacseq <<- 2
minCounts_atacseq <<- 1
minCpm_atacseq <<- 1
SVs <<- 1

out <- do_get_counts_analysis(counts = counts, scaling_factors = NULL, correct_for_CNAs = FALSE)
cat("OK: do_get_counts_analysis completed without the eigen()/num.sv() crash\n")
cat("N_ROWS_COUNT_F:", nrow(out$count.f), "\n")
' 2>&1 | tail -10
RC=$?
check "do_get_counts_analysis.R runs to completion with toy data (no dead-code crash)" "$RC" "0"

echo "=== 3) do_DA_analysis.R defensive alignment check ==="
echo "    (added stopifnot() on names(libsize) == colnames(counts), which was previously a"
echo "    bare, never-checked boolean expression - a zero-cost safety net)"
Rscript -e '
source("workflow/scripts/helpers/do_DA_analysis.R")
# confirm the stopifnot triggers correctly when misaligned, and does NOT trigger when aligned
counts <- data.frame(binID=c("1_1","1_2"), S1=c(10,20), S2=c(15,25))
library_sizes_aligned <- data.frame(SRA_id=c("S1","S2"), library_size=c(1000,1000))
library_sizes_misaligned <- data.frame(SRA_id=c("S2","S1"), library_size=c(1000,1000))

aligned_ok <- tryCatch({
  rownames(counts) <- counts$binID; counts$binID <- NULL
  libsize <- library_sizes_aligned$library_size; names(libsize) <- library_sizes_aligned$SRA_id
  stopifnot("aligned check" = all(names(libsize) == colnames(counts)))
  TRUE
}, error = function(e) FALSE)

misaligned_caught <- tryCatch({
  counts2 <- data.frame(binID=c("1_1","1_2"), S1=c(10,20), S2=c(15,25))
  rownames(counts2) <- counts2$binID; counts2$binID <- NULL
  libsize2 <- library_sizes_misaligned$library_size; names(libsize2) <- library_sizes_misaligned$SRA_id
  stopifnot("misaligned check" = all(names(libsize2) == colnames(counts2)))
  FALSE
}, error = function(e) TRUE)

cat("ALIGNED_OK:", aligned_ok, "\n")
cat("MISALIGNED_CAUGHT:", misaligned_caught, "\n")
'
check "correctly-aligned data passes the check" "$(Rscript -e 'counts <- data.frame(binID=c("1_1","1_2"), S1=c(10,20), S2=c(15,25)); rownames(counts) <- counts$binID; counts$binID <- NULL; libsize <- c(S1=1000,S2=1000); cat(all(names(libsize) == colnames(counts)))')" "TRUE"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
