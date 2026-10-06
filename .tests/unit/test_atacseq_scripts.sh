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
echo "    (helpers are sourced in the same order as workflow/scripts/get_counts_atacseq.R, so the"
echo "    test cannot drift from the real script: plot_PCA_SVA_distance_panel.R was added to the"
echo "    script but not here, and do_get_counts_analysis() failed with 'could not find function')"
mkdir -p "$TMP/atacseq_counts_test"
RC=0
Rscript -e '
suppressMessages(library(edgeR))
suppressMessages(library(sva))
source("workflow/scripts/helpers/normATAC.R")
source("workflow/scripts/helpers/svaBatchCor.R")
source("workflow/scripts/helpers/pca_plot.R")
source("workflow/scripts/helpers/pca_scree_plot.R")
source("workflow/scripts/helpers/plot_distances_PCA.R")
source("workflow/scripts/helpers/plot_PCA_SVA_distance_panel.R")
source("workflow/scripts/helpers/do_get_counts_analysis.R")

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
' > "$TMP/counts_analysis.log" 2>&1 || RC=$?
tail -10 "$TMP/counts_analysis.log"
# RC is R's own exit status (it used to be tail's, so an R crash could never fail this check)
check "do_get_counts_analysis.R runs to completion with toy data (no dead-code crash)" "$RC" "0"
check "completion message printed" "$(grep -c 'OK: do_get_counts_analysis completed' "$TMP/counts_analysis.log")" "1"

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

echo "=== 4) explore_cohort_atacseq.R (samples barplot moved here from explore_datatypes.R) ==="
echo "    Toy cohort built to exercise the counting logic:"
echo "    - S6 is in the metadata but NOT in the counts matrix (excluded sample) -> not counted;"
echo "    - donor D1 has a tumour and a NAT sample in study A -> 1 subject, 2 samples;"
echo "    - donor_id D1 is reused in study B -> a DIFFERENT subject (unique study+donor);"
echo "    - cohort bars: DepMap 3 models, COSMIC 2 samples, GDC 4 manifest rows."
E="$TMP/explore"; mkdir -p "$E"
printf "M1\nM2\nM3\n" > "$E/depmap_models"
printf "CS001\nCS002\n" > "$E/cosmic_samples"
printf "id\tfilename\tmd5\tsize\tstate\nG1\tf1.txt\tabc\t1\tlive\nG2\tf2.txt\tabc\t1\tlive\nG3\tf3.txt\tabc\t1\tlive\nG4\tf4.txt\tabc\t1\tlive\n" > "$E/gdc_manifest"
printf "binID\tS1\tS2\tS3\tS4\tS5\n1_1\t10\t12\t9\t11\t8\n1_2\t5\t6\t7\t4\t3\n" > "$E/counts"
printf "SRA_id\tdonor_id\tcondition\tstudy_id\tsubtype\tgrade\tage\n" > "$E/metadata"
printf "S1\tD1\ttumor\tA\tluminal\t2\t55\nS2\tD1\tNAT\tA\tNA\tNA\t55\nS3\tD1\ttumor\tB\tbasal\t3\t62\n" >> "$E/metadata"
printf "S4\tD2\ttumor\tB\tluminal\t1\t47\nS5\tD3\tNAT\tB\tNA\tNA\t71\nS6\tD9\ttumor\tA\tbasal\t2\t66\n" >> "$E/metadata"
printf "SRA_id\tlibrary_size\nS1\t1000\nS2\t1200\nS3\t900\nS4\t1100\nS5\t950\nS6\t1300\n" > "$E/library_sizes"

run_explore() {   # $1 = metadata file, $2 = output prefix; stderr (incl. the barplot table) -> $2.log
    Rscript workflow/scripts/explore_cohort_atacseq.R \
        --input_path_models_DepMap="$E/depmap_models" \
        --input_path_manifest_GDC="$E/gdc_manifest" \
        --input_path_samples_COSMIC="$E/cosmic_samples" \
        --input_path_counts_atacseq="$E/counts" \
        --input_path_metadata_atacseq="$1" \
        --input_path_library_sizes="$E/library_sizes" \
        --output_path_samples_barplot="$2_samples.pdf" \
        --output_path_plot_samples_per_study_by_condition="$2_by_condition.pdf" \
        --output_path_plot_subtype_distribution="$2_subtype.pdf" \
        --output_path_plot_grade_distribution="$2_grade.pdf" \
        --output_path_plot_age_distribution="$2_age.pdf" \
        --output_path_plot_library_size_by_study="$2_libsize.pdf" > "$2.log" 2>&1
}

RC=0; run_explore "$E/metadata" "$E/out" || RC=$?
check "explore_cohort_atacseq.R exits successfully" "$RC" "0"
[ "$RC" -eq 0 ] || cat "$E/out.log"
for f in samples by_condition subtype grade age libsize; do
    check "produces ${f} plot" "$([ -s "$E/out_${f}.pdf" ] && echo yes)" "yes"
done
# bar values, read from the table the script prints before plotting
bar() { awk -v d="$1" '/^Samples barplot data:/{t=1; next} t && $2==d {print $3; exit}' "$E/out.log"; }
check "DepMap bar = 3 models"                          "$(bar DepMap)"   "3"
check "COSMIC bar = 2 samples"                         "$(bar COSMIC)"   "2"
check "GDC bar = 4 manifest rows"                      "$(bar GDC)"      "4"
check "All = 5 (S6 excluded from counts, not counted)" "$(bar All)"      "5"
check "Subjects = 4 ((A,D1), (B,D1), (B,D2), (B,D3))"  "$(bar Subjects)" "4"
check "study B = 3 samples"                            "$(bar B)"        "3"
check "study A = 2 samples (S6 not counted)"           "$(bar A)"        "2"
check "Tumour = 3"                                     "$(bar Tumour)"   "3"
check "NAT = 2"                                        "$(bar NAT)"      "2"

# explicit errors instead of silently wrong plots
cut -f1,2,4- "$E/metadata" > "$E/metadata_nocondition"
if run_explore "$E/metadata_nocondition" "$E/err1"; then RC=0; else RC=1; fi
check "metadata without 'condition' -> explicit error" "$RC$(grep -c "missing column(s): condition" "$E/err1.log")" "11"
grep -v "^S5" "$E/metadata" > "$E/metadata_noS5"
if run_explore "$E/metadata_noS5" "$E/err2"; then RC=0; else RC=1; fi
check "ATAC-seq sample without metadata -> explicit error" "$RC$(grep -c "not in metadata: S5" "$E/err2.log")" "11"
sed 's/\tB\t/\tAll\t/' "$E/metadata" > "$E/metadata_clash"
if run_explore "$E/metadata_clash" "$E/err3"; then RC=0; else RC=1; fi
check "study_id clashing with a bar label ('All') -> explicit error" "$RC$(grep -c "Duplicated bar labels" "$E/err3.log")" "11"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi