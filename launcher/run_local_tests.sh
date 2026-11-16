#!/usr/bin/env bash
# run_local_tests.sh - run the tests and Snakemake (dry-runs / toy runs) on the WORKSTATION, without SLURM.
#
# Purpose
#   Check code changes locally before pushing them and running them on the cluster.
#   There is no real data on the workstation, so Snakemake can only be run here as a
#   dry-run (-n) or on the toy fixtures of .tests/integration/.
#   launcher/launcher.sh is the CLUSTER launcher (SLURM, cluster paths): it is not used here.
#
#   Everything runs in a TEMPORARY COPY of the code (mktemp -d), never in the repository:
#   Snakemake creates resources/, data/, results/, logs/ and .snakemake/ at parse time, and the
#   dry-runs need empty placeholder files in resources/. None of these must appear in
#   ~/Desktop/SelectionCallerCodeBase (see CLAUDE.md, section 2).
#
# Modes (one or more; default: unit dryrun)
#   unit     every .tests/unit/*.sh (toy data, no Snakemake DAG except where a test runs one itself)
#   dryrun   snakemake -n on:
#              - the integration toy configs (.tests/integration/config/config_dataset.yaml and
#                config_ML.yaml), with the toy fixtures copied into resources/ as in CI;
#              - every production dataset config (config/config_dataset*_*.yaml), with empty
#                placeholders for the source files that no rule produces (as in CI).
#              config/config_ML.yaml cannot be dry-run here (it needs the real feature matrix).
#   toy-ml   real (non-dry) run of workflow/snakefile_ML on the toy ML config (a few minutes).
#
# Options
#   --keep       do not delete the temporary copy at the end (its path is printed), to inspect outputs/logs
#   --cores N    cores for toy-ml (default 2)
#
# Inputs   the code of this repository (workflow/, config/, .tests/, launcher/)
# Outputs  a pass/fail summary on stdout; per-step logs in <tmp>/_local_test_logs/ (kept with --keep or on failure)
#
# Usage
#   bash launcher/run_local_tests.sh                    # unit + dryrun
#   bash launcher/run_local_tests.sh dryrun             # only the dry-runs
#   bash launcher/run_local_tests.sh unit dryrun toy-ml --keep
#
# Requires  ~/miniforge3 with the conda env snakemake-env (same name as on the cluster), Rscript,
#           and the tools/packages each unit test declares in its own header.

set -uo pipefail
# NOTE: deliberately NOT using `set -e`: every step is run, its exit status recorded, and a
# summary printed at the end.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
MODES=()
KEEP_TMP=0
TOY_CORES=2
while [ $# -gt 0 ]; do
    case "$1" in
        unit|dryrun|toy-ml) MODES+=("$1"); shift ;;
        --keep)             KEEP_TMP=1; shift ;;
        --cores)            TOY_CORES="${2:?--cores needs a value}"; shift 2 ;;
        -h|--help)          sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "ERROR: unknown argument '$1' (valid: unit, dryrun, toy-ml, --keep, --cores N; see --help)"; exit 2 ;;
    esac
done
[ ${#MODES[@]} -eq 0 ] && MODES=(unit dryrun)

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
if [ ! -f ~/miniforge3/etc/profile.d/conda.sh ]; then
    echo "ERROR: ~/miniforge3/etc/profile.d/conda.sh not found: cannot activate the snakemake-env conda env."
    exit 1
fi
source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-env || { echo "ERROR: cannot activate conda env 'snakemake-env' (create it from launcher/snakemake_env.yml)."; exit 1; }
command -v snakemake > /dev/null 2>&1 || { echo "ERROR: 'snakemake' not found after activating snakemake-env."; exit 1; }

# ---------------------------------------------------------------------------
# Temporary copy of the code (never run inside the repository)
# ---------------------------------------------------------------------------
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/selectioncaller_local_tests.XXXXXX")"
# Guard: an empty WORK_DIR would make the rsync below write into "/".
if [ -z "$WORK_DIR" ] || [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: could not create a temporary directory in '${TMPDIR:-/tmp}' (does it exist and is it writable?)."
    exit 1
fi
LOG_DIR="$WORK_DIR/_local_test_logs"
mkdir -p "$LOG_DIR"
# Only code is copied: no .git, no _samples (format examples, must not be used as data),
# no _from_cluster (outputs fetched from the cluster), no pipeline output folders.
rsync -a \
    --exclude '.git/' --exclude '_samples/' --exclude '_from_cluster/' \
    --exclude '/data/' --exclude '/resources/' --exclude '/results/' --exclude '/analyses/' \
    --exclude '/logs/' --exclude '/logs_pilot/' --exclude '.snakemake/' --exclude '__pycache__/' \
    "$REPO_ROOT/" "$WORK_DIR/"
cd "$WORK_DIR" || { echo "ERROR: cannot cd into $WORK_DIR"; exit 1; }

# ---------------------------------------------------------------------------
# Step runner and summary
# ---------------------------------------------------------------------------
STEP_NAMES=()
STEP_STATUS=()
run_step() {
    # run_step <name> <command...>: runs the command in $WORK_DIR, output to its own log file
    local name="$1"; shift
    local log_file="$LOG_DIR/$(echo "$name" | tr ' /:' '___').log"
    printf '  %-70s ' "$name"
    ( cd "$WORK_DIR" && "$@" ) > "$log_file" 2>&1
    local status=$?
    STEP_NAMES+=("$name")
    if [ "$status" -eq 0 ]; then
        STEP_STATUS+=("PASS"); echo "PASS"
    else
        STEP_STATUS+=("FAIL"); echo "FAIL (exit $status, log: $log_file)"
        tail -n 15 "$log_file" | sed 's/^/      | /'
    fi
}

# Same list of empty placeholders as the CI job "Dry-run (production config)"
# (.github/workflows/ci.yml): source files that no rule produces must exist for --dry-run,
# even though their content is never read.
make_production_placeholders() {
    mkdir -p resources/DepMap/models resources/DepMap/probabilities resources/DepMap/scores
    touch resources/DepMap/models/Model.csv
    touch resources/DepMap/probabilities/CRISPRGeneDependency.csv
    touch resources/DepMap/scores/CRISPRGeneEffect.csv
    mkdir -p resources/IntOGen/IntOGen-Drivers-20240920
    touch resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv
    mkdir -p resources/COSMIC/classification/Cosmic_Classification_Tsv_v103_GRCh38
    mkdir -p resources/COSMIC/samples/Cosmic_Sample_Tsv_v103_GRCh38
    mkdir -p resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38
    mkdir -p resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38
    touch resources/COSMIC/classification/Cosmic_Classification_Tsv_v103_GRCh38/Cosmic_Classification_v103_GRCh38.tsv
    touch resources/COSMIC/samples/Cosmic_Sample_Tsv_v103_GRCh38/Cosmic_Sample_v103_GRCh38.tsv
    touch resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38/Cosmic_GenomeScreensMutant_v103_GRCh38.tsv
    touch resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38/Cosmic_NonCodingVariants_v103_GRCh38.tsv
    mkdir -p resources/GDC/CNAs/clinical.cohort.2026-01-28
    touch resources/GDC/CNAs/clinical.cohort.2026-01-28/clinical.tsv
    touch resources/GDC/CNAs/gdc_sample_sheet.2026-01-28.tsv
    touch resources/GDC/CNAs/gdc_manifest.2026-01-28.170055.txt
    mkdir -p resources/GDC/CNAs/files/dummy
    touch resources/GDC/CNAs/files/dummy/dummy.seg.v2.txt
    mkdir -p resources/fragile_sites/origin resources/ohnologs_and_paralogs/origin resources/haploinsufficiency resources/epigenomes
    touch resources/fragile_sites/origin/fragile_sites_li_hg38.bed
    touch resources/ohnologs_and_paralogs/origin/OhnologList.RData resources/ohnologs_and_paralogs/origin/ParalogList.RData
    mkdir -p resources/protein_interactions/origin
    touch resources/protein_interactions/origin/HumanPCgenes_Subunits_allchr.RData resources/protein_interactions/origin/H_sapiens_interfacesHQ_genename.RData
    touch resources/haploinsufficiency/haploinsufficiency.csv
    mkdir -p resources/atacseq/files/coverage resources/atacseq/files/cnvkit
    touch resources/epigenomes/mapping_EID_tt_raw.tsv
}

# Each dry-run gets its own clean copy of the code, so placeholders / toy fixtures / the
# folders Snakemake creates at parse time cannot leak from one dry-run into another.
fresh_copy() {
    local dest="$1"
    rsync -a --exclude '_local_test_logs/' --exclude '_dryrun_*/' "$WORK_DIR/" "$dest/"
}

dryrun_production() {
    local snakefile="$1" config="$2" dest="$3"
    fresh_copy "$dest" && cd "$dest" && make_production_placeholders \
        && snakemake -n -s "$snakefile" --configfile "$config" --rerun-triggers mtime input
}

dryrun_toy_dataset() {
    local dest="$1"
    # Toy fixtures copied into the paths expected by the rules, as in the CI integration job
    fresh_copy "$dest" && cd "$dest" && mkdir -p resources \
        && cp -r .tests/integration/data/resources/. ./resources/ \
        && Rscript .tests/integration/make_toy_GTEx_interactor_inputs.R \
               resources/reference_data/gencode/gencode.v1.annotation.gff3.gz config/mapping_GTEx_TCGA resources \
        && snakemake -n -s workflow/snakefile_dataset_bin_level \
               --configfile .tests/integration/config/config_dataset.yaml --rerun-triggers mtime input
}

dryrun_toy_ml() {
    local dest="$1"
    fresh_copy "$dest" && cd "$dest" \
        && snakemake -n -s workflow/snakefile_ML --configfile .tests/integration/config/config_ML.yaml \
               --rerun-triggers mtime input
}

run_toy_ml() {
    local dest="$1"
    fresh_copy "$dest" && cd "$dest" \
        && snakemake -s workflow/snakefile_ML --configfile .tests/integration/config/config_ML.yaml \
               --cores "$TOY_CORES" -p
}

# ---------------------------------------------------------------------------
# Run the selected modes
# ---------------------------------------------------------------------------
echo "Repository:      $REPO_ROOT"
echo "Temporary copy:  $WORK_DIR"
echo "Snakemake:       $(snakemake --version)"
echo "Modes:           ${MODES[*]}"

for mode in "${MODES[@]}"; do
    case "$mode" in
        unit)
            echo; echo "=== Unit tests (.tests/unit/) ==="
            for test_script in .tests/unit/*.sh; do
                run_step "unit: $(basename "$test_script")" bash "$test_script"
            done
            ;;
        dryrun)
            echo; echo "=== Snakemake dry-runs ==="
            run_step "dry-run: toy dataset config"  dryrun_toy_dataset "$WORK_DIR/_dryrun_toy_dataset"
            run_step "dry-run: toy ML config"       dryrun_toy_ml      "$WORK_DIR/_dryrun_toy_ml"
            for config in config/config_dataset*_*.yaml; do
                tag="$(basename "$config" .yaml)"
                run_step "dry-run: $config" dryrun_production workflow/snakefile_dataset_bin_level "$config" "$WORK_DIR/_dryrun_$tag"
            done
            # config/config_ML.yaml is not dry-run here: snakefile_ML refuses at parse time when the
            # configured feature matrix does not exist, and real matrices exist only on the cluster.
            ;;
        toy-ml)
            echo; echo "=== Toy run of snakefile_ML (real execution, $TOY_CORES cores) ==="
            run_step "toy run: snakefile_ML" run_toy_ml "$WORK_DIR/_toyrun_ML"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
N_FAIL=0
for status in "${STEP_STATUS[@]}"; do [ "$status" = "FAIL" ] && N_FAIL=$((N_FAIL + 1)); done
echo
echo "=== Summary: $(( ${#STEP_STATUS[@]} - N_FAIL )) passed, $N_FAIL failed (of ${#STEP_STATUS[@]}) ==="
for i in "${!STEP_NAMES[@]}"; do
    [ "${STEP_STATUS[$i]}" = "FAIL" ] && echo "  FAIL  ${STEP_NAMES[$i]}"
done

if [ "$KEEP_TMP" -eq 1 ] || [ "$N_FAIL" -gt 0 ]; then
    echo "Temporary copy kept (logs in $LOG_DIR). Delete it when done: rm -rf $WORK_DIR"
else
    rm -rf "$WORK_DIR"
fi
[ "$N_FAIL" -eq 0 ]
