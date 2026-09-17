#!/bin/bash
#SBATCH --job-name=snakemake_single
#SBATCH --partition=general
#SBATCH --output=/hpcnfs/scratch/MS/ieo7429/SelectionCaller/snakemake_%j.out
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=600G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ieo7429@ieo.it

set -uo pipefail
# NOTE: deliberately NOT using `set -e` here - several steps below check $? explicitly
# (dry-run, then the main run), and `set -e` would abort the script the instant any of
# them fails, before we get a chance to inspect/report the exit code ourselves.

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
SNAKEFILE=${1:?"Usage: sbatch $0 <snakefile> <configfile> <cores>"}
CONFIGFILE=${2:?"Usage: sbatch $0 <snakefile> <configfile> <cores>"}
CORES=${3:?"Usage: sbatch $0 <snakefile> <configfile> <cores>"}

LOCAL_BASE="/hpcnfs/scratch/MS/ieo7429/SelectionCaller"

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
SECONDS=0   # bash builtin, used for elapsed-time reporting at the end

log()      { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
section()  { printf '\n=== %s ===\n' "$*"; }
hr()       { printf '%s\n' "--------------------------------------------------------------------"; }

elapsed_human() {
    local s=$SECONDS
    printf '%02dh:%02dm:%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60))
}

# Always print a closing summary, however the script exits (success, failure, or
# cancellation) - makes the .out log easy to skim without hunting for the exit code.
FINAL_STATUS=1
finish() {
    hr
    log "Job finished (exit status: $FINAL_STATUS, elapsed: $(elapsed_human))"
    hr
    exit "$FINAL_STATUS"
}
trap finish EXIT

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
hr
log "Job ID:        ${SLURM_JOB_ID:-N/A}"
log "Job name:      ${SLURM_JOB_NAME:-N/A}"
log "Node:          ${SLURM_NODELIST:-N/A}"
log "CPUs assigned: ${SLURM_CPUS_PER_TASK:-N/A} (--cores requested: $CORES)"
log "Memory:        ${SLURM_MEM_PER_NODE:-N/A} MB"
log "Working dir:   $LOCAL_BASE"
log "Snakefile:     $SNAKEFILE"
log "Configfile:    $CONFIGFILE"
hr

cd "$LOCAL_BASE" || { log "ERROR: cannot cd into $LOCAL_BASE"; FINAL_STATUS=1; exit; }

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
section "Activating environment"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-env

if ! command -v snakemake > /dev/null 2>&1; then
    log "ERROR: 'snakemake' not found on PATH after activating snakemake-env"
    FINAL_STATUS=1
    exit
fi
log "snakemake: $(snakemake --version)"

if [ ! -f "$SNAKEFILE" ]; then
    log "ERROR: Snakefile not found: $SNAKEFILE"
    FINAL_STATUS=1
    exit
fi
if [ ! -f "$CONFIGFILE" ]; then
    log "ERROR: config file not found: $CONFIGFILE"
    FINAL_STATUS=1
    exit
fi

# ---------------------------------------------------------------------------
# Unlock (in case a previous run left a stale lock) + dry-run
# ---------------------------------------------------------------------------
section "Dry-run"

snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIGFILE" --unlock

log "Checking the DAG before committing any compute..."
snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIGFILE" --rerun-incomplete --dry-run
EXIT_STATUS_DRY_RUN=$?

if [ "$EXIT_STATUS_DRY_RUN" -ne 0 ]; then
    log "ERROR: dry-run failed (exit status $EXIT_STATUS_DRY_RUN) - not starting the real run."
    FINAL_STATUS=$EXIT_STATUS_DRY_RUN
    exit
fi
log "Dry-run OK."

# ---------------------------------------------------------------------------
# Main run
# ---------------------------------------------------------------------------
section "Main pipeline run"

snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIGFILE" --unlock

log "Launching Snakemake with $CORES cores..."
snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    --cores "$CORES" \
    --use-conda \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds
EXIT_STATUS_FINAL_RUN=$?

if [ "$EXIT_STATUS_FINAL_RUN" -eq 0 ]; then
    log "Pipeline completed successfully."
else
    log "ERROR: pipeline exited with status $EXIT_STATUS_FINAL_RUN."
fi

FINAL_STATUS=$EXIT_STATUS_FINAL_RUN
