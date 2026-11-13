#!/bin/bash
#SBATCH --job-name=snakemake_single
#SBATCH --partition=general
#SBATCH --output=/hpcnfs/scratch/MS/ieo7429/SelectionCaller/snakemake_%j.out
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ieo7429@ieo.it
# NOTE: the #SBATCH values above are only defaults - override them per run on the
# command line (sbatch --cpus-per-task=... --mem=... --time=... run_snakemake.sbatch ...),
# so the same launcher can be sized for each config.

set -uo pipefail
# NOTE: deliberately NOT using `set -e` - exit codes of the dry-run and main run are
# checked explicitly below.

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
USAGE="Usage: sbatch [sbatch opts] $0 <snakefile> <configfile> [cores] [target_rule]"
SNAKEFILE=${1:?"$USAGE"}
CONFIGFILE=${2:?"$USAGE"}
CORES="${3:-${SLURM_CPUS_PER_TASK:-1}}"   # default: exactly what SLURM gave us
TARGET_RULE="${4:-}"                       # e.g. "select_best_model" (--until); empty = rule all

# Rerun triggers: by default only missing/older inputs or a changed input set cause
# reruns, so edits to scripts/Snakefile/conda env do NOT rerun finished jobs.
# To restore Snakemake's full default:  RERUN_TRIGGERS="mtime params input software-env code"
# To force a rule that must be redone:  EXTRA_ARGS="--forcerun train_models"
RERUN_TRIGGERS="${RERUN_TRIGGERS:-mtime input}"
EXTRA_ARGS="${EXTRA_ARGS:-}"

LOCAL_BASE="/hpcnfs/scratch/MS/ieo7429/SelectionCaller"

# Unbuffered Python output: per-fold progress lines reach the logs in real time.
export PYTHONUNBUFFERED=1

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
SECONDS=0
log()     { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
section() { printf '\n=== %s ===\n' "$*"; }
hr()      { printf '%s\n' "--------------------------------------------------------------------"; }
elapsed_human() { local s=$SECONDS; printf '%02dh:%02dm:%02ds' $((s/3600)) $(((s%3600)/60)) $((s%60)); }

FINAL_STATUS=1
finish() { hr; log "Job finished (exit status: $FINAL_STATUS, elapsed: $(elapsed_human))"; hr; exit "$FINAL_STATUS"; }
trap finish EXIT

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
hr
log "Job ID:         ${SLURM_JOB_ID:-N/A}"
log "Job name:       ${SLURM_JOB_NAME:-N/A}"
log "Node:           ${SLURM_NODELIST:-N/A}"
log "CPUs assigned:  ${SLURM_CPUS_PER_TASK:-N/A} (--cores used: $CORES)"
log "Memory:         ${SLURM_MEM_PER_NODE:-N/A} MB"
log "Working dir:    $LOCAL_BASE"
log "Snakefile:      $SNAKEFILE"
log "Configfile:     $CONFIGFILE"
log "Target rule:    ${TARGET_RULE:-<none - running full 'rule all'>}"
log "Rerun triggers: $RERUN_TRIGGERS"
log "Extra args:     ${EXTRA_ARGS:-<none>}"
hr

if [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "$CORES" -gt "$SLURM_CPUS_PER_TASK" ]; then
    log "WARNING: --cores ($CORES) > CPUs assigned by SLURM ($SLURM_CPUS_PER_TASK): jobs will compete for CPUs."
fi

cd "$LOCAL_BASE" || { log "ERROR: cannot cd into $LOCAL_BASE"; FINAL_STATUS=1; exit; }

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
section "Activating environment"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-env
command -v snakemake > /dev/null 2>&1 || { log "ERROR: 'snakemake' not found after activating snakemake-env"; FINAL_STATUS=1; exit; }
log "snakemake: $(snakemake --version)"
[ -f "$SNAKEFILE" ]  || { log "ERROR: Snakefile not found: $SNAKEFILE";  FINAL_STATUS=1; exit; }
[ -f "$CONFIGFILE" ] || { log "ERROR: config file not found: $CONFIGFILE"; FINAL_STATUS=1; exit; }

# Snapshot of the config actually used by this job (traceability if the file is edited later)
mkdir -p logs/config_snapshots
cp "$CONFIGFILE" "logs/config_snapshots/$(basename "$CONFIGFILE" .yaml)_${SLURM_JOB_ID:-local}.yaml"
log "Config snapshot saved in logs/config_snapshots/"

# Common options: identical for dry-run and main run, so the dry-run shows exactly
# what the real run will do (previously the dry-run lacked --use-conda).
# $RERUN_TRIGGERS and $EXTRA_ARGS are intentionally unquoted (word splitting).
COMMON=(--snakefile "$SNAKEFILE" --configfile "$CONFIGFILE" --use-conda --rerun-incomplete)
UNTIL=()
[ -n "$TARGET_RULE" ] && UNTIL=(--until "$TARGET_RULE")

# ---------------------------------------------------------------------------
# Unlock + dry-run
# ---------------------------------------------------------------------------
section "Dry-run"
snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIGFILE" --unlock

log "Checking the DAG before committing any compute..."
snakemake "${COMMON[@]}" --rerun-triggers $RERUN_TRIGGERS $EXTRA_ARGS ${UNTIL[@]+"${UNTIL[@]}"} --dry-run
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
log "Launching Snakemake with $CORES cores..."
snakemake "${COMMON[@]}" --rerun-triggers $RERUN_TRIGGERS $EXTRA_ARGS ${UNTIL[@]+"${UNTIL[@]}"} \
    --cores "$CORES" \
    --keep-going \
    --printshellcmds
EXIT_STATUS_FINAL_RUN=$?

if [ "$EXIT_STATUS_FINAL_RUN" -eq 0 ]; then
    log "Pipeline completed successfully."
else
    log "ERROR: pipeline exited with status $EXIT_STATUS_FINAL_RUN."
fi
FINAL_STATUS=$EXIT_STATUS_FINAL_RUN