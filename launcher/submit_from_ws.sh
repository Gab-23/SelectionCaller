#!/usr/bin/env bash
# submit_from_ws.sh - submit cluster jobs FROM THE WORKSTATION (ssh + sbatch on login.ieo.it).
#
# Purpose
#   Launch, from ~/Desktop/SelectionCallerCodeBase, either
#     - Snakemake through the existing cluster launcher (launcher/launcher.sh, unchanged), or
#     - a pilot/analysis script through `sbatch --wrap` with the pilot conda env,
#   on the cluster checkout /hpcnfs/scratch/MS/ieo7429/SelectionCaller.
#   It never modifies files on the cluster: it only reads git state, checks the queue and calls sbatch.
#
# Safety checks before submitting (it refuses if one fails)
#   1. the cluster checkout is at the same commit as the local HEAD (git rev-parse HEAD on both sides);
#   2. the cluster working tree has no modified tracked files;
#   (a local working tree with uncommitted changes only gives a WARNING: those changes are NOT on the cluster)
#   3. Snakemake mode only: never two Snakemake jobs on the same working directory at the same time.
#      If a job named snakemake* of this user is pending/running, the new job is chained to it with
#      --dependency=afterok:<id> (if the previous job fails, the new one stays pending with reason
#      DependencyNeverSatisfied: cancel it with scancel).
#
# Usage
#   launcher/submit_from_ws.sh snakemake <snakefile> <configfile> [cores] [target_rule] [options]
#   launcher/submit_from_ws.sh pilot "<command run in the project folder>" [options]
#
#   Options (both modes)
#     --after <jobid>          add --dependency=afterok:<jobid> (repeatable)
#     --sbatch "<opts>"        extra sbatch options, e.g. "--mem=200G --time=72:00:00 --partition=medium"
#     --dry                    run every check and print the sbatch command, but do not submit
#   Options (snakemake mode; exported to launcher/launcher.sh)
#     --rerun-triggers "<t>"   RERUN_TRIGGERS (launcher default: "mtime input")
#     --extra-args "<args>"    EXTRA_ARGS, e.g. "--forcerun train_models"
#   Options (pilot mode)
#     --job-name <name>        SLURM job name (default: pilot); log in logs_pilot/<name>_<jobid>.out
#
#   If [cores] is given and --sbatch does not set --cpus-per-task, --cpus-per-task=<cores> is added
#   (launcher.sh's own default is 20; --cores must match the CPUs of the job).
#
# Examples
#   launcher/submit_from_ws.sh snakemake workflow/snakefile_dataset_bin_level config/config_dataset_LUAD.yaml 20 combine_dataset
#   launcher/submit_from_ws.sh pilot "python pilot_omics_background.py --help" --job-name pilot_LUAD --sbatch "--cpus-per-task=8 --mem=64G"
#
# Outputs  the job IDs and the path of the job logs on the cluster (fetch them with launcher/fetch_from_cluster.sh).
#
# Environment overrides (testing only): CLUSTER_HOST, CLUSTER_DIR.
# The ssh password is asked once: one shared ssh connection (ControlMaster) is reused for every step.

set -euo pipefail

CLUSTER_HOST="${CLUSTER_HOST:-ieo7429@login.ieo.it}"
CLUSTER_DIR="${CLUSTER_DIR:-/hpcnfs/scratch/MS/ieo7429/SelectionCaller}"
CONDA_SH='~/miniforge3/etc/profile.d/conda.sh'
PILOT_ENV="$CLUSTER_DIR/.snakemake/conda/56a4e7b4dc24521ebcee3d03f716b893_"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

die() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
USAGE="Usage: $0 snakemake <snakefile> <configfile> [cores] [target_rule] [options]
       $0 pilot \"<command>\" [options]          (see the header of this script)"
[ $# -ge 1 ] || die "$USAGE"
MODE="$1"; shift

POSITIONAL=()
AFTER_IDS=()
SBATCH_EXTRA=""
DRY=0
RERUN_TRIGGERS_ARG=""
EXTRA_ARGS_ARG=""
PILOT_JOB_NAME="pilot"
while [ $# -gt 0 ]; do
    case "$1" in
        --after)          AFTER_IDS+=("${2:?--after needs a job id}"); shift 2 ;;
        --sbatch)         SBATCH_EXTRA="${2?--sbatch needs a value}"; shift 2 ;;
        --dry)            DRY=1; shift ;;
        --rerun-triggers) RERUN_TRIGGERS_ARG="${2:?--rerun-triggers needs a value}"; shift 2 ;;
        --extra-args)     EXTRA_ARGS_ARG="${2?--extra-args needs a value}"; shift 2 ;;
        --job-name)       PILOT_JOB_NAME="${2:?--job-name needs a value}"; shift 2 ;;
        --*)              die "unknown option '$1'. $USAGE" ;;
        *)                POSITIONAL+=("$1"); shift ;;
    esac
done
for job_id in ${AFTER_IDS[@]+"${AFTER_IDS[@]}"}; do
    [[ "$job_id" =~ ^[0-9]+$ ]] || die "--after expects a numeric SLURM job id, got '$job_id'."
done

case "$MODE" in
    snakemake)
        [ ${#POSITIONAL[@]} -ge 2 ] && [ ${#POSITIONAL[@]} -le 4 ] || die "$USAGE"
        SNAKEFILE="${POSITIONAL[0]}"
        CONFIGFILE="${POSITIONAL[1]}"
        CORES="${POSITIONAL[2]:-}"
        TARGET_RULE="${POSITIONAL[3]:-}"
        [ -z "$CORES" ] || [[ "$CORES" =~ ^[0-9]+$ ]] || die "cores must be an integer, got '$CORES'."
        # The paths must exist in the local checkout (same commit as the cluster, checked below)
        [ -f "$REPO_ROOT/$SNAKEFILE" ]  || die "snakefile '$SNAKEFILE' not found in $REPO_ROOT (give it relative to the project folder)."
        [ -f "$REPO_ROOT/$CONFIGFILE" ] || die "config file '$CONFIGFILE' not found in $REPO_ROOT (give it relative to the project folder)."
        [ "$PILOT_JOB_NAME" = "pilot" ] || die "--job-name is for pilot mode only (snakemake jobs are named snakemake_<target>_<config>)."
        ;;
    pilot)
        [ ${#POSITIONAL[@]} -eq 1 ] || die "pilot mode takes exactly one quoted command. $USAGE"
        PILOT_COMMAND="${POSITIONAL[0]}"
        [ -z "$RERUN_TRIGGERS_ARG$EXTRA_ARGS_ARG" ] || die "--rerun-triggers / --extra-args are for snakemake mode only."
        [[ "$PILOT_JOB_NAME" =~ ^[A-Za-z0-9_.-]+$ ]] || die "--job-name may contain only letters, digits, '_', '.', '-' (got '$PILOT_JOB_NAME')."
        ;;
    *) die "first argument must be 'snakemake' or 'pilot'. $USAGE" ;;
esac

# ---------------------------------------------------------------------------
# One shared ssh connection (password typed once)
# ---------------------------------------------------------------------------
SSH_CONTROL_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sc_ssh.XXXXXX")"
[ -n "$SSH_CONTROL_DIR" ] && [ -d "$SSH_CONTROL_DIR" ] || die "cannot create a temporary directory for the ssh control socket."
SSH_OPTS=(-o ControlMaster=auto -o "ControlPath=$SSH_CONTROL_DIR/%r@%h:%p" -o ControlPersist=300)
cleanup() {
    ssh "${SSH_OPTS[@]}" -O exit "$CLUSTER_HOST" > /dev/null 2>&1 || true
    rm -rf "$SSH_CONTROL_DIR"
}
trap cleanup EXIT

# remote <script>: runs a bash script on the cluster (script sent on stdin, so no quoting issues)
remote() { ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" 'bash -s' <<< "$1"; }

# ---------------------------------------------------------------------------
# Check 1-2: same commit, clean cluster working tree
# ---------------------------------------------------------------------------
LOCAL_HEAD="$(git -C "$REPO_ROOT" rev-parse HEAD)"
LOCAL_DIRTY="$(git -C "$REPO_ROOT" status --porcelain --untracked-files=no)"
if [ -n "$LOCAL_DIRTY" ]; then
    echo "WARNING: the local working tree has uncommitted changes to tracked files; they are NOT on the cluster:"
    echo "$LOCAL_DIRTY" | sed 's/^/    /'
fi

echo "Connecting to $CLUSTER_HOST ..."
REMOTE_STATE="$(remote "cd $(printf '%q' "$CLUSTER_DIR") || exit 3
echo \"HEAD \$(git rev-parse HEAD)\"
git status --porcelain --untracked-files=no | sed 's/^/DIRTY /'
squeue -h -u \"\$USER\" -t PENDING,RUNNING -o 'JOB %i %j %T' 2>/dev/null
")" || die "could not read the state of $CLUSTER_DIR on $CLUSTER_HOST (ssh failed or folder missing)."

REMOTE_HEAD="$(echo "$REMOTE_STATE" | awk '$1 == "HEAD" {print $2}')"
REMOTE_DIRTY="$(echo "$REMOTE_STATE" | sed -n 's/^DIRTY //p')"
[ -n "$REMOTE_HEAD" ] || die "could not read git HEAD of the cluster checkout. Output was:
$REMOTE_STATE"

if [ "$REMOTE_HEAD" != "$LOCAL_HEAD" ]; then
    die "the cluster checkout is at ${REMOTE_HEAD:0:7}, the local HEAD is ${LOCAL_HEAD:0:7}.
       Push the local commits (git push), then on the cluster run:
         cd $CLUSTER_DIR && git pull
       and submit again."
fi
if [ -n "$REMOTE_DIRTY" ]; then
    die "the cluster working tree has modified tracked files (code must not be edited on the cluster):
$(echo "$REMOTE_DIRTY" | sed 's/^/    /')
       Resolve them on the cluster (e.g. move the edits into the workstation copy, push, pull), then submit again."
fi
echo "OK: cluster checkout at ${LOCAL_HEAD:0:7} (same as local HEAD), no modified tracked files."

# ---------------------------------------------------------------------------
# Check 3 (snakemake only): chain after snakemake jobs already queued/running
# ---------------------------------------------------------------------------
DEPENDENCY_IDS=(${AFTER_IDS[@]+"${AFTER_IDS[@]}"})
if [ "$MODE" = "snakemake" ]; then
    QUEUED_SNAKEMAKE="$(echo "$REMOTE_STATE" | awk '$1 == "JOB" && $3 ~ /^snakemake/ {print $2" "$3" "$4}')"
    if [ -n "$QUEUED_SNAKEMAKE" ]; then
        echo "Snakemake jobs already queued/running in this working directory (the new job will wait for them, afterok):"
        echo "$QUEUED_SNAKEMAKE" | sed 's/^/    /'
        while read -r job_id _; do DEPENDENCY_IDS+=("$job_id"); done <<< "$QUEUED_SNAKEMAKE"
    fi
fi

# ---------------------------------------------------------------------------
# Build the sbatch command
# ---------------------------------------------------------------------------
SBATCH_ARGS=(--parsable)
if [ ${#DEPENDENCY_IDS[@]} -gt 0 ]; then
    SBATCH_ARGS+=("--dependency=afterok:$(IFS=:; echo "${DEPENDENCY_IDS[*]}")")
fi
# $SBATCH_EXTRA is intentionally word-split (several options in one string)
read -r -a SBATCH_EXTRA_ARRAY <<< "$SBATCH_EXTRA"
SBATCH_ARGS+=(${SBATCH_EXTRA_ARRAY[@]+"${SBATCH_EXTRA_ARRAY[@]}"})

if [ "$MODE" = "snakemake" ]; then
    JOB_NAME="snakemake_${TARGET_RULE:-all}_$(basename "$CONFIGFILE" .yaml)"
    SBATCH_ARGS+=("--job-name=$JOB_NAME")
    if [ -n "$CORES" ] && [[ " $SBATCH_EXTRA " != *"--cpus-per-task"* ]] && [[ " $SBATCH_EXTRA " != *" -c "* ]]; then
        SBATCH_ARGS+=("--cpus-per-task=$CORES")
    fi
    # Environment variables read by launcher.sh (sbatch exports the submitting environment by default)
    ENV_PREFIX=""
    [ -n "$RERUN_TRIGGERS_ARG" ] && ENV_PREFIX+="RERUN_TRIGGERS=$(printf '%q' "$RERUN_TRIGGERS_ARG") "
    [ -n "$EXTRA_ARGS_ARG" ]     && ENV_PREFIX+="EXTRA_ARGS=$(printf '%q' "$EXTRA_ARGS_ARG") "
    LAUNCHER_ARGS=("$SNAKEFILE" "$CONFIGFILE")
    [ -n "$CORES" ] && LAUNCHER_ARGS+=("$CORES")
    [ -n "$TARGET_RULE" ] && { [ -n "$CORES" ] || die "give [cores] when giving [target_rule] (positional arguments)."; LAUNCHER_ARGS+=("$TARGET_RULE"); }
    SUBMIT_LINE="${ENV_PREFIX}sbatch $(printf '%q ' "${SBATCH_ARGS[@]}") launcher/launcher.sh $(printf '%q ' "${LAUNCHER_ARGS[@]}")"
    LOG_HINT="$CLUSTER_DIR/snakemake_<jobid>.out (launcher.sh default, unless --output was given in --sbatch)"
else
    # logs_pilot/ must already exist on the cluster: this script does not create folders there
    SBATCH_ARGS+=("--job-name=$PILOT_JOB_NAME" "--output=$CLUSTER_DIR/logs_pilot/%x_%j.out")
    WRAPPED=". $CONDA_SH && conda activate $(printf '%q' "$PILOT_ENV") && cd $(printf '%q' "$CLUSTER_DIR") && export PYTHONUNBUFFERED=1 && $PILOT_COMMAND"
    SUBMIT_LINE="[ -d logs_pilot ] || { echo 'ERROR: logs_pilot/ does not exist in $CLUSTER_DIR' >&2; exit 4; }
sbatch $(printf '%q ' "${SBATCH_ARGS[@]}") --wrap=$(printf '%q' "$WRAPPED")"
    LOG_HINT="$CLUSTER_DIR/logs_pilot/${PILOT_JOB_NAME}_<jobid>.out"
fi

echo
echo "Command on the cluster (in $CLUSTER_DIR):"
echo "$SUBMIT_LINE" | sed 's/^/    /'
if [ "$DRY" -eq 1 ]; then
    echo "--dry given: nothing submitted."
    exit 0
fi

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
SUBMIT_OUTPUT="$(remote "cd $(printf '%q' "$CLUSTER_DIR") || exit 3
$SUBMIT_LINE")" || die "sbatch failed on the cluster. Output:
$SUBMIT_OUTPUT"
JOB_ID="$(echo "$SUBMIT_OUTPUT" | tail -n 1 | cut -d';' -f1)"
[[ "$JOB_ID" =~ ^[0-9]+$ ]] || die "unexpected sbatch output (no job id):
$SUBMIT_OUTPUT"

echo
echo "Submitted job $JOB_ID${DEPENDENCY_IDS[0]+ (after: ${DEPENDENCY_IDS[*]})}"
echo "Log on the cluster: ${LOG_HINT//<jobid>/$JOB_ID}"
echo "Follow it with:"
echo "    launcher/cluster_status.sh --watch          (queue, progress, finished jobs)"
echo "    launcher/cluster_status.sh --job $JOB_ID     (accounting and log tail of this job)"
echo "Wait for it and fetch outputs with:"
echo "    launcher/fetch_from_cluster.sh --wait $JOB_ID <paths relative to the project folder>"
echo "JOB_ID=$JOB_ID"
