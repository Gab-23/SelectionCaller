#!/usr/bin/env bash
# fetch_from_cluster.sh - copy selected cluster outputs to the workstation for analysis.
#
# Purpose
#   Optionally wait for a SLURM job to end, then copy the requested files/folders (given relative
#   to the project folder /hpcnfs/scratch/MS/ieo7429/SelectionCaller) into
#   ~/Desktop/SelectionCallerCodeBase/_from_cluster/, keeping the same relative paths
#   (e.g. analyses/v4/LUAD/1000000bp/x.tsv -> _from_cluster/analyses/v4/LUAD/1000000bp/x.tsv).
#   _from_cluster/ is ignored by git (.git/info/exclude, like _samples/) and is NOT the pipeline's
#   data/ or results/: nothing here is read by Snakemake.
#
# Safety (the workstation disk is almost full)
#   - before copying, the total size on the cluster is measured (du) and the copy is refused if it is
#     larger than --max-size (default 1G) or would leave less than 20 GB free on the local disk;
#   - raw/binary files are excluded unless --include-binary: *.bam *.bai *.cram *.sra *.fastq.gz *.fq.gz
#     *.cool *.mcool *.hic *.pairs.gz *.pairsam.gz;
#   - only paths inside the project folder are accepted (no absolute paths, no '..').
#   It never modifies anything on the cluster (read-only: squeue, sacct, du, rsync as sender).
#
# Usage
#   launcher/fetch_from_cluster.sh [options] <path relative to the project folder> [...]
#   Options
#     --wait <jobid>      wait until the job leaves the queue (checked every 60 s), print its final state,
#                         and also fetch its log (snakemake_<jobid>.out, logs_pilot/*_<jobid>.out)
#     --max-size <size>   maximum total size, with suffix K/M/G (default 1G)
#     --include-binary    do not exclude BAM/SRA/FASTQ/cool/hic/pairs files
#     --dry               measure and list what would be copied, copy nothing
#
# Examples
#   launcher/fetch_from_cluster.sh --wait 123456 analyses/v4/LUAD/1000000bp
#   launcher/fetch_from_cluster.sh --dry results/ML
#
# Outputs  files under _from_cluster/<same relative path>; existing local copies are updated (rsync).
#          They are temporary working copies: remove them after the analysis with launcher/clean_fetched.sh
#          (the originals on the cluster are never touched).
# Environment overrides (testing only): CLUSTER_HOST, CLUSTER_DIR.

set -euo pipefail
# Numbers are printed and parsed with a decimal point whatever the workstation locale (it_IT uses ',').
export LC_ALL=C

CLUSTER_HOST="${CLUSTER_HOST:-ieo7429@login.ieo.it}"
CLUSTER_DIR="${CLUSTER_DIR:-/hpcnfs/scratch/MS/ieo7429/SelectionCaller}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_DEST="$REPO_ROOT/_from_cluster"
MIN_FREE_LOCAL_BYTES=$((20 * 1024 * 1024 * 1024))
POLL_SECONDS=60
BINARY_PATTERNS=('*.bam' '*.bai' '*.cram' '*.sra' '*.fastq.gz' '*.fq.gz' '*.cool' '*.mcool' '*.hic' '*.pairs.gz' '*.pairsam.gz')

die() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
USAGE="Usage: $0 [--wait <jobid>] [--max-size 1G] [--include-binary] [--dry] <path relative to the project folder> [...]"
WAIT_JOB=""
MAX_SIZE="1G"
INCLUDE_BINARY=0
DRY=0
REMOTE_PATHS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --wait)           WAIT_JOB="${2:?--wait needs a job id}"; shift 2 ;;
        --max-size)       MAX_SIZE="${2:?--max-size needs a value}"; shift 2 ;;
        --include-binary) INCLUDE_BINARY=1; shift ;;
        --dry)            DRY=1; shift ;;
        -h|--help)        sed -n '2,36p' "$0"; exit 0 ;;
        --*)              die "unknown option '$1'. $USAGE" ;;
        *)                REMOTE_PATHS+=("$1"); shift ;;
    esac
done
[ -z "$WAIT_JOB" ] || [[ "$WAIT_JOB" =~ ^[0-9]+$ ]] || die "--wait expects a numeric SLURM job id, got '$WAIT_JOB'."
[ ${#REMOTE_PATHS[@]} -gt 0 ] || [ -n "$WAIT_JOB" ] || die "nothing to fetch. $USAGE"
for remote_path in ${REMOTE_PATHS[@]+"${REMOTE_PATHS[@]}"}; do
    case "$remote_path" in
        /*)                      die "'$remote_path' is absolute: give paths relative to $CLUSTER_DIR." ;;
        ..|../*|*/..|*/../*)     die "'$remote_path' contains '..': only paths inside $CLUSTER_DIR are allowed." ;;
    esac
done
# Size limit in bytes
if [[ "$MAX_SIZE" =~ ^([0-9]+)([KMG]?)$ ]]; then
    case "${BASH_REMATCH[2]}" in
        K) MAX_BYTES=$(( BASH_REMATCH[1] * 1024 )) ;;
        M) MAX_BYTES=$(( BASH_REMATCH[1] * 1024 * 1024 )) ;;
        G) MAX_BYTES=$(( BASH_REMATCH[1] * 1024 * 1024 * 1024 )) ;;
        *) MAX_BYTES=${BASH_REMATCH[1]} ;;
    esac
else
    die "--max-size must be an integer with optional K/M/G suffix (e.g. 500M, 2G), got '$MAX_SIZE'."
fi

# ---------------------------------------------------------------------------
# One shared ssh connection (password typed once)
# ---------------------------------------------------------------------------
SSH_CONTROL_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sc_ssh.XXXXXX")"
[ -n "$SSH_CONTROL_DIR" ] && [ -d "$SSH_CONTROL_DIR" ] || die "cannot create a temporary directory for the ssh control socket."
[[ "$SSH_CONTROL_DIR" != *" "* ]] || die "the temporary directory '$SSH_CONTROL_DIR' contains spaces: set TMPDIR to a path without spaces."
SSH_OPTS=(-o ControlMaster=auto -o "ControlPath=$SSH_CONTROL_DIR/%r@%h:%p" -o ControlPersist=300 -o ServerAliveInterval=60)
cleanup() {
    ssh "${SSH_OPTS[@]}" -O exit "$CLUSTER_HOST" > /dev/null 2>&1 || true
    rm -rf "$SSH_CONTROL_DIR"
}
trap cleanup EXIT
remote() { ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" 'bash -s' <<< "$1"; }
REMOTE_DIR_Q="$(printf '%q' "$CLUSTER_DIR")"

# ---------------------------------------------------------------------------
# Wait for the job, then add its log to the paths to fetch
# ---------------------------------------------------------------------------
if [ -n "$WAIT_JOB" ]; then
    echo "Waiting for job $WAIT_JOB to leave the queue (checked every ${POLL_SECONDS}s; Ctrl-C stops waiting, the job keeps running)..."
    while true; do
        QUEUE_STATE="$(remote "squeue -h -j $WAIT_JOB -o %T 2>/dev/null || true")" \
            || die "ssh failed while checking job $WAIT_JOB."
        [ -z "$QUEUE_STATE" ] && break
        printf '  [%s] job %s: %s\n' "$(date '+%H:%M:%S')" "$WAIT_JOB" "$QUEUE_STATE"
        sleep "$POLL_SECONDS"
    done
    FINAL_STATE="$(remote "sacct -j $WAIT_JOB -X -n -o State 2>/dev/null | head -n 1 | tr -d ' '" || true)"
    echo "Job $WAIT_JOB ended with state: ${FINAL_STATE:-<unknown: sacct returned nothing>}"
    if [ "$FINAL_STATE" != "COMPLETED" ]; then
        echo "WARNING: job $WAIT_JOB did not end COMPLETED: outputs may be missing or partial; read its log first."
    fi
    # Job logs: launcher.sh writes snakemake_<id>.out, submit_from_ws.sh pilot jobs write logs_pilot/<name>_<id>.out
    JOB_LOGS="$(remote "cd $REMOTE_DIR_Q && ls -1 snakemake_${WAIT_JOB}.out logs_pilot/*_${WAIT_JOB}.out 2>/dev/null || true")"
    if [ -n "$JOB_LOGS" ]; then
        while read -r log_path; do REMOTE_PATHS+=("$log_path"); done <<< "$JOB_LOGS"
    else
        echo "WARNING: no log of job $WAIT_JOB found (looked for snakemake_${WAIT_JOB}.out and logs_pilot/*_${WAIT_JOB}.out)."
    fi
fi
[ ${#REMOTE_PATHS[@]} -gt 0 ] || { echo "Nothing to fetch."; exit 0; }

# ---------------------------------------------------------------------------
# Measure the size on the cluster (same exclusions as the copy)
# ---------------------------------------------------------------------------
EXCLUDE_ARGS=()
if [ "$INCLUDE_BINARY" -eq 0 ]; then
    for pattern in "${BINARY_PATTERNS[@]}"; do EXCLUDE_ARGS+=("--exclude=$pattern"); done
fi
# Quoted for the remote shell; empty when nothing is excluded (an empty '' would be read by du as a file name)
EXCLUDE_ARGS_Q=""
[ ${#EXCLUDE_ARGS[@]} -gt 0 ] && EXCLUDE_ARGS_Q="$(printf '%q ' "${EXCLUDE_ARGS[@]}")"
SIZE_OUTPUT="$(remote "cd $REMOTE_DIR_Q || exit 3
missing=0
for p in $(printf '%q ' "${REMOTE_PATHS[@]}"); do
    [ -e \"\$p\" ] || { echo \"MISSING \$p\"; missing=1; }
done
[ \$missing -eq 0 ] || exit 5
du -sb $EXCLUDE_ARGS_Q $(printf '%q ' "${REMOTE_PATHS[@]}")
")" || die "could not measure the requested paths on the cluster:
$(echo "$SIZE_OUTPUT" | sed -n 's/^MISSING /    not found: /p')"

echo "Size on the cluster (binary files $( [ "$INCLUDE_BINARY" -eq 1 ] && echo included || echo excluded )):"
echo "$SIZE_OUTPUT" | awk '{printf "    %10.1f MB  %s\n", $1/1048576, $2}'
TOTAL_BYTES="$(echo "$SIZE_OUTPUT" | awk '{s += $1} END {print s + 0}')"
mkdir -p "$LOCAL_DEST"
FREE_LOCAL_BYTES="$(df --output=avail -B1 "$LOCAL_DEST" | tail -n 1 | tr -d ' ')"
printf 'Total: %.1f MB (limit %s); free on the local disk: %.1f GB\n' \
    "$(echo "$TOTAL_BYTES / 1048576" | bc -l)" "$MAX_SIZE" "$(echo "$FREE_LOCAL_BYTES / 1073741824" | bc -l)"

if [ "$TOTAL_BYTES" -gt "$MAX_BYTES" ]; then
    die "the requested files total more than --max-size $MAX_SIZE. Fetch narrower paths, or raise --max-size if you really need them all."
fi
if [ $(( FREE_LOCAL_BYTES - TOTAL_BYTES )) -lt "$MIN_FREE_LOCAL_BYTES" ]; then
    die "copying would leave less than 20 GB free on the local disk. Free space first or fetch less."
fi

# ---------------------------------------------------------------------------
# Copy (rsync --relative keeps the path below the project folder)
# ---------------------------------------------------------------------------
RSYNC_SOURCES=()
for remote_path in "${REMOTE_PATHS[@]}"; do
    RSYNC_SOURCES+=("$CLUSTER_HOST:$CLUSTER_DIR/./$remote_path")
done
# The ssh options contain no spaces (checked when the control directory is created), so they can be
# passed to rsync -e as one plain string.
RSYNC_ARGS=(-a --relative --human-readable --itemize-changes -e "ssh ${SSH_OPTS[*]}" ${EXCLUDE_ARGS[@]+"${EXCLUDE_ARGS[@]}"})
[ "$DRY" -eq 1 ] && RSYNC_ARGS+=(--dry-run)

rsync "${RSYNC_ARGS[@]}" "${RSYNC_SOURCES[@]}" "$LOCAL_DEST/"

if [ "$DRY" -eq 1 ]; then
    echo "--dry given: nothing copied (the list above is what would be copied)."
else
    echo "Done: files are under $LOCAL_DEST/<same relative path>."
fi
