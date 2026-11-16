#!/usr/bin/env bash
# cluster_status.sh - mirror on the workstation the state of your jobs on the cluster (read-only).
#
# Purpose
#   One view, from the workstation, of:
#     - the cluster checkout commit vs the local HEAD (are you running the code you think?);
#     - squeue: your pending/running jobs (state, elapsed / time limit, CPUs, memory, reason, dependency);
#     - for each running Snakemake job: the last progress line of its log ("N of M steps (P%) done");
#     - sacct: your jobs that ended in the last --since hours (state, exit code, elapsed, CPUs, memory).
#   With --job <id>: details of one job (sacct with the per-step peak memory MaxRSS, and the last
#   lines of its log).
#   It only runs read-only commands on the cluster (git rev-parse, squeue, sacct, grep/tail of logs).
#
# Usage
#   launcher/cluster_status.sh                  # one snapshot
#   launcher/cluster_status.sh --watch [sec]    # refresh every <sec> seconds (default 60), Ctrl-C to stop
#   launcher/cluster_status.sh --since 48       # finished jobs of the last 48 hours (default 24)
#   launcher/cluster_status.sh --job 123456 [--lines 40]
#
# Logs looked up: snakemake_<id>.out (launcher.sh) and logs_pilot/*_<id>.out (pilot jobs of submit_from_ws.sh).
# Environment overrides (testing only): CLUSTER_HOST, CLUSTER_DIR.
# The ssh password is asked once: one shared ssh connection is reused (also across --watch refreshes).

set -euo pipefail
export LC_ALL=C

CLUSTER_HOST="${CLUSTER_HOST:-ieo7429@login.ieo.it}"
CLUSTER_DIR="${CLUSTER_DIR:-/hpcnfs/scratch/MS/ieo7429/SelectionCaller}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

die() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
WATCH_SECONDS=""
SINCE_HOURS=24
JOB_ID=""
LOG_LINES=30
while [ $# -gt 0 ]; do
    case "$1" in
        --watch)
            if [[ "${2:-}" =~ ^[0-9]+$ ]]; then WATCH_SECONDS="$2"; shift 2; else WATCH_SECONDS=60; shift; fi ;;
        --since)   SINCE_HOURS="${2:?--since needs a number of hours}"; shift 2 ;;
        --job)     JOB_ID="${2:?--job needs a job id}"; shift 2 ;;
        --lines)   LOG_LINES="${2:?--lines needs a number}"; shift 2 ;;
        -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
        *)         die "unknown argument '$1' (see --help)." ;;
    esac
done
[[ "$SINCE_HOURS" =~ ^[0-9]+$ ]] || die "--since expects an integer number of hours, got '$SINCE_HOURS'."
[[ "$LOG_LINES" =~ ^[0-9]+$ ]]   || die "--lines expects an integer, got '$LOG_LINES'."
[ -z "$JOB_ID" ] || [[ "$JOB_ID" =~ ^[0-9]+$ ]] || die "--job expects a numeric SLURM job id, got '$JOB_ID'."
[ -z "$WATCH_SECONDS" ] || [ "$WATCH_SECONDS" -ge 10 ] || die "--watch interval must be at least 10 seconds (do not load the login node)."

# ---------------------------------------------------------------------------
# One shared ssh connection (password typed once)
# ---------------------------------------------------------------------------
SSH_CONTROL_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sc_ssh.XXXXXX")"
[ -n "$SSH_CONTROL_DIR" ] && [ -d "$SSH_CONTROL_DIR" ] || die "cannot create a temporary directory for the ssh control socket."
SSH_OPTS=(-o ControlMaster=auto -o "ControlPath=$SSH_CONTROL_DIR/%r@%h:%p" -o ControlPersist=600 -o ServerAliveInterval=60)
cleanup() {
    ssh "${SSH_OPTS[@]}" -O exit "$CLUSTER_HOST" > /dev/null 2>&1 || true
    rm -rf "$SSH_CONTROL_DIR"
}
trap cleanup EXIT
remote() { ssh "${SSH_OPTS[@]}" "$CLUSTER_HOST" 'bash -s' <<< "$1"; }
REMOTE_DIR_Q="$(printf '%q' "$CLUSTER_DIR")"

# ---------------------------------------------------------------------------
# Remote scripts (read-only)
# ---------------------------------------------------------------------------
# Overview: commit, queue, progress of running snakemake jobs, recently finished jobs
OVERVIEW_SCRIPT="cd $REMOTE_DIR_Q || exit 3
echo \"=== Cluster checkout: \$(git rev-parse --short HEAD)\$( [ -n \"\$(git status --porcelain --untracked-files=no)\" ] && echo '  (WARNING: modified tracked files)')\"
echo
echo '=== Queue (squeue) ==='
squeue -u \"\$USER\" -o '%.10i %.45j %.9T %.11M %.11l %.5C %.7m %.22R %E'
echo
echo '=== Progress of running Snakemake jobs (last progress line of the log) ==='
found=0
for job_id in \$(squeue -h -u \"\$USER\" -t RUNNING -o '%i %j' | awk '\$2 ~ /^snakemake/ {print \$1}'); do
    found=1
    log=snakemake_\${job_id}.out
    if [ -f \"\$log\" ]; then
        progress=\$(grep -E '[0-9]+ of [0-9]+ steps' \"\$log\" | tail -n 1)
        errors=\$(grep -c -E '^Error in rule' \"\$log\" || true)
        echo \"  \$job_id: \${progress:-<no progress line yet>}  (rule errors so far: \$errors)\"
    else
        echo \"  \$job_id: log \$log not found\"
    fi
done
[ \$found -eq 1 ] || echo '  (none)'
echo
echo '=== Ended in the last $SINCE_HOURS h (sacct) ==='
sacct -u \"\$USER\" -X -S now-${SINCE_HOURS}hours -s CD,F,TO,CA,OOM,NF \\
      -o JobID%10,JobName%45,State%12,ExitCode%8,Elapsed%11,AllocCPUS%5,ReqMem%8,End%19
"

# One job: accounting with per-step peak memory, then the tail of its log
JOB_SCRIPT="cd $REMOTE_DIR_Q || exit 3
echo '=== Job $JOB_ID (sacct, all steps) ==='
sacct -j $JOB_ID -o JobID%16,JobName%30,State%12,ExitCode%8,Start%19,Elapsed%11,AllocCPUS%5,ReqMem%8,MaxRSS%10 --units=G
echo
log=\$(ls -1 snakemake_${JOB_ID}.out logs_pilot/*_${JOB_ID}.out 2>/dev/null | head -n 1)
if [ -n \"\$log\" ]; then
    echo \"=== Last $LOG_LINES lines of \$log ===\"
    tail -n $LOG_LINES \"\$log\"
    errors=\$(grep -E '^Error in rule' \"\$log\" | sort | uniq -c || true)
    [ -z \"\$errors\" ] || { echo; echo '=== Failed rules in the log ==='; echo \"\$errors\"; }
else
    echo 'No log found (looked for snakemake_${JOB_ID}.out and logs_pilot/*_${JOB_ID}.out).'
fi
"

show_once() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $CLUSTER_HOST:$CLUSTER_DIR   (local HEAD: $(git -C "$REPO_ROOT" rev-parse --short HEAD))"
    if [ -n "$JOB_ID" ]; then
        remote "$JOB_SCRIPT" || die "ssh/remote command failed."
    else
        remote "$OVERVIEW_SCRIPT" || die "ssh/remote command failed."
    fi
}

if [ -z "$WATCH_SECONDS" ]; then
    show_once
else
    while true; do
        OUTPUT="$(show_once 2>&1)" || { echo "$OUTPUT"; exit 1; }
        clear
        echo "$OUTPUT"
        echo
        echo "(refresh every ${WATCH_SECONDS}s, Ctrl-C to stop)"
        sleep "$WATCH_SECONDS"
    done
fi
