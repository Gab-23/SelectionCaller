#!/usr/bin/env bash
# clean_fetched.sh - remove LOCAL copies of cluster outputs fetched with launcher/fetch_from_cluster.sh.
#
# Purpose
#   Files fetched into ~/Desktop/SelectionCallerCodeBase/_from_cluster/ are temporary working copies
#   (the workstation disk is almost full). This script deletes them once analysed.
#   It NEVER touches the cluster: it makes no ssh connection, the originals stay permanent in
#   /hpcnfs/scratch/MS/ieo7429/SelectionCaller and can be fetched again at any time.
#
# Usage
#   launcher/clean_fetched.sh [--yes] [<path relative to _from_cluster/> ...]
#     no path   remove everything under _from_cluster/
#     <path>    remove only these (same relative paths as on the cluster, e.g. analyses/v4/LUAD/1000000bp)
#     --yes     do not ask for confirmation
#
# Before deleting it prints what will be removed (number of files, size) and asks for confirmation.

set -euo pipefail
export LC_ALL=C

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOCAL_DEST="$REPO_ROOT/_from_cluster"

die() { echo "ERROR: $*" >&2; exit 1; }

ASSUME_YES=0
TARGET_PATHS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --yes)     ASSUME_YES=1; shift ;;
        -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
        --*)       die "unknown option '$1' (valid: --yes)." ;;
        *)         TARGET_PATHS+=("$1"); shift ;;
    esac
done

[ -d "$LOCAL_DEST" ] || { echo "Nothing to clean: $LOCAL_DEST does not exist."; exit 0; }

# Resolve every target and make sure it lies inside _from_cluster/ (no absolute paths, no '..', no symlink escapes)
TO_REMOVE=()
if [ ${#TARGET_PATHS[@]} -eq 0 ]; then
    while IFS= read -r -d '' entry; do TO_REMOVE+=("$entry"); done \
        < <(find "$LOCAL_DEST" -mindepth 1 -maxdepth 1 -print0)
else
    for target in "${TARGET_PATHS[@]}"; do
        case "$target" in
            /*)                  die "'$target' is absolute: give paths relative to _from_cluster/." ;;
            ..|../*|*/..|*/../*) die "'$target' contains '..': only paths inside _from_cluster/ are allowed." ;;
        esac
        full_path="$LOCAL_DEST/$target"
        [ -e "$full_path" ] || die "'$target' not found under $LOCAL_DEST (nothing removed)."
        resolved="$(realpath "$full_path")"
        [[ "$resolved" == "$(realpath "$LOCAL_DEST")"/* ]] || die "'$target' resolves outside _from_cluster/ ($resolved): refused."
        TO_REMOVE+=("$full_path")
    done
fi
[ ${#TO_REMOVE[@]} -gt 0 ] || { echo "Nothing to clean: $LOCAL_DEST is empty."; exit 0; }

# Preview
N_FILES="$(find "${TO_REMOVE[@]}" -type f | wc -l)"
TOTAL_SIZE="$(du -sch "${TO_REMOVE[@]}" | tail -n 1 | cut -f1)"
echo "Local copies to remove (the originals on the cluster are NOT touched):"
for entry in "${TO_REMOVE[@]}"; do
    printf '    %8s  %s\n' "$(du -sh "$entry" | cut -f1)" "${entry#"$LOCAL_DEST"/}"
done
echo "Total: $N_FILES files, $TOTAL_SIZE"

if [ "$ASSUME_YES" -eq 0 ]; then
    read -r -p "Remove them? [y/N] " answer
    [[ "$answer" =~ ^[yY]$ ]] || { echo "Nothing removed."; exit 0; }
fi

rm -rf -- "${TO_REMOVE[@]}"
# Remove folders left empty by a partial clean (only inside _from_cluster/)
find "$LOCAL_DEST" -mindepth 1 -type d -empty -delete
echo "Removed $N_FILES files ($TOTAL_SIZE). Fetch them again with launcher/fetch_from_cluster.sh if needed."
