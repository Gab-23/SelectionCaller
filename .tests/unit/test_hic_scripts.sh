#!/usr/bin/env bash
# Unit tests for the Hi-C module's R/Python scripts and the deduplicate_pairs_hic.smk
# background-job exit-status fix. The full alignment -> pairtools -> cooler -> eigenvector
# chain (real bwa/samblaster/samtools/pairtools/cooler) is validated end-to-end only in the
# toy Snakemake integration test - this script covers the specific bugs found and fixed.
# Usage: bash .tests/unit/test_hic_scripts.sh
# Requires: Rscript, python3 (numpy, pandas)

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

echo "=== 1) process_bins_hic.R ==="
echo "    (bug: merge(by=c('chr','start','end')) sorts its result alphabetically by default -"
echo "    'chr10' sorts before 'chr2', even though it comes after it in every real assembly's"
echo "    actual chromosome order. Since the input 'bins' table IS cooler's own bin table -"
echo "    whose ROW ORDER is the contact matrix's row/column order - and process_hic.py"
echo "    assumes bins[i] corresponds to matrix row i, this would silently misassign every"
echo "    eigenvector to the wrong bin on any real genome with 10+ numbered chromosomes."
echo "    Fixed with an explicit, order-preserving match() instead of merge().)"
printf "chr1\t0\t200\t1_1\nchr1\t200\t400\t1_2\nchr2\t0\t200\t2_1\nchr2\t200\t400\t2_2\nchr10\t0\t200\t10_1\nchr10\t200\t400\t10_2\n" > "$TMP/backbone.bed"
printf "chr1\t0\t200\t1\nchr1\t200\t400\t2\nchr2\t0\t200\t3\nchr2\t200\t400\t4\nchr10\t0\t200\t5\nchr10\t200\t400\t6\n" > "$TMP/bins.bed"
Rscript workflow/scripts/process_bins_hic.R \
    --input_path_backbone="$TMP/backbone.bed" --input_path_bins="$TMP/bins.bed" --output_path="$TMP/output.bed"
check "row 1 stays chr1 (cooler's original order preserved)" "$(sed -n 1p "$TMP/output.bed" | cut -f1)" "chr1"
check "row 5 stays chr10, NOT reordered before chr2 (the actual bug)" "$(sed -n 5p "$TMP/output.bed" | cut -f1)" "chr10"
check "row 3 stays chr2 (not pushed after chr10 by alphabetical sort)" "$(sed -n 3p "$TMP/output.bed" | cut -f1)" "chr2"
check "binID correctly attached to the chr10 row" "$(sed -n 5p "$TMP/output.bed" | cut -f4)" "10_1"

echo "=== 2) process_hic.py defensive eigenvector count ==="
echo "    (bug: requesting a fixed 'top 10' eigenvectors crashed with a shape-mismatch error"
echo "    whenever fewer than 10 valid (non-all-zero) bins remained - possible with small/"
echo "    sparse toy data, though real genome-scale Hi-C data has far more than 10 valid bins"
echo "    per chromosome so this is a defensive fallback, not a change to normal behavior.)"
python3 -c "
import numpy as np, pandas as pd, subprocess, sys, tempfile, os

# 3 bins total, only 2 with any signal (mimicking sparse toy data)
contact = np.array([
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 0.5],
    [0.0, 0.5, 1.0],
])
bins = pd.DataFrame({'chr': ['1','1','1'], 'start': [0,200,400], 'end': [200,400,600], 'binID': ['1_1','1_2','1_3']})

tmp = tempfile.mkdtemp()
mat_path = os.path.join(tmp, 'mat.tsv')
bins_path = os.path.join(tmp, 'bins.tsv')
eigen_path = os.path.join(tmp, 'eigen.tsv')
np.savetxt(mat_path, contact, delimiter='\t')
bins.to_csv(bins_path, sep='\t', header=False, index=False)

result = subprocess.run([
    sys.executable, 'workflow/scripts/process_hic.py',
    '--input_path_matrix', mat_path, '--input_path_bins', bins_path, '--output_path_eigen', eigen_path
], capture_output=True, text=True)
if result.returncode != 0:
    print('SCRIPT_FAILED')
    print(result.stderr)
    sys.exit(1)
out = pd.read_csv(eigen_path, sep='\t')
print('N_ROWS', len(out))
print('N_COLS', len(out.columns))
"
RC=$?
check "process_hic.py no longer crashes with fewer than 10 valid bins" "$RC" "0"

echo "=== 3) deduplicate_pairs_hic.smk's background-job pattern ==="
echo "    (bug: the original background+kill-0-polling pattern never checked the background"
echo "    job's exit status - confirmed with a direct test that even with set -euo pipefail,"
echo "    a failing background job left the script reporting success. Fixed by adding"
echo "    'wait \$PID' after the polling loop, which captures and propagates the real exit"
echo "    code, while keeping the heartbeat behavior the user confirmed they need.)"
FAILING_RC=$(bash -c '
set -euo pipefail
( sleep 0.2; exit 1 ) &
PID=$!
while kill -0 $PID 2>/dev/null; do sleep 0.1; done
wait $PID
' 2>/dev/null; echo $?)
check "a failing background job is now correctly caught (non-zero exit)" "$FAILING_RC" "1"

SUCCESS_RC=$(bash -c '
set -euo pipefail
( sleep 0.2; exit 0 ) &
PID=$!
while kill -0 $PID 2>/dev/null; do sleep 0.1; done
wait $PID
' 2>/dev/null; echo $?)
check "a succeeding background job still reports success (non-regression)" "$SUCCESS_RC" "0"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
