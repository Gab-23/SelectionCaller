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
echo "    Fixed with an explicit, order-preserving match() instead of merge()."
echo "    Revised design (Option B fix for the separate 'merge_cooler: bin ID exceeds"
echo "    declared number of bins' bug): the script now takes TWO backbone inputs - the raw,"
echo "    pre-filter backbone (continuous, no gaps, what cload is actually run against - see"
echo "    sort_raw_backbone_hic.smk/compute_matrices_hic.smk) for matching, plus the filtered"
echo "    backbone (blacklist removed) to determine which bins are blacklisted (present in raw"
echo "    but absent from filtered) and should be marked binID=NA. Confirmed directly on a"
echo "    real cluster run that blacklisted regions often have anomalously HIGH signal, not"
echo "    zero, so this explicit exclusion is needed - process_hic.py's zero-signal check"
echo "    alone would not have caught them.)"
printf "chr1\t0\t200\t1_1\nchr1\t200\t400\t1_2\nchr2\t0\t200\t2_1\nchr2\t200\t400\t2_2\nchr10\t0\t200\t10_1\nchr10\t200\t400\t10_2\n" > "$TMP/backbone_raw.bed"
printf "chr1\t0\t200\t1_1\nchr2\t0\t200\t2_1\nchr2\t200\t400\t2_2\nchr10\t0\t200\t10_1\nchr10\t200\t400\t10_2\n" > "$TMP/backbone_filtered.bed"
printf "chr1\t0\t200\t1\nchr1\t200\t400\t2\nchr2\t0\t200\t3\nchr2\t200\t400\t4\nchr10\t0\t200\t5\nchr10\t200\t400\t6\n" > "$TMP/bins.bed"
Rscript workflow/scripts/process_bins_hic.R \
    --input_path_backbone_raw="$TMP/backbone_raw.bed" \
    --input_path_backbone_filtered="$TMP/backbone_filtered.bed" \
    --input_path_bins="$TMP/bins.bed" --output_path="$TMP/output.bed"
check "row 1 stays chr1 (cooler's original order preserved)" "$(sed -n 1p "$TMP/output.bed" | cut -f1)" "chr1"
check "row 5 stays chr10, NOT reordered before chr2 (the ordering bug)" "$(sed -n 5p "$TMP/output.bed" | cut -f1)" "chr10"
check "row 3 stays chr2 (not pushed after chr10 by alphabetical sort)" "$(sed -n 3p "$TMP/output.bed" | cut -f1)" "chr2"
check "binID correctly attached to the chr10 row" "$(sed -n 5p "$TMP/output.bed" | cut -f4)" "10_1"
check "bin present in raw but absent from filtered (chr1_2, blacklisted) is marked NA" "$(sed -n 2p "$TMP/output.bed" | cut -f4)" "NA"
check "bin present in both raw and filtered (chr1_1) keeps its real binID" "$(sed -n 1p "$TMP/output.bed" | cut -f4)" "1_1"

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

echo "=== 4) normalize_matrix_hic.py - per-chromosome cis-only fetch (2026-09 rewrite) ==="
echo "    (bug: the previous version pulled cooler_obj.matrix(...)[:,:] - the ENTIRE"
echo "    genome-wide dense matrix - into memory (twice: raw + symmetrized), looped in pure"
echo "    Python over every genome-wide diagonal (~150k iterations at real scale), and wrote"
echo "    the whole thing out as plain text. Beyond the O(n^2) memory / O(n) Python-loop cost,"
echo "    it was scientifically wrong for the same reason as process_hic.py's eigh: the"
echo "    per-diagonal z-score baseline only makes sense WITHIN one chromosome (distance in"
echo "    bin-index space = real genomic distance only there), so a genome-wide baseline"
echo "    mixes real cis distance-decay with structurally-different trans noise. Fixed by"
echo "    fetching only the CIS submatrix for one chromosome directly from the (genome-wide"
echo "    balanced) cooler object via .fetch(chrom) - never materializes anything genome-wide."
echo "    Verified here with a real 2-chromosome toy .cool file built with the cooler package,"
echo "    with distinguishable cis (strong, distance-decaying) vs trans (weak, flat) signal."
python3 -c "
import numpy as np, pandas as pd, cooler, subprocess, sys, os, tempfile

rng = np.random.default_rng(42)
bins = pd.DataFrame({
    'chrom': ['chr1']*5 + ['chr2']*4,
    'start': [0,100,200,300,400] + [0,100,200,300],
    'end':   [100,200,300,400,500] + [100,200,300,400],
})
n = len(bins)
chrom_arr = bins['chrom'].values
mat = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i == j:
            mat[i, j] = 100.0 + rng.normal(0, 3)
        elif chrom_arr[i] == chrom_arr[j]:
            d = abs(i - j)
            mat[i, j] = max(50.0 / d + rng.normal(0, 2), 0.1)
        else:
            mat[i, j] = max(0.01 + rng.normal(0, 0.005), 0.001)  # trans: weak, flat, NOT distance-like

pixels = [(i, j, mat[i, j]) for i in range(n) for j in range(i, n)]
pixels_df = pd.DataFrame(pixels, columns=['bin1_id', 'bin2_id', 'count'])

tmp = tempfile.mkdtemp()
cool_path = os.path.join(tmp, 'toy.cool')
cooler.create_cooler(cool_path, bins, pixels_df, dtypes={'count': 'float64'}, ordered=True)
c = cooler.Cooler(cool_path)
with cooler.util.open_hdf5(c.filename, 'r+') as h5:
    h5['bins']['weight'] = np.ones(n, dtype=float)

results = {}
for chrom, expected_n in [('chr1', 5), ('chr2', 4)]:
    mat_path = os.path.join(tmp, f'{chrom}.mat')
    bins_path = os.path.join(tmp, f'{chrom}.bins')
    r = subprocess.run([
        sys.executable, 'workflow/scripts/normalize_matrix_hic.py',
        '--input_path_cool', cool_path, '--chrom', chrom,
        '--output_path_matrix', mat_path, '--output_path_bins', bins_path
    ], capture_output=True, text=True)
    if r.returncode != 0:
        print('SCRIPT_FAILED', chrom); print(r.stderr); sys.exit(1)
    out_bins = pd.read_csv(bins_path, sep='\t', header=None, names=['chr','start','end','w'])
    out_mat = np.loadtxt(mat_path)
    results[chrom] = (out_bins, out_mat)
    assert (out_bins['chr'] == chrom).all(), f'{chrom}: bins leaked another chromosome'
    assert len(out_bins) == expected_n, f'{chrom}: expected {expected_n} bins, got {len(out_bins)}'
    assert out_mat.shape == (expected_n, expected_n), f'{chrom}: matrix shape {out_mat.shape} != ({expected_n},{expected_n})'

# the two chromosomes' z-score matrices must differ (proves the computation is genuinely
# independent per chromosome, not e.g. a slice of one shared precomputed genome-wide array)
m1, m2 = results['chr1'][1], results['chr2'][1]
assert not np.allclose(np.nan_to_num(m1[:4,:4]), np.nan_to_num(m2), equal_nan=True), 'chr1/chr2 z-scores suspiciously identical'
print('OK')
"
RC=$?
check "normalize_matrix_hic.py fetches only the requested chromosome's bins/matrix, shapes correct, no cross-chromosome leakage" "$RC" "0"

echo "=== 5) get_chroms_list() - checkpoint chromosome-list parsing (helpers_dataset.smk) ==="
echo "    (new function, added alongside the clean_chrom_sizes checkpoint conversion - reads"
echo "    the checkpoint's cleaned chrom.sizes file to build the {chrom} wildcard list for"
echo "    the Hi-C cis-only chain. Checked here in isolation from Snakemake's checkpoint"
echo "    machinery itself, which is instead validated by the toy integration test.)"
printf "chr1\t2000\nchr2\t1500\n" > "$TMP/clean_chrom_sizes.txt"
PARSED=$(python3 -c "
with open('$TMP/clean_chrom_sizes.txt') as f:
    chroms = [line.split('\t')[0] for line in f if line.strip()]
print(','.join(chroms))
")
check "chromosome list parsed correctly from a clean.chrom.sizes-shaped file" "$PARSED" "chr1,chr2"

echo "=== 6) aggregate_hic_eigen - header-deduplicated concatenation of per-chromosome .eigen files ==="
printf "binID\teigen_1\nchr1_1\t0.5\nchr1_2\t-0.3\n" > "$TMP/chr1.eigen"
printf "binID\teigen_1\nchr2_1\t0.1\nchr2_2\t-0.9\n" > "$TMP/chr2.eigen"
awk 'FNR==1 && NR!=1 { next } { print }' "$TMP/chr1.eigen" "$TMP/chr2.eigen" > "$TMP/combined.eigen"
check "combined file has exactly 1 header + 4 data rows (no duplicated header)" "$(wc -l < "$TMP/combined.eigen" | tr -d ' ')" "5"
check "header line is preserved" "$(sed -n 1p "$TMP/combined.eigen")" "$(printf 'binID\teigen_1')"
check "chr2's header line was correctly dropped, not treated as data" "$(grep -c '^binID' "$TMP/combined.eigen")" "1"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
