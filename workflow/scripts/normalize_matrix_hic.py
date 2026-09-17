import sys
import cooler
import warnings
import numpy as np
import pandas as pd
import argparse

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_cool', dest='input_path_cool', required=True, type=str, help='')
parser.add_argument('--chrom', dest='chrom', required=True, type=str, help='')
parser.add_argument('--output_path_matrix', dest='output_path_matrix', required=True, type=str, help='')
parser.add_argument('--output_path_bins', dest='output_path_bins', required=True, type=str, help='')

args = parser.parse_args()

COOL_PATH = args.input_path_cool
CHROM = args.chrom
OUTFILE_MATRIX_ZSCORE = args.output_path_matrix
OUTFILE_BINS = args.output_path_bins

print(f"loading Cooler object, fetching cis-only region for {CHROM}")
cooler_obj = cooler.Cooler(COOL_PATH)

# NOTE (fix, see conversation): the previous version pulled cooler_obj.matrix(...)[:,:] -
# the ENTIRE genome-wide dense matrix - into memory (twice: raw + symmetrized), looped in
# pure Python over every one of its ~150k diagonals, and wrote the whole thing out as a
# plain-text file. Two separate problems: (1) computationally infeasible at genome scale
# (hundreds of GB, per-diagonal Python loop over ~150k iterations); (2) scientifically wrong
# for the SAME reason as process_hic.py's eigh - "distance between bin i and bin j" in
# diag-index space only corresponds to real genomic distance WITHIN one chromosome, so a
# genome-wide z-score baseline mixes cis decay (real, distance-dependent) with trans noise
# (structurally different) into the same per-diagonal mean/std. Fetching the cis-only
# submatrix for one chromosome via cooler's own indexed on-disk fetch (never materializing
# anything genome-wide) fixes both at once - the .cool file was already balanced genome-wide
# upstream (balance_cooler.smk), which is the one part of ICE correction that legitimately
# benefits from genome-wide fitting of technical bias factors; this step only takes the
# CIS SLICE of that already-balanced matrix.
cooler_bins = cooler_obj.bins().fetch(CHROM)[:].reset_index(drop=True)
count_matrix_iced = cooler_obj.matrix(balance=True).fetch(CHROM)
count_matrix_iced = count_matrix_iced.astype(float)
count_matrix_iced_symm = np.maximum(count_matrix_iced, count_matrix_iced.T)  # symmetrize

test = count_matrix_iced[10:25, 10:25]
test_symm = count_matrix_iced_symm[10:25, 10:25]

print(pd.DataFrame(test))
print(np.allclose(test, test.T, equal_nan=True))
print(pd.DataFrame(test_symm))
print(np.allclose(test_symm, test_symm.T, equal_nan=True))

print(f"computing diagonal mean and std for {CHROM}")
n_rows = count_matrix_iced.shape[0]

diag_mean_arr = np.empty(n_rows, dtype=float)
diag_std_arr = np.empty(n_rows, dtype=float)

for K in range(n_rows):  # now O(n_chrom_bins), not O(n_genome_bins) - a single chromosome's
                          # bin count (thousands, not ~150k) makes this loop cheap again
    diagonal = np.diag(count_matrix_iced_symm, K)
    diag_mean_arr[K] = np.nanmean(diagonal)
    diag_std_arr[K] = np.nanstd(diagonal)

idx = np.arange(n_rows)
distance_matrix = np.abs(idx[:, None] - idx[None, :])

print("computing z-scores")
count_matrix_zscore = (count_matrix_iced_symm - diag_mean_arr[distance_matrix]) / diag_std_arr[distance_matrix]

test_zscore = count_matrix_zscore[10:25, 10:25]
print(pd.DataFrame(test_zscore))
print(np.allclose(test_zscore, test_zscore.T, equal_nan=True))

print("saving files")

fmt_matrix = "%.4f"
fmt_bins = ["%s", "%.0f", "%.0f", "%.6f"]

np.savetxt(OUTFILE_MATRIX_ZSCORE, count_matrix_zscore, fmt=fmt_matrix, delimiter="\t")
np.savetxt(OUTFILE_BINS, cooler_bins, fmt=fmt_bins, delimiter="\t")

print("Done!")
exit(0)