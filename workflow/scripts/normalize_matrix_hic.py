import sys
import cooler
import warnings
import numpy as np
import pandas as pd
import argparse

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_cool', dest='input_path_cool', required=True, type=str, help='')
parser.add_argument('--output_path_matrix', dest='output_path_matrix', required=True, type=str, help='')
parser.add_argument('--output_path_bins', dest='output_path_bins', required=True, type=str, help='')

args = parser.parse_args()

COOL_PATH = args.input_path_cool
OUTFILE_MATRIX_ZSCORE = args.output_path_matrix
OUTFILE_BINS = args.output_path_bins

print("loading Cooler object")
cooler_obj = cooler.Cooler(COOL_PATH)

print("extracting bins and matrix")
cooler_bins = cooler_obj.bins()[:] # extract bins
cooler_matrix = cooler_obj.matrix(balance = True, sparse = False) # take the balanced (balance = True) matrix
count_matrix_iced = cooler_matrix[:,:]; count_matrix_iced = count_matrix_iced.astype(float)
count_matrix_iced_symm = np.maximum(count_matrix_iced, count_matrix_iced.T) # symmetrize the matrix

test = count_matrix_iced[10:25,10:25]
test_symm = count_matrix_iced_symm[10:25,10:25]

print(pd.DataFrame(test))
print(np.allclose(test, test.T, equal_nan = True))
print(pd.DataFrame(test_symm))
print(np.allclose(test_symm, test_symm.T, equal_nan = True))

print("computing diagonal mean and std")
n_rows = count_matrix_iced.shape[0]

diag_mean_arr = np.empty(n_rows, dtype = float); diag_std_arr = np.empty(n_rows, dtype = float)

for K in range(n_rows): # iterate over diagonals
    
    diagonal = np.diag(count_matrix_iced_symm, K)
    diag_mean = np.nanmean(diagonal) # take the mean
    diag_std = np.nanstd(diagonal) # take the sd

    diag_mean_arr[K] = diag_mean
    diag_std_arr[K] = diag_std

idx = np.arange(n_rows)
distance_matrix = np.abs(idx[:,None] - idx[None,:])

# the idea is that idx[:,None] makes a row vector with idx, idx[None,:] makes a column vector with idx and the subtraction gives
# a symmetrical matrix where each cell represents the distance between two bins

print("computing z-scores")
count_matrix_zscore = (count_matrix_iced_symm - diag_mean_arr[distance_matrix]) / diag_std_arr[distance_matrix] # compute Z-score normalization

test_zscore = count_matrix_zscore[10:25,10:25]
print(pd.DataFrame(test_zscore))
print(np.allclose(test_zscore, test_zscore.T, equal_nan = True))

print("saving files")

fmt_matrix = "%.4f"
fmt_bins = ["%s","%.0f","%.0f","%.6f"]

np.savetxt(OUTFILE_MATRIX_ZSCORE, count_matrix_zscore, fmt=fmt_matrix, delimiter="\t") # save the output
np.savetxt(OUTFILE_BINS, cooler_bins, fmt=fmt_bins, delimiter="\t") # save the output

print("Done!")
exit(0)