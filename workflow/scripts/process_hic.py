import sys
import argparse
import warnings
import numpy
import pandas
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_matrix', dest='input_path_matrix', required=True, type=str)
parser.add_argument('--input_path_bins',   dest='input_path_bins',   required=True, type=str)
parser.add_argument('--output_path_eigen', dest='output_path_eigen', required=True, type=str)
#parser.add_argument('--output_path_insulation', dest='output_path_insulation', required=True, type=str, help='')

args = parser.parse_args()

CONTACT_MATRIX = args.input_path_matrix
BINS_FILE      = args.input_path_bins
OUTFILE_EIGEN  = args.output_path_eigen
#OUTFILE_INSULATION = args.output_path_insulation

# Load bins (chr, start, end, binID) — cooler-indexed
bins = pandas.read_csv(BINS_FILE, sep="\t", header=None, names=["chr", "start", "end", "binID"])

contact = numpy.loadtxt(CONTACT_MATRIX, delimiter="\t") # load contact matrix
print("NaNs/infs before cleaning:", numpy.isnan(contact).sum(), numpy.isinf(contact).sum())
contact = numpy.nan_to_num(contact, nan=0.0, posinf=0.0, neginf=0.0) # set NaNs to 0
contact = (contact + contact.T) / 2 # force matrix to be symmetric
print("NaNs/infs after cleaning:", numpy.isnan(contact).sum(), numpy.isinf(contact).sum())

# Sanity check
assert contact.shape[0] == len(bins), \
    f"Matrix rows ({contact.shape[0]}) != bins ({len(bins)})"

row_sums = contact.sum(axis=1) # compute rowSums
bad_bins = row_sums == 0 # get all-0 columns
valid = ~bad_bins # get valid bins
print(f"Bad bins: {bad_bins.sum()} / {len(bad_bins)}")
contact_valid = contact[numpy.ix_(valid, valid)] # subset matrix
eigenvalues, eigenvectors = numpy.linalg.eigh(contact_valid)  # ascending, so reverse
n_valid = contact_valid.shape[0]
n_eigen = min(10, n_valid)  # fix: with fewer than 10 valid bins (possible for a small toy
                            # chromosome, or if an unusually large fraction of bins are
                            # dropped as "bad"), eigh() can only return as many
                            # eigenvectors as there are valid bins - requesting a fixed 10
                            # unconditionally caused a shape-mismatch crash when fewer were
                            # available. Real genome-scale Hi-C data has far more than 10
                            # valid bins per chromosome, so this only ever engages as a
                            # defensive fallback, not a change to normal behavior.
top10 = eigenvectors[:, ::-1][:, :n_eigen] # get top eigenvectors (up to 10)
n_bins = contact.shape[0] # get number of bins
full_eigvec = numpy.full((n_bins, n_eigen), numpy.nan) # initialize with NaNs
full_eigvec[valid, :] = top10 # fill in eigenvectors

# Assemble output: binID + eigenvectors
eigen_cols = [f"eigen_{i+1}" for i in range(n_eigen)]
df_eigen = pandas.DataFrame(full_eigvec, columns=eigen_cols)
df_out   = pandas.concat([bins[["binID"]].reset_index(drop=True), df_eigen], axis=1)
df_out.to_csv(OUTFILE_EIGEN, sep="\t", index=False, na_rep="NA")

