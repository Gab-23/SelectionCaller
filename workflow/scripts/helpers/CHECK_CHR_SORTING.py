import numpy as np

def CHECK_CHR_SORTING(binIDs):

    chromosomes = [int(elem.split("_")[0]) for elem in binIDs.tolist()]
    bin_nums = [int(elem.split("_")[1]) for elem in binIDs.tolist()]

    # NOTE (fix): the original computed, per unique chromosome, a "chromosomes_f" sublist that
    # was ALWAYS composed of the same repeated value (since it was extracted via a mask that
    # already filters for that one chromosome) - so "chromosomes_f == sorted(chromosomes_f)" was
    # trivially True regardless of whether different chromosomes were contiguous or interleaved
    # in the original array. Confirmed with a direct test: an array with chr1/chr2 fully
    # interleaved was incorrectly reported as sorted. This check gates GENERATE_KFOLD_IDXS's
    # naive contiguous-row-index KFold split, which silently produces spatially meaningless
    # folds if chromosomes aren't genuinely grouped in contiguous, increasing blocks - exactly
    # the failure this check exists to catch. Per your confirmation this case never occurs with
    # your actual (already correctly sorted) data, but is worth guarding against regardless.
    if chromosomes != sorted(chromosomes):
        return False

    unique_chromosomes = sorted(set(chromosomes))

    for chr in unique_chromosomes:
        idxs = [idx for idx, elem in enumerate(chromosomes) if elem == chr]
        bin_nums_f = [bin_nums[idx] for idx in idxs]
        if bin_nums_f != sorted(bin_nums_f):
            return False

    return True
