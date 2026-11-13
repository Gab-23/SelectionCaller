"""Interleaved blocked K-fold with a buffer (spatial CV for autocorrelated genomic bins).

Why: contiguous KFold(shuffle=False) on a coordinate-sorted chromosome makes each fold a
compact chunk (often a whole arm). For targets dominated by arm-level levels (CNAs) this asks
the model to extrapolate to an arm it never saw, and hyperparameter search then favours
near-constant predictors. Here every fold samples the whole chromosome instead:

  * bins are grouped in blocks of BLOCK_BINS consecutive bins (per chromosome);
  * blocks are assigned to folds in rotation (block b -> fold b % N_SPLITS), so the test
    blocks of each fold are spread along the chromosome;
  * training excludes every bin within BUFFER_BINS of a test bin of the same chromosome
    (the buffer), so the model never trains on the immediate neighbours of what it predicts.

Distances are computed from the bin numbers in the binIDs ("{chrom}_{bin_num}"), so the
function does not rely on row order. Returns {fold: (train_idx, test_idx)} with POSITIONAL,
sorted indices - same format as GENERATE_KFOLD_IDXS / CHROMOSOME_KFOLD_IDXS.
"""
import numpy as np


def PARSE_BINIDS(binIDs):
    binIDs = np.asarray(binIDs, dtype=str)
    parts = np.char.partition(binIDs, "_")
    chroms, sep, nums = parts[:, 0], parts[:, 1], parts[:, 2]
    if np.any(sep != "_") or np.any(chroms == ""):
        raise ValueError("PARSE_BINIDS: binIDs must have the format '{chrom}_{bin_num}'.")
    try:
        nums = nums.astype(int)
    except ValueError:
        raise ValueError("PARSE_BINIDS: the part after '_' in binIDs must be an integer bin number.")
    return chroms, nums


def INTERLEAVED_KFOLD_IDXS(binIDs, N_SPLITS, BLOCK_BINS, BUFFER_BINS, MIN_TRAIN_FRACTION=0.2):
    if N_SPLITS < 2:
        raise ValueError("INTERLEAVED_KFOLD_IDXS: N_SPLITS must be >= 2.")
    if BLOCK_BINS < 1 or BUFFER_BINS < 0:
        raise ValueError("INTERLEAVED_KFOLD_IDXS: BLOCK_BINS must be >= 1 and BUFFER_BINS >= 0.")

    chroms, nums = PARSE_BINIDS(binIDs)
    n = len(nums)
    if len(set(zip(chroms, nums))) != n:
        raise ValueError("INTERLEAVED_KFOLD_IDXS: duplicated binIDs.")

    fold_of = ((nums - 1) // BLOCK_BINS) % N_SPLITS            # rotation of blocks over folds
    idxs = {}
    for k in range(N_SPLITS):
        is_test = fold_of == k
        dist = np.full(n, np.inf)
        for c in np.unique(chroms):
            on_c = np.where(chroms == c)[0]
            t = np.sort(nums[on_c][is_test[on_c]])
            if t.size == 0:
                continue
            x = nums[on_c]
            j = np.searchsorted(t, x)
            left = np.abs(x - t[np.clip(j - 1, 0, t.size - 1)])
            right = np.abs(t[np.clip(j, 0, t.size - 1)] - x)
            dist[on_c] = np.minimum(left, right)
        test_idx = np.where(is_test)[0]
        train_idx = np.where(~is_test & (dist > BUFFER_BINS))[0]
        if test_idx.size == 0:
            raise ValueError(f"INTERLEAVED_KFOLD_IDXS: fold {k} has no test bins (too few bins for "
                             f"N_SPLITS={N_SPLITS} x BLOCK_BINS={BLOCK_BINS}).")
        if train_idx.size < MIN_TRAIN_FRACTION * n:
            raise ValueError(f"INTERLEAVED_KFOLD_IDXS: fold {k} keeps only {train_idx.size}/{n} training bins "
                             f"with BUFFER_BINS={BUFFER_BINS}; reduce the buffer or increase BLOCK_BINS.")
        idxs[str(k)] = (train_idx, test_idx)
    return idxs