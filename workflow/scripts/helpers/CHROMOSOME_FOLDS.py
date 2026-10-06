"""Leave-chromosome(s)-out splitting (BLOCKING == "chromosome").

Every split keeps whole chromosomes together: a chromosome is never in both the training and
the evaluation side. Indices returned are POSITIONAL (0..n-1) with respect to the binIDs passed
in, so they can be used on the matching X / y rows (same convention as the other K-fold helpers
is ASSUMED - check against GENERATE_KFOLD_IDXS / SUBSET_DATA_BY_IDXS before wiring in).

binIDs have the format "{chrom}_{bin_num}" (see resolve_chroms in helpers_ML.py).
"""
import numpy as np
from sklearn.model_selection import GroupKFold

CHROMOSOME_BLOCKING = "chromosome"


def CHROMS_FROM_BINIDS(binIDs):
    """Chromosome label of every bin, as a numpy array of strings."""
    binIDs = np.asarray(binIDs, dtype=str)
    chroms = np.array([b.split("_")[0] for b in binIDs])
    if np.any(chroms == "") or not all("_" in b for b in binIDs):
        raise ValueError("CHROMS_FROM_BINIDS: binIDs must have the format '{chrom}_{bin_num}'.")
    return chroms


def CHROMOSOME_TRAIN_TEST_SPLIT(binIDs, TEST_SIZE, SEED):
    """Held-out test set made of WHOLE chromosomes.

    Chromosomes are visited in a seeded random order and added to the test set until its
    fraction of bins reaches TEST_SIZE; the last addition is kept only if it brings the
    fraction closer to TEST_SIZE. At least one chromosome ends up on each side.
    Returns (train_idx, test_idx), positional and sorted.
    """
    if not 0 < TEST_SIZE < 1:
        raise ValueError(f"CHROMOSOME_TRAIN_TEST_SPLIT: TEST_SIZE must be in (0, 1), got {TEST_SIZE}.")
    chroms = CHROMS_FROM_BINIDS(binIDs)
    uniq, counts = np.unique(chroms, return_counts=True)
    if len(uniq) < 2:
        raise ValueError("CHROMOSOME_TRAIN_TEST_SPLIT: need at least 2 chromosomes, found "
                         f"{len(uniq)} - BLOCKING = 'chromosome' is only valid in genome_wide mode.")
    n = len(chroms)
    order = np.random.default_rng(SEED).permutation(len(uniq))
    test_chroms, frac = [], 0.0
    for i in order[:-1]:                                   # never put ALL chromosomes in test
        new_frac = frac + counts[i] / n
        if frac > 0 and abs(new_frac - TEST_SIZE) > abs(frac - TEST_SIZE):
            break
        test_chroms.append(uniq[i]); frac = new_frac
        if frac >= TEST_SIZE:
            break
    is_test = np.isin(chroms, test_chroms)
    return np.where(~is_test)[0], np.where(is_test)[0]


def CHROMOSOME_KFOLD_IDXS(binIDs, N_SPLITS):
    """K folds of whole chromosomes (GroupKFold, groups = chromosome; fold sizes balanced by
    number of bins). With N_SPLITS == number of chromosomes this is exact leave-one-chromosome-out.
    Returns {fold: (train_idx, val_idx)}, positional and sorted.
    """
    chroms = CHROMS_FROM_BINIDS(binIDs)
    n_chroms = len(np.unique(chroms))
    if n_chroms < N_SPLITS:
        raise ValueError(f"CHROMOSOME_KFOLD_IDXS: {N_SPLITS} folds requested but only {n_chroms} "
                         "chromosome(s) available - each fold needs at least one whole chromosome.")
    splitter = GroupKFold(n_splits=N_SPLITS)
    dummy_X = np.zeros((len(chroms), 1))
    return {fold: (np.sort(tr), np.sort(va))
            for fold, (tr, va) in enumerate(splitter.split(dummy_X, groups=chroms))}


def _TAKE(obj, idx):
    """Positional row subset for DataFrame / Series / array-like (keeps pandas index labels)."""
    return obj.iloc[idx] if hasattr(obj, "iloc") else np.asarray(obj)[idx]


def CHROMOSOME_SPLIT_DATASET(DATASET, TARGET_COL, COLS_TO_EXCLUDE, TEST_SIZE, SEED, SPLIT_DATASET_FUN):
    """Drop-in replacement for SPLIT_DATASET(...) when BLOCKING == "chromosome".

    Reuses SPLIT_DATASET's own X / y / binIDs construction (EXIT_EARLY=True: same target and
    excluded columns, no blocking/sorting), then holds out whole chromosomes as the test set.
    Returns the same 6-tuple as SPLIT_DATASET:
    X_train, X_test, y_train, y_test, binIDs_train, binIDs_test
    """
    X, y, binIDs = SPLIT_DATASET_FUN(
        DATASET,
        TARGET_COL = TARGET_COL,
        COLS_TO_EXCLUDE = COLS_TO_EXCLUDE,
        EXIT_EARLY = True
    )
    train_idx, test_idx = CHROMOSOME_TRAIN_TEST_SPLIT(binIDs, TEST_SIZE, SEED)
    return (_TAKE(X, train_idx), _TAKE(X, test_idx),
            _TAKE(y, train_idx), _TAKE(y, test_idx),
            _TAKE(binIDs, train_idx), _TAKE(binIDs, test_idx))