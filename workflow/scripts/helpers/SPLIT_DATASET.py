import numpy as np
import sklearn as sk

def SPLIT_DATASET(DATASET, TARGET_COL, COLS_TO_EXCLUDE = None, EXIT_EARLY = False, BLOCKING = None, SEED = None, TEST_SIZE = None, CHECK_SORTING_FUN = None, BLOCK_DATASET_FUN = None):

    if COLS_TO_EXCLUDE is None:
        COLS_TO_EXCLUDE = []

    cols_to_drop = ["binID", TARGET_COL] + list(COLS_TO_EXCLUDE)  # binID is always excluded from
                                                                    # X (it's a row identifier, not
                                                                    # a feature) - not something
                                                                    # that varies by config, unlike
                                                                    # TARGET_COL/COLS_TO_EXCLUDE

    X = DATASET.loc[:, ~DATASET.columns.isin(cols_to_drop)] # not DATASET.columns is in cols_to_drop
    y = DATASET.loc[:, TARGET_COL]
    binIDs = DATASET["binID"]

    if EXIT_EARLY:
        return([X, y, binIDs])

    if BLOCKING == "noblocking":

        X_train, X_test, \
        y_train, y_test, \
        binIDs_train, binIDs_test = sk.model_selection.train_test_split(X, y, binIDs, test_size = TEST_SIZE, random_state = SEED)

    else:
        if not CHECK_SORTING_FUN(binIDs):
            raise ValueError(f"BLOCKING != 'noblocking' and {CHECK_SORTING_FUN} = False")

        X_train, X_test, \
        y_train, y_test, \
        binIDs_train, binIDs_test = BLOCK_DATASET_FUN(X, y, binIDs, int(BLOCKING))

    return([X_train, X_test, y_train, y_test, binIDs_train, binIDs_test])