import numpy as np
import sklearn.model_selection as sk

def SPLIT_DATASET(DATASET, BLOCKING, SEED, TEST_SIZE, CHECK_SORTING_FUN):

    X = DATASET.loc[:, (DATASET.columns != "logFC") & (DATASET.columns != "binID")]
    y = DATASET.loc[:, "logFC"]
    binIDs = DATASET["binID"]

    if BLOCKING == "no_blocking":

        X_train, X_test, \
        y_train, y_test, \
        binIDs_train, binIDs_test = sk.model_selection.train_test_split(X, y, binIDs, test_size = TEST_SIZE, random_state = SEED)

    else:
        if not CHECK_SORTING_FUN(binIDs):
            raise ValueError(f"BLOCKING != 'no_blocking' and {CHECK_SORTING_FUN} = False")

        X_train, X_test, \
        y_train, y_test, \
        binIDs_train, binIDs_test = BLOCK_DATASET(X, y, binIDs, int(BLOCKING))

    return([X_train, X_test, y_train, y_test, binIDs_train, binIDs_test])
