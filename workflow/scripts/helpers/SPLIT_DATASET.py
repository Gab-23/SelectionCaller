def SPLIT_DATASET(DATASET, BLOCKING, SEED, TEST_SIZE):

    X = DATASET.loc[:, DATASET.columns != "logFC"]
    X = DATASET.loc[:, DATASET.columns != "binID"]
    y = DATASET.loc[:, "logFC"]
    binIDs = DATASET[["binID"]]

    if BLOCKING = "no_blocking":

        X_train, X_test, 
        y_train, y_test, 
        binIDs_train, binIDs_test = sk.model_selection.train_test_split(X, y, bin_IDs, test_size = TEST_SIZE, random_state = SEED)

    else:
        if CHECK_CHR_SORTING(X, binIDs):
            raise ValueError("BLOCKING != 'no_blocking' and CHECK_CHR_SORTING = False")

        X_train, X_test, 
        y_train, y_test, 
        binIDs_train, binIDs_test = BLOCK_DATASET(X, y, binIDs, BLOCKING)

    return([X_train, X_test, y_train, y_test, binIDs_train, binIDs_test])
