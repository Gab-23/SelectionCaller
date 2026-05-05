import sklearn as sk

def GENERATE_KFOLD_IDXS(BLOCKING, OUTER_CV, X, binIDs):
    # this function assumes that if BLOCKING = "no_blocking", X_train should still at least chromosome sorted
    # This is because the function BLOCK_DATASET should keep this sorting and KFold generates sequential indexes
    if BLOCKING != "no_blocking" and not CHECK_CHR_SORTING(binIDs):
        raise ValueError("BLOCKING != 'no_blocking' and CHECK_CHR_SORTING = False")

    kf = sk.model_selection.KFold(n_splits=OUTER_CV)
    idxs = dict() # dictionary over which I can iterate
    
    for idx, (train_idx, test_idx) in enumerate(kf.split(X)):
        idxs[str(idx)] = (train_idx, test_idx)

    return(idxs)