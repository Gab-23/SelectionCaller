import sklearn as sk

def GENERATE_KFOLD_IDXS(BLOCKING, OUTER_CV, X_train):

    if BLOCKING == "no_blocking":
           kf = KFold(n_splits=OUTER_CV)
           idxs = kf.split(X_train) # dictionary over which I can iterate
       else:
           idxs = BLOCKED_KFOLD(X_train, OUTER_CV)

    return(idxs)