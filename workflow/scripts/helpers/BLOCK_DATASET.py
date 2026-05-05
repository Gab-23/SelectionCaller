import numpy as np
import math as mh
import sklearn.model_selection as sk

def BLOCK_DATASET(X, y, binIDs, BLOCKING):
    ### The function assumes that the dataset is at least chromosome sorted!
    
    RATIO = 2.5

    chromosomes = [int(elem.split("_")[0]) for elem in binIDs.tolist()]
    bin_nums = [int(elem.split("_")[1]) for elem in binIDs.tolist()]

    idxs_train_list = list()
    idxs_test_list = list()

    unique_chromosomes = sorted([int(x) for x in list(set(chromosomes))])

    for chr in unique_chromosomes:

        idxs = np.array([idx for idx, elem in enumerate(chromosomes) if elem == chr])

        mask = np.tile(np.repeat(np.array([True,False]), [BLOCKING,int(BLOCKING/RATIO)], axis = 0), mh.ceil(len(idxs) / BLOCKING))
        mask = mask[0:len(idxs)]
        not_mask = ~mask

        idxs_train = idxs[mask]
        idxs_test = idxs[not_mask]

        idxs_train_list.append(idxs_train)
        idxs_test_list.append(idxs_test)

    flat_idxs_train = [idx for idx_list in idxs_train_list for idx in idx_list]
    flat_idxs_test = [idx for idx_list in idxs_test_list for idx in idx_list]

    X_train = X.iloc[flat_idxs_train]
    y_train = y.iloc[flat_idxs_train]
    binIDs_train = binIDs.iloc[flat_idxs_train]

    X_test = X.iloc[flat_idxs_test]
    y_test = y.iloc[flat_idxs_test]
    binIDs_test = binIDs.iloc[flat_idxs_test]

    return([X_train, X_test, y_train, y_test, binIDs_train, binIDs_test])