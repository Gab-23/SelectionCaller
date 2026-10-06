import numpy as np
import pandas as pd

def STRINGENCY_TO_PROMINENCE(x_abs, stringency):
    std = np.std(x_abs)
    if std == 0:
        return 0, 0
    prominence = np.quantile(x_abs, 1 - stringency) # compute the stringency-th quantile, which is just K * std(x_abs), where K is prominence / std
    K = prominence / std
    return K, prominence