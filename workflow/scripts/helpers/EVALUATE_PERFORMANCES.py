import numpy as np
import sklearn as sk

def EVALUATE_PERFORMANCES(REGRESSOR, X_list, y_list, metrics_list, CV):

    X_train, X_test = X_list
    y_train, y_test = y_list
    
    performances = {}

    for metric in metrics:

        mean_metric = np.mean(sk.model_selection.cross_val_score(REGRESSOR, X_train, y_train, cv = CV, scoring = metric))
        performances[metric] = mean_metric

    return(performances)