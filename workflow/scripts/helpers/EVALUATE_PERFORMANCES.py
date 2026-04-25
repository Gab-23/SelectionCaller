import numpy as np
import sklearn as sk

def adj_r2(REGRESSOR, X, y):
    r2 = sk.metrics.r2_score(y_true =  y, y_pred = REGRESSOR.predict(X))
    n = X.shape[0]
    p = X.shape[1]
    return 1-(1-r2)*(n-1)/(n-p-1)

def EVALUATE_PERFORMANCES(REGRESSOR, X, y, metrics_list):

    performances = {}

    for metric in metrics_list:

        if metric == "adj_r2":
            evaluator = adj_r2
        else:
            evaluator = sk.metrics.get_scorer(metric)

        performances[metric] = evaluator(REGRESSOR, X, y)
        
    return(performances)