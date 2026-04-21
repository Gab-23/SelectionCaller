import sklearn as sk
import skopt as sko 
def OPTIMIZE_HYPERPARAMS(REGRESSOR, HYPERPARAM_SPACE, OPTIMIZER, X_train, y_train, CV, SEED, N_JOBS, N_ITER):
    
    if OPTIMIZER == "GridSearchCV":
        
        opti = sk.model_selection.GridSearchCV(estimator = REGRESSOR, search_spaces = HYPERPARAM_SPACE,
                                               n_iter = N_ITER, cv = CV, random_state = SEED,
                                               refit = True, n_jobs = N_JOBS, verbose = 0)

    elif OPTIMIZER == "RandomSearchCV":
        
        opti = sk.model_selection.RandomizedSearchCV(estimator = REGRESSOR, search_spaces = HYPERPARAM_SPACE,
                                                     n_iter = N_ITER, cv = CV, random_state = SEED,
                                                     refit = True, n_jobs = N_JOBS, verbose = 0)
 
    elif OPTIMIZER == "BayesSearchCV":
        
        opti = sko.BayesSearchCV(estimator = REGRESSOR, search_spaces = HYPERPARAM_SPACE,
                                 n_iter = N_ITER, cv = CV, random_state = SEED,
                                 refit = True, n_jobs = N_JOBS, verbose = 0)
    
    else:
        pass
    
    opti.fit(X_train, y_train)

    best_model = opti.best_estimator_
    best_hyperparams = opti.best_params_
    
    return([best_model, best_hyperparams])
