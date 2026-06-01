import sklearn as sk
import skopt as sko

def OPTIMIZE_HYPERPARAMS(ESTIMATOR, HYPERPARAM_SPACE, OPTIMIZER, X, y, SEED, N_JOBS, N_ITER, CV):
    
    HYPERPARAM_SPACE = HYPERPARAM_SPACE[OPTIMIZER]
    # fix space for BayesSearch    
    if OPTIMIZER == "GridSearchCV":
        
        opti = sk.model_selection.GridSearchCV(estimator = ESTIMATOR, 
                                               param_grid = HYPERPARAM_SPACE,
                                               n_iter = N_ITER, 
                                               cv = CV,
                                               refit = True, 
                                               n_jobs = N_JOBS, 
                                               verbose = 0)

    elif OPTIMIZER == "RandomSearchCV":
        
        opti = sk.model_selection.RandomizedSearchCV(estimator = ESTIMATOR, 
                                                     param_distributions = HYPERPARAM_SPACE,
                                                     n_iter = N_ITER, 
                                                     cv = CV, 
                                                     random_state = SEED,
                                                     refit = True, 
                                                     n_jobs = N_JOBS, 
                                                     verbose = 0)
 
    elif OPTIMIZER == "BayesSearchCV":
        
        opti = sko.BayesSearchCV(estimator = ESTIMATOR, 
                                 search_spaces = HYPERPARAM_SPACE,
                                 n_iter = N_ITER, 
                                 cv = CV, 
                                 random_state = SEED,
                                 refit = True,
                                 n_jobs = N_JOBS,
                                 verbose = 0)
    
    else:

        pass
    
    opti.fit(X, y)
    best_model = opti.best_estimator_
    best_hyperparams = opti.best_params_
    
    return([best_model, best_hyperparams])
