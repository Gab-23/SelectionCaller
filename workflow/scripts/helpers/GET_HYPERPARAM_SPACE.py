import scipy as sp

def GET_HYPERPARAM_SPACE(MODEL):
    
    hyperparams_dictionary = {
        "LinearRegression":{"GridSearchCV": {},
                            "RandomSearchCV":{},
                            "BayesSearchCV":{}
        },
        "SVR":{"GridSearchCV": {},
               "RandomSearchCV":{},
               "BayesSearchCV":{}
        },
        "KNeighborsRegressor":{"GridSearchCV": {},
                               "RandomSearchCV":{},
                               "BayesSearchCV":{}
        },
        "XGBRegressor":{"GridSearchCV": {'model__learning_rate': [0.01, 0.05, 0.1, 0.2],
                                         'model__max_depth': [5, 7, 10],
                                         'model__min_child_weight': [20, 30, 40],
                                         'model__gamma': [0.25, 0.5, 1.0],
                                         'model__subsample': [0.75, 0.85, 1.0],
                                         'model__colsample_bytree': [0.7, 0.85, 1.0],
                                         'model__colsample_bylevel': [0.7, 0.85, 1.0],
                                         'model__reg_alpha': [0.001, 0.1, 10],
                                         'model__reg_lambda': [0.1, 10, 100],
                                         'model__n_estimators': [200, 250, 300]
                                         },
                        "RandomSearchCV":{'model__learning_rate': sp.stats.uniform(0.01, 0.29),
                                          'model__max_depth': sp.stats.randint(5, 10),
                                          'model__min_child_weight': sp.stats.randint(20, 40),
                                          'model__gamma': sp.stats.loguniform(0.25, 1.0),
                                          'model__subsample': sp.stats.uniform(0.75, 0.25),
                                          'model__colsample_bytree': sp.stats.uniform(0.7, 0.3),
                                          'model__colsample_bylevel': sp.stats.uniform(0.7, 0.3),
                                          'model__reg_alpha': sp.stats.loguniform(1e-3, 10),
                                          'model__reg_lambda': sp.stats.loguniform(0.1, 100),
                                          'model__n_estimators': sp.stats.randint(200, 300)
                                         },
                        "BayesSearchCV":{'model__learning_rate': sp.stats.uniform(0.01, 0.29),
                                         'model__max_depth': sp.stats.randint(5, 10),
                                         'model__min_child_weight': sp.stats.randint(20, 40),
                                         'model__gamma': sp.stats.loguniform(0.25, 1.0),
                                         'model__subsample': sp.stats.uniform(0.75, 0.25),
                                         'model__colsample_bytree': sp.stats.uniform(0.7, 0.3),
                                         'model__colsample_bylevel': sp.stats.uniform(0.7, 0.3),
                                         'model__reg_alpha': sp.stats.loguniform(1e-3, 10),
                                         'model__reg_lambda': sp.stats.loguniform(0.1, 100),
                                         'model__n_estimators': sp.stats.randint(200, 300)
                                        }
        }
    }

    hyperparams = hyperparams_dictionary[MODEL]

    return hyperparams