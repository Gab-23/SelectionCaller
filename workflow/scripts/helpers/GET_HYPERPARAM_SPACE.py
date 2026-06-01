import scipy as sp
import skopt as sko 

def GET_HYPERPARAM_SPACE(MODEL):
    
    hyperparams_dictionary = {
        "LinearRegression": {
            "GridSearchCV": {
                "model__fit_intercept": [True, False],
                "model__positive": [True, False]
            },
            "RandomSearchCV": {
                "model__fit_intercept": [True, False],
                "model__positive": [True, False]
            },
            "BayesSearchCV": {
                "model__fit_intercept": sko.space.Categorical([True, False]),
                "model__positive": sko.space.Categorical([True, False])
            }
        },
        "SVR": {
            "GridSearchCV": {
                "model__kernel": ["rbf", "linear", "poly"],
                "model__C": [0.1, 1, 10, 100],
                "model__gamma": ["scale", "auto"],
                "model__epsilon": [0.01, 0.1, 0.2, 0.5]
            },
            "RandomSearchCV": {
                "model__kernel": ["rbf", "linear", "poly"],
                "model__C": sp.stats.loguniform(0.1, 100),
                "model__gamma": sp.stats.loguniform(1e-4, 1),
                "model__epsilon": sp.stats.uniform(0.01, 0.5)
            },
            "BayesSearchCV": {
                "model__kernel": sko.space.Categorical(["rbf", "linear", "poly"]),
                "model__C": sko.space.Real(0.1, 100, prior="log-uniform"),
                "model__gamma": sko.space.Real(1e-4, 1, prior="log-uniform"),
                "model__epsilon": sko.space.Real(0.01, 0.5)
            }
        },
        "KNeighborsRegressor": {
            "GridSearchCV": {
                "model__n_neighbors": [3, 5, 7, 9, 11, 15],
                "model__weights": ["uniform", "distance"],
                "model__p": [1, 2]
            },
            "RandomSearchCV": {
                "model__n_neighbors": sp.stats.randint(3, 30),
                "model__weights": ["uniform", "distance"],
                "model__p": [1, 2]
            },
            "BayesSearchCV": {
                "model__n_neighbors": sko.space.Integer(3, 30),
                "model__weights": sko.space.Categorical(["uniform", "distance"]),
                "model__p": sko.space.Categorical([1, 2])
            }
        },
        "XGBRegressor":{
            "GridSearchCV": {
                'model__learning_rate': [0.01, 0.05, 0.1, 0.2],
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
            "RandomSearchCV":{
                'model__learning_rate': sp.stats.uniform(0.01, 0.29),
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
            "BayesSearchCV": {
                "model__learning_rate": sko.space.Real(0.01, 0.30),
                "model__max_depth": sko.space.Integer(5, 10),
                "model__min_child_weight": sko.space.Integer(20, 40),
                "model__gamma": sko.space.Real(0.25, 1.0, prior="log-uniform"),
                "model__subsample": sko.space.Real(0.75, 1.0),
                "model__colsample_bytree": sko.space.Real(0.7, 1.0),
                "model__colsample_bylevel": sko.space.Real(0.7, 1.0),
                "model__reg_alpha": sko.space.Real(1e-3, 10, prior="log-uniform"),
                "model__reg_lambda": sko.space.Real(0.1, 100, prior="log-uniform"),
                "model__n_estimators": sko.space.Integer(200, 300)
            }
        },
        "RandomForestRegressor": {
            "GridSearchCV": {
                "model__n_estimators": [100, 200, 300],
                "model__max_depth": [None, 5, 10, 20],
                "model__min_samples_split": [2, 5, 10],
                "model__min_samples_leaf": [1, 2, 4],
                "model__max_features": ["sqrt", "log2"],
                "model__bootstrap": [True, False]
            },
            "RandomSearchCV": {
                "model__n_estimators": sp.stats.randint(100, 500),
                "model__max_depth": sp.stats.randint(3, 30),
                "model__min_samples_split": sp.stats.randint(2, 20),
                "model__min_samples_leaf": sp.stats.randint(1, 10),
                "model__max_features": ["sqrt", "log2", None],
                "model__bootstrap": [True, False]
            },
            "BayesSearchCV": {
                "model__n_estimators": sko.space.Integer(100, 500),
                "model__max_depth": sko.space.Integer(3, 30),
                "model__min_samples_split": sko.space.Integer(2, 20),
                "model__min_samples_leaf": sko.space.Integer(1, 10),
                "model__max_features": sko.space.Categorical(["sqrt", "log2", None]),
                "model__bootstrap": sko.space.Categorical([True, False])
            }
        }
    }

    hyperparams = hyperparams_dictionary[MODEL]

    return hyperparams