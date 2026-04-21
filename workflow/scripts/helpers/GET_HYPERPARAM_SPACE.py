import scipy as sp

def GET_HYPERPARAM_SPACE(MODEL):
    
    hyperparams_dictionary = {
        "LinearRegression":{},
        "SVR":{},
        "KNeighborsRegressor":{},
        "XGBRegressor":{}
    }

    hyperparams = hyperparams_dictionary[MODEL]

    return hyperparams
