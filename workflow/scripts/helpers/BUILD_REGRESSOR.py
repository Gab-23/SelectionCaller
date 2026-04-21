from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from xgboost import XGBRegressor

def BUILD_MODEL(MODEL):
    
    models_dictionary = {
        "LinearRegression":LinearRegression,
        "SVR":SVR,
        "KNeighborsRegressor":KNeighborsRegressor,
        "XGBRegressor":XGBRegressor
    }

    model = models_dictionary[MODEL]()

    return model