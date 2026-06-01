from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

def BUILD_MODEL(MODEL):

    if MODEL == "XGBRegressor":
        model = XGBRegressor(n_jobs=1)
    elif MODEL == "KNeighborsRegressor":
        model = KNeighborsRegressor(n_jobs=1)
    elif MODEL == "LinearRegression":
        model = LinearRegression()
    elif MODEL == "SVR":
        model = SVR()
    elif MODEL == "RandomForestRegressor":
        model = RandomForestRegressor(n_jobs=1)

    return model