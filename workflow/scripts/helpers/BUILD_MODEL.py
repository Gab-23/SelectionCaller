from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor

def BUILD_MODEL(MODEL, N_JOBS=1):
    # NOTE (fix): n_jobs was hardcoded to 1 for every parallelizable model - parametrized per
    # your request. Note that 1 remains a sensible default: this model sits inside an outer
    # hyperparameter search / cross-validation that already parallelizes across folds and
    # candidates, so letting the model itself also spawn threads risks oversubscription: N_JOBS
    # is now a real, callable-controlled parameter rather than a silently fixed value.

    if MODEL == "XGBRegressor":
        model = XGBRegressor(n_jobs=N_JOBS)
    elif MODEL == "KNeighborsRegressor":
        model = KNeighborsRegressor(n_jobs=N_JOBS)
    elif MODEL == "LinearRegression":
        model = LinearRegression()
    elif MODEL == "SVR":
        model = SVR()
    elif MODEL == "RandomForestRegressor":
        model = RandomForestRegressor(n_jobs=N_JOBS)
    else:
        raise ValueError(f"Unknown MODEL: {MODEL}")

    return model
