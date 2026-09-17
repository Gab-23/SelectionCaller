import pandas as pd

def GET_OUT_OF_FOLD_PREDICTIONS(MODEL, DATA, binIDs):
    predictions = MODEL.predict(DATA)
    return pd.DataFrame({"prediction": predictions}, index=binIDs.values)