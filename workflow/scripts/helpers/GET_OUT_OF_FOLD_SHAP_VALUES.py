def GET_OUT_OF_FOLD_SHAP_VALUES(EXPLAINER, DATA, binIDs):
    values = EXPLAINER.shap_values(DATA)
    return pd.DataFrame(values, index=binIDs.values, columns=DATA.columns)