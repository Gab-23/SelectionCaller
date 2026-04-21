def AGGREGATE_SELECTED_FEATURES(cols, selected_features_dict):

    rows = list(selected_features_dict.keys())

    selected_feats_df = pd.DataFrame(np.zeros((len(rows), len(cols))), index = rows, columns = cols)

    for k,v in selected_features_dict.items():
        for feature in v:
            selected_feats_df.loc[k,feature] = 1

    return(selected_feats_df)