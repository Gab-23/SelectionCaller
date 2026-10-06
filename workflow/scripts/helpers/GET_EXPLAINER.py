import numpy as np
import pandas as pd
import shap as sh

class _PreprocessingExplainerWrapper:
    def __init__(self, explainer, preprocessor):
        self._explainer = explainer
        self._preprocessor = preprocessor
        self.expected_value = explainer.expected_value

    def shap_values(self, X):
        X_transformed = self._preprocessor.transform(X)
        return self._explainer.shap_values(X_transformed)


def GET_EXPLAINER(MODEL, IS_TREE_BASED, IS_LINEAR, BACKGROUND_DATA):
    if IS_TREE_BASED:
        explainer = sh.TreeExplainer(MODEL.named_steps["model"], approximate=False)  # no background provided, so feature_perturbation = "tree_path_dependent"
    elif IS_LINEAR:
        preprocessor = MODEL[:-1]
        background_transformed = preprocessor.transform(BACKGROUND_DATA)
        raw_explainer = sh.LinearExplainer(MODEL.named_steps["model"], masker=background_transformed)
        explainer = _PreprocessingExplainerWrapper(raw_explainer, preprocessor)
    else:
        explainer = sh.Explainer(MODEL.predict, BACKGROUND_DATA)
        # NOTE (speed-up): the base value is the mean model output over the masker's background
        # sample. It was obtained as explainer(BACKGROUND_DATA).base_values[0], i.e. by computing
        # full SHAP values for the whole training set (~4x the SHAP work actually needed, for
        # one number). Computed directly here; verified equal to machine precision.
        background = pd.DataFrame(explainer.masker.data, columns=BACKGROUND_DATA.columns)
        explainer.expected_value = float(np.mean(MODEL.predict(background)))
    return explainer