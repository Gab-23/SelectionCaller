import shap as sh

def GET_EXPLAINER(MODEL, IS_TREE_BASED, IS_LINEAR, BACKGROUND_DATA):
    # NOTE (fix): this branch was never finished - shap.LinearExplainer(...) passed a bare
    # Ellipsis as its only argument, but LinearExplainer requires two positional arguments
    # (model, masker) - confirmed directly: TypeError, missing 1 required positional argument
    # 'masker'. A masker/background dataset was never threaded through to this function at all,
    # so this is a genuine missing parameter, not just a typo - added BACKGROUND_DATA (expected
    # to be the training data used to fit MODEL, e.g. X_train_train) as a new argument.

    if IS_TREE_BASED:
        explainer = sh.TreeExplainer(MODEL.named_steps["model"], approximate=False)
    elif IS_LINEAR:
        explainer = sh.LinearExplainer(MODEL.named_steps["model"], masker=BACKGROUND_DATA)
    else:
        explainer = sh.Explainer(MODEL.named_steps["model"], BACKGROUND_DATA)

    return(explainer)
