import shap as sh

def GET_EXPLAINER(MODEL, IS_TREE_BASED, IS_LINEAR):

    if IS_TREE_BASED:
        explainer = sh.TreeExplainer(MODEL.named_steps["model"], approximate=False)
    elif IS_LINEAR:
        explainer = sh.LinearExplainer(...)
    else:
        explainer = sh.Explainer(MODEL.named_steps["model"], approximate=False)

    return(explainer)