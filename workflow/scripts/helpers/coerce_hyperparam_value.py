import numpy as np
import pandas as pd


def coerce_hyperparam_value(v):
    # Turns one consensus hyperparameter value (as written by select_best_model.R) into the
    # Python type the model expects: "TRUE"/"FALSE" -> bool, "250" -> int, "0.1" -> float,
    # anything else (e.g. "rbf", "depthwise") -> str.
    # NOTE (fix): non-string input is handled explicitly. Previously a value pandas had already
    # parsed went straight into int(): True became 1 (LinearRegression rejects fit_intercept=1)
    # and every fractional float was TRUNCATED, e.g. learning_rate 0.1 -> 0 and subsample 0.85
    # -> 0, so an XGBRegressor chosen as best model was refitted with learning_rate = 0 and
    # silently predicted a constant (no error raised).
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return int(v) if float(v).is_integer() else float(v)
    v = str(v).strip()
    if v.upper() in ("TRUE", "FALSE"):
        return v.upper() == "TRUE"
    try:
        return int(v)          # only succeeds for integer literals: int("0.1") raises
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v


def read_consensus_hyperparams(path):
    # Reads the consensus hyperparameters file (2 columns, no header: name, value) into a
    # {name: value} dict ready for model.set_params(), with the "model__" pipeline prefix removed.
    # dtype=str + keep_default_na=False: values reach coerce_hyperparam_value exactly as R wrote
    # them. Without it pandas guesses the column type from the whole file: a file with only
    # "TRUE" becomes a bool column, one with only numbers a float column (see the NOTE above).
    df = pd.read_csv(path, header=None, sep="\t", dtype=str, keep_default_na=False)
    df.columns = ["hyperparam_name", "value"]
    names = [name.replace("model__", "") for name in df.hyperparam_name.tolist()]
    return {k: coerce_hyperparam_value(v) for k, v in zip(names, df.value)}