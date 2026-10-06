def coerce_hyperparam_value(v):
    if v in ("TRUE", "FALSE"):
        return v == "TRUE"
    try:
        return int(v)
    except (ValueError, TypeError):
        pass
    try:
        return float(v)
    except (ValueError, TypeError):
        pass
    return v