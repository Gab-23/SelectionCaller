import numpy as np
import pandas as pd

def SMOOTH_SIGNAL(values, smoothing_size, center=True, min_periods=1):
    return pd.Series(values).rolling(window=smoothing_size, center=center, min_periods=min_periods).mean().values