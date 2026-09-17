# The calculation in your comment is attempting to calculate an approximate "robust" version of a z-score, 
# but with a center based off using the median instead of the mean 
# and a spread based off using median absolute deviation from the median instead of the standard deviation.
from scipy.stats import norm

def COMPUTE_MAD_Z_PVALS(residuals):
    median_R = residuals.median()
    mad_R = (residuals - median_R).abs().median()
    z_robust = (residuals - median_R) / (1.4826 * mad_R)
    abs_z = z_robust.abs()
    p_values = 2 * (1 - norm.cdf(abs_z))
    return p_values, z_robust