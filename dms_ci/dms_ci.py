# import beta distribution
from scipy.stats import beta
import numpy as np

def dms_ci(p, n, alpha = 0.95):
    """Provides confidence intervals for DMS-MaPseq data.

    Parameters
    ----------

    p : array_like
        Number of mutations for each position.

    n : array_like
        Number of reads for each position.

    alpha : float, optional
        Significance level of the confidence interval. Default is 0.05.

    Returns
    -------

    low : array_like
        Lower confidence interval.

    high : array_like
        Upper confidence interval.


    """
    # formatting
    p = np.array(p)
    n = np.array(n)
    
    # compute the distribution beta parameters
    a = p*n
    b = (1-p)*n
    
    # compute the confidence interval
    low, high = beta.ppf(0.5 - alpha/2, a, b), beta.ppf(0.5 + alpha/2, a, b)
    
    # if p = 0, the bounds are 0
    if hasattr(low, '__iter__'):
        low[p == 0] = 0
        high[p == 0] = 0
        
    # cap the confidence interval to the [0, 1] range. 
    low = np.clip(low, 0, 1)
    if hasattr(low, '__iter__'):
        low[np.isnan(low)] = 0
    high = np.clip(high, 0, 1)
    if hasattr(high, '__iter__'):
        high[np.isnan(high)] = 1
    
    return low, high