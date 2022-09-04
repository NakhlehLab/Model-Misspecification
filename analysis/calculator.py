from scipy.stats import chisquare
import numpy as np
import math

def internal_time(max_prob):
    t = -np.log(1.5*(1-max_prob))
    return t

def compute_exp(t):
    expct = []
    x = math.exp(-t)
    expct.append(1-2/3*x)
    expct.append(1/3*x)
    expct.append(1/3*x)
    return expct

# X2, delta degree_of_freedom
# obs, exps, ddof
def compute_Pvalue(obs, exps, ddof):
    chi2, pvalue = chisquare(obs, exps, ddof=ddof)
    return chi2, pvalue
    

