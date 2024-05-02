# calculations.py
import numpy as np
from scipy.stats import binom, norm

def logit(p):
    """
    Calculate the logit transformation of a probability.
    
    Args:
    p (float): Probability value (0 < p < 1).
    
    Returns:
    float: The logit of the probability.
    """
    return np.log(p / (1 - p))

def inv_logit(l):
    """
    Calculate the logistic function, the inverse of the logit transformation.
    
    Args:
    l (float): Logit value.
    
    Returns:
    float: Probability corresponding to the logit value.
    """
    return np.exp(l) / (1 + np.exp(l))

def calculate_estimates(x11, x10, x01, x00, prevalence, confidence=0.95):
    """
    Calculate diagnostic test performance metrics including sensitivity, specificity, PPV, NPV, and their confidence intervals.
    
    Args:
    x11 (int): Count of true positives.
    x10 (int): Count of false positives.
    x01 (int): Count of false negatives.
    x00 (int): Count of true negatives.
    prevalence (float): Disease prevalence in the population.
    confidence (float, optional): Confidence level for interval estimates, default is 0.95.
    
    Returns:
    dict: Dictionary with Sensitivity, Specificity, PPV, NPV, and their confidence intervals.
    """
    # Calculate sensitivity and specificity
    Se = x11 / (x11 + x01)  # Sensitivity
    Sp = x00 / (x10 + x00)  # Specificity

    # Calculate positive and negative predictive values
    PPV = Se * prevalence / (Se * prevalence + (1 - Sp) * (1 - prevalence))
    NPV = Sp * (1 - prevalence) / ((1 - Se) * prevalence + Sp * (1 - prevalence))

    n1 = x11 + x01  # Total true condition positive cases
    n0 = x10 + x00  # Total true condition negative cases

    # Calculate confidence intervals for Sensitivity and Specificity using the Clopper-Pearson exact method
    Se_CI = binom.interval(confidence=confidence, n=n1, p=Se)
    Sp_CI = binom.interval(confidence=confidence, n=n0, p=Sp)
    Se_CI = (Se_CI[0] / n1, Se_CI[1] / n1)  # Adjust CI to proportion format
    Sp_CI = (Sp_CI[0] / n0, Sp_CI[1] / n0)

    # Calculate the logit transformations for PPV and NPV
    logit_PPV = logit(PPV)
    logit_NPV = logit(NPV)

    # Calculate variance of logit(PPV) and logit(NPV) using the delta method
    Var_logit_PPV = (1 / (Se * (1 - Se) * n1) + 1 / ((1 - Sp) * Sp * n0))
    Var_logit_NPV = (1 / ((1 - Se) * Se * n1) + 1 / (Sp * (1 - Sp) * n0))

    # Calculate confidence intervals for logit(PPV) and logit(NPV)
    z = norm.ppf(1 - (1-confidence) / 2)  # Z-score for the specified confidence level
    logit_PPV_CI = (logit_PPV - z * np.sqrt(Var_logit_PPV), logit_PPV + z * np.sqrt(Var_logit_PPV))
    logit_NPV_CI = (logit_NPV - z * np.sqrt(Var_logit_NPV), logit_NPV + z * np.sqrt(Var_logit_NPV))

    # Convert logit confidence intervals back to the probability scale
    PPV_CI = (inv_logit(logit_PPV_CI[0]), inv_logit(logit_PPV_CI[1]))
    NPV_CI = (inv_logit(logit_NPV_CI[0]), inv_logit(logit_NPV_CI[1]))

    # Return all calculated metrics and their confidence intervals in a dictionary
    return {
        'Sensitivity': Se, 'Specificity': Sp,
        'Sensitivity CI': Se_CI, 'Specificity CI': Sp_CI,
        'PPV': PPV, 'NPV': NPV,
        'PPV CI': PPV_CI, 'NPV CI': NPV_CI
    }
