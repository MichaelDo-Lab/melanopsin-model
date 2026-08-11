# Standard library imports

# Third-party imports
import numpy as np

def sanityCheck():
    print("nonlinearities Sanity Checked")

def weib(x, k=15.0, lam=0.35):
    """
    Description
    -----------
    CDF of the weibull distribution https://en.wikipedia.org/wiki/Weibull_distribution

    Parameters
    ----------
    x : float
        Input (just a number; this function is not vectorized).
    k : float
        Shape parameter. Increasing k makes the graph more steep.
    lambda : float
        Scale parameter. Scales x axis.

    Returns
    -------
    g : float
        Output of the Weibull CDF clipped to avoid tiny negative rounding
    """ 

    f = 1 - np.exp(- (x/lam)**k)
    g = np.clip(1.0 - f, 0.0, 1.0) # clip to avoid tiny negative rounding or >1
    return g

def sigmoid(x, k = 50, L = 0.99, x0 = 0.1):
    """
    Description
    -----------
    3pl sigmoid for adding nonlinearity to the model. Full description found at https://en.wikipedia.org/wiki/Logistic_function

    Parameters
    ----------
    x : float
        Input (just a number; this function is not vectorized).
    k : float
        Logistic growth rate (steepness of the curve).
    L : float
        Carrying capacity, or the maximum height of the curve (its asymptote).
    x0 : float
        "Midpoint" of the curve, or its inflection point, or the value halfway to maximum.

    Returns
    -------
    Sigmoidal value at point x (float).
    """
    a = k*(x-x0)

    if a<0:
      return L * np.exp(a) / (1+np.exp(a))
    else:
      return L / (1 + np.exp(-a))