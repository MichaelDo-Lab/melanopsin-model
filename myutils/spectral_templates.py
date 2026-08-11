# Third-party imports
import numpy as np


def gov(wlen, lambdamax):
    """
    Description
    -----------
    Calculate the Govardovskii nomogram with beta band for opsins.

    Parameters
    ----------
    wlen : array, shape (1,)
        Array of wavelengths over which to calculate the nomogram
    lambdamax : float
        Peak wavelength for nomogram.

    Returns
    -------
    nomo : array, shape (1,)
        The computed nomogram, same shape as wlen.
    """
    nomoA = 69.7
    nomoB = 28
    nomoC = -14.9
    nomoD = 0.674
    nomo_a = 0.8795+0.0459*np.exp(-1*(lambdamax-300)**2/11940)
    nomo_b = 0.922
    nomo_c = 1.104
    beta_mb = 189 + 0.315 * lambdamax
    beta_b = -40.5 + 0.195 * lambdamax

    nomo = 1/((np.exp(nomoA*(nomo_a-(lambdamax/wlen)))) + (np.exp(nomoB*(nomo_b - (lambdamax/wlen)))) + (np.exp(nomoC*(nomo_c - (lambdamax/wlen))) + nomoD)) + 0.26 * np.exp(-((wlen-beta_mb)/beta_b)**2)

    return nomo


def customNomogram(wlen, params):
    """
    Description
    -----------
    Nomogram with fully customizable parameters. Equation is of the form of Govardovskii 2000 

    Parameters
    ----------
    wlen : array shape (1,) 
        Wlens over which to calculate nomogram
    params: python dict
        Dictionary of parameter values. Must have the keys: 'true_lambda_max','lambda_max_param','A','B','C','D','a','b','c' and numerical values for each 

    Output
    ------
    nomo : array shape (1,) 
        Return values of the nomogram.
    """

    # Check params dictionary for required values 
    required = ['true_lambda_max','lambda_max_param','A','B','C','D','a','b','c']
    missing = [k for k in required if k not in params]
    if missing:
        raise ValueError(f"Missing parameters: {missing}")

    # Extract parameters
    [lambdamax_true, lambdamax_param, nomoA, nomoB, nomoC, nomoD, nomo_a, nomo_b, nomo_c] = \
        [params[k] for k in [
            'true_lambda_max',
            'lambda_max_param',
            'A','B','C','D','a','b','c'
        ]]

    # Beta band parameters are based on true lambda max, not the parameter value 
    beta_mb = 189 + 0.315 * lambdamax_true
    beta_b = -40.5 + 0.195 * lambdamax_true

    # Find the peak value, then divide by this to normalize peak to 1
    wlen_range = np.linspace(400,600,num=10000)
    raw = lambda wlen : 1 / ((np.exp(nomoA*(nomo_a-(lambdamax_param/wlen)))) + (np.exp(nomoB*(nomo_b - (lambdamax_param/wlen)))) + (np.exp(nomoC*(nomo_c - (lambdamax_param/wlen))) + nomoD)) + 0.26 * np.exp(-((wlen-beta_mb)/beta_b)**2)
    normval = np.max(raw(wlen_range))

    nomo = (1 / ((np.exp(nomoA*(nomo_a-(lambdamax_param/wlen)))) + (np.exp(nomoB*(nomo_b - (lambdamax_param/wlen)))) + (np.exp(nomoC*(nomo_c - (lambdamax_param/wlen))) + nomoD)) + 0.26 * np.exp(-((wlen-beta_mb)/beta_b)**2)) / normval

    return nomo

def getNomogramParams(model_params, pigment):
    """
    Description
    -----------
    Extracts the parameters for the nomogram of the pigment of interest.

    Parameters
    ----------
    model_params : phython dict
        Parameter dictionary of the model.
    pigment: string
        Either 'R', 'M', or 'E'.

    Output
    ------
    Nomogram parameters for a specific pigment.
    """
    keys = [
        'true_lambda_max',
        'lambda_max_param',
        'A','B','C','D','a','b','c'
    ]

    return {k: model_params[f"{pigment}_{k}"] for k in keys}