import numpy as np

def lightTimingsCheck(lightTimings):
    if len(lightTimings) % 2 != 0:
        raise ValueError("Length of 'lightTimings' is not even. 'lightTimings' must be a list of tuples or an even lengthed array-like object.")
    if isinstance(lightTimings, np.ndarray) != True:
        raise TypeError("Incorrect data type for 'lightTimings'. 'lightTimings' must be a list of tuples or an even lengthed array-like object.")
    
def check_array_length(arr, name, expected_length):
    """
    Checks if a numpy array has the expected length.
    Raises a ValueError if not.
    
    Parameters:
        arr: np.ndarray
        name: str, name of the array for error messages
        expected_length: int, required length
    """
    
    if len(arr) != expected_length:
        raise ValueError(f"{name} must have length {expected_length}, but got {len(arr)}")
    
def check_array_bounds(arr, name, min_val=0, max_val=1):
    """
    Checks if all elements in a 1D NumPy array are within [min_val, max_val].
    Raises a ValueError if any value is out of bounds.
    
    Parameters:
        arr : np.ndarray
        name : str, name of the array for error messages
        min_val : minimum allowed value (inclusive)
        max_val : maximum allowed value (inclusive)
    """
    if np.any(arr < min_val) or np.any(arr > max_val):
        raise ValueError(f"{name} generated a value outside [{min_val}, {max_val}]. Check 'rate' parameter. Large rates can cause the model to blow up.")

def intWlenLtChecker(intensity, wlen, lightTimings):
    # Define expected lengths
    expected_lengths = {'intensity': len(intensity), 'wlen': len(intensity), 'lightTimings': 2*len(intensity)}
    
    # Check all arrays
    check_array_length(intensity, 'intensity', expected_lengths['intensity'])
    check_array_length(wlen, 'wlen', expected_lengths['wlen'])
    check_array_length(lightTimings, 'lightTimings', expected_lengths['lightTimings'])

def is_array_like(x):
    return isinstance(x, (list, np.ndarray))

def is_number(x):
    return isinstance(x, (int, float, np.integer, np.floating))

def intensityChecker(intensities):
    if is_array_like(intensities) != True:
        raise ValueError("'intensities' is wrong data type. 'intensities' must be a (1,) list or numpy array.")
    
def wlenChecker(wlen):
    if is_array_like(wlen) != True:
        raise ValueError("'wlen' is wrong data type. 'wlen' must be a list or numpy array.")
    
def tiChecker(ti):
    if is_number(ti) != True:
        raise ValueError("'ti' is wrong data type. 'ti' must be an int or float.")
    if ti < 0:
        raise ValueError("'ti' can not be negative.")

def stateBoundsChecker(R, M, E, O, Rprime, Mprime, Eprime, idx=None):
    """
    Checks that state variables are within [0, 1].
    Raises a ValueError if any value is out of bounds.

    Parameters
    ----------
    R, M, E, O, Rprime, Mprime, Eprime : np.ndarray
        State variable arrays (fraction per timepoint).
    idx : int, optional
        If given, only the single timepoint at this index is checked (e.g. use
        idx=i+1 inside the time loop to validate only the current step). If None,
        the entire arrays are checked.
    """
    if idx is not None:
        check_array_bounds(np.atleast_1d(R[idx]), "R")
        check_array_bounds(np.atleast_1d(M[idx]), "M")
        check_array_bounds(np.atleast_1d(E[idx]), "E")
        check_array_bounds(np.atleast_1d(O[idx]), "O")
        check_array_bounds(np.atleast_1d(Rprime[idx]), "Rprime")
        check_array_bounds(np.atleast_1d(Mprime[idx]), "Mprime")
        check_array_bounds(np.atleast_1d(Eprime[idx]), "Eprime")
    else:
        check_array_bounds(R, "R")
        check_array_bounds(M, "M")
        check_array_bounds(E, "E")
        check_array_bounds(O, "O")
        check_array_bounds(Rprime, "Rprime")
        check_array_bounds(Mprime, "Mprime")
        check_array_bounds(Eprime, "Eprime")  
