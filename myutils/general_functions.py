# Standard library imports
import os
import re

# Third-party imports
import numpy as np
from scipy.interpolate import interp1d
from PyPDF2 import PdfMerger

def sanityCheck():
    print("general_functions.py Sanity Checked")

def interpolateData(data, newLength):
    """
    Description
    -----------
    interpolateData() will interpolate a single 1 dimensional array up to a certain number of points, assuming that (# desired points) > (len(data)).

    Parameters
    ----------
    data : array-like, shape (1,)
       Array of data destined to be interpolated.
    newLength : int
        Desired number of data points that "data" will be interpolated up to.

    Returns
    -------
    New interperated data array
    """
    # Create a linearly spaced array x_old representing the original indices of the input array
    x_old = np.linspace(0, 1, len(data))

    # Create a linearly spaced array x_new representing the new indices for the resampled array
    x_new = np.linspace(0, 1, newLength)

    # Create an interpolation function f based on the original array and indices
    f = interp1d(x_old, data, kind='linear')

    # Use the interpolation function to calculate the values of the resampled array at the new indices
    return f(x_new)

def combinePDFs(directory = r'.\Figures (temp)', filename="combined.pdf"):
    """
    Description
    -----------
    Combines all of the pdf's in a target directory into a single pdf ordered by number.

    Parameters
    ----------
    directory : string
        The directory containing the pdf's to be merged. Default argument is a file that comes along with the github.
    filename : string
        Name of the new combined .pdf file.

    Returns
    -------
    Nothing. But! This function does produce the combined pdf inside of the target directory.
    """
    # Create a PdfMerger object
    merger = PdfMerger()

    # Get all PDF files in the directory
    pdf_files = [f for f in os.listdir(directory) if f.endswith('.pdf')]

    # Sort files by the numerical value in their names
    def extract_number(filename):
        match = re.search(r'\d+', filename)  # Find the first number in the filename
        return int(match.group()) if match else 0  # Return 0 if no number is found

    pdf_files.sort(key=lambda x: extract_number(x))

    # Merge each PDF file
    for pdf in pdf_files:
        pdf_path = os.path.join(directory, pdf)
        print(f"Adding {pdf}...")
        merger.append(pdf_path)

    # Write the combined PDF to a file
    output_path = os.path.join(directory, filename)
    merger.write(output_path)
    merger.close()

    print(f"All PDFs combined into {output_path}")

def measureModelIR(current, start_times, search_window, interstim, rate):
    """
    Description
    -----------
    Finds the maximum value in a given time interval after a stimulus. Can be used to find many peaks over many intervals.

    Parameters
    ----------
    current : numpy array, shape (1,)
        Activity from the model.
    start_times : numpy array, shape (1,)
        Times the stimulus comes on.
    search_window : float
        Time window (in seconds) after the entry given from 'start_times' that we are going to use to search for baseline period and transient responses (search window).
    interstim : float
        Time between stimuli.
    rate : float
        Simulation timestep.

    Returns
    -------
    transient, persistent, prestimbaselines : tuple
        transient : numpy array, shape (1,)
            Array of transient responses predicted from the model.
        persistent : numpy array, shape (1,)
            Array of the persistent responses predicted from the model.
        prestimbaselines : numpy array, shape (1,)
            Array of the current the step before the stimulus goes. We call this the baseline value.
    """
    start_indices = [int(np.floor(start_time/rate)) for start_time in start_times]
    searchindices = int(np.ceil(search_window/rate))
    interstim_indices = int(np.floor(interstim/rate))

    transient = [0]*len(start_indices)

    persistent = [0]*len(start_indices)

    prestimbaselines = [0]*len(start_indices)

    darkbaseline = np.mean(current[0:searchindices])

    for i in range(len(start_indices)):

        prestimbaseline = current[start_indices[i] - 1] # changed 20250530. Was taking the mean in a window before, but this gets messed up for probe flashes

        #print(prestimbaseline)

        transient[i] = np.max(current[start_indices[i]:start_indices[i] + searchindices]) - prestimbaseline

        prestimbaselines[i] = prestimbaseline

        #print(transient[i])

        persistent[i] = np.mean(current[start_indices[i] - searchindices + interstim_indices:start_indices[i] + interstim_indices]) - darkbaseline # Measure persistent activity right before the next pulse

    return transient, persistent, prestimbaselines

def removeNanCells(dataset_list):
    """
    Description
    -----------
    Remove cells from a data frame that have NaN values.

    Parameters
    ----------
    dataset_list : python list (1,)
        The python list generated from getDataset().

    Returns
    -------
    dataset_list_no_nan : python list (1,)
        Dataset list with all non-NaN entries.
        len(dataset_list) >= len(dataset_list_no_nan).
    """

    dataset_list_no_nan = []
    for dataset in dataset_list:
        dataset_no_nan = []
        for cell in dataset:
            if np.isnan(cell).any() == False:
                dataset_no_nan.append(cell)
        dataset_list_no_nan.append(dataset_no_nan)

    return dataset_list_no_nan

def syncTimeseries(data_synch, times_synch, data_ref, times_ref, method='linear'):
    """
    Description
    -----------
    Interpolates data_synch to match the timebase of data_ref.

    Parameters
    ----------
        times_synch : shape (1,) array
            Time points for data_synch (to be synced).
        data_synch : shape (1,) array
            Data to be interpolated.
        times_ref : shape (1,) array
            Time points for data_ref (target timebase).
        data_ref : shape (1,) array
            Reference data (unchanged).
        method : shape (1,) array
            Interpolation method (e.g., 'linear', 'nearest', 'cubic').

    Returns
    -------
        synced_data : shape (1,) array
            data_synch interpolated to times_ref (1D array).
        times_ref : shape (1,) array
            The unchanged timebase.
        data_ref : shape (1,) array
            The unchanged reference data.
    """
    # Convert inputs to numpy arrays and ensure 1D
    times_synch = np.asarray(times_synch).ravel()
    data_synch = np.asarray(data_synch).ravel()
    times_ref = np.asarray(times_ref).ravel()
    data_ref = np.asarray(data_ref).ravel()

    # Create interpolator for data_synch
    interpolator = interp1d(times_synch, data_synch, kind=method, bounds_error=False, fill_value=np.nan)

    # Interpolate data_synch onto the timebase of data_ref
    synced_data = interpolator(times_ref)

    return synced_data