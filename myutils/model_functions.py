# Standard library imports
import os
import warnings

# Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# tqdm.auto selects the notebook widget in Jupyter and a plain console bar
# everywhere else (including a frozen GUI), so this works in every environment.
from tqdm.auto import tqdm

# Local package imports
from .app_paths import manuscript_dir, manuscript_recordings_dir, manuscript_spectra_dir
from .spectral_templates import gov, customNomogram, getNomogramParams
from .nonlinearities import weib, sigmoid
from .error_functions import *

# Repo `data/manuscript/recordings/` directory, for example recordings used by setupModelRun.
# Resolved via app_paths so it is correct both from source and when frozen.
_RECORDINGS_DIR = str(manuscript_recordings_dir())
_SPECTRA_DIR = str(manuscript_spectra_dir())
_IR_DATA_DIR = str(manuscript_dir() / "irradiance response measurements")


def _pulse_starts_skipping_old_indices(stim_start, interstim, n_old, skip_indices):
    """
    Onset times for pulses that lived on a uniform grid stim_start + k * interstim,
    k = 0 .. n_old - 1, omitting placeholder slots skip_indices (dark = no interval).
    """
    skip = frozenset(skip_indices)
    return np.array(
        [stim_start + k * interstim for k in range(n_old) if k not in skip],
        dtype=float,
    )


def sanityCheckMelanopsinFunction():
    print("Sanity Checked for melanopsin_function")

def powerSpectrumToPhotonFlux(power, wlen):
    """
    Description
    -----------
    Converts a power spectrum to a photon flux.

    Parameters
    ----------
    power : array-like, shape (1,)
        Array of power values in microwatts.
    wlen : array-like, shape (1,)
        Array of wavelengths in nanometers. Must have the same length as 'power'.

    Returns
    -------
    nphotons : array-like, shape (1,)
        Array of photon counts. Same shape as power and wlen.
    """

    if len(power) != len(wlen):
        raise ValueError('power and wlen are not the same length.')

    # Set power measurements < 0 to 0
    power[power < 0] = 0

    h = 6.626*10**(-34) # Planck's constant
    c = 2.9979*10**8 # Speed of light
    single_photon_energy = [h*c / (w*10**(-9)) for w in wlen]
    nphotons = [(power_uw * 10**-6) / ev * 10**-8 for power_uw, ev in zip(power, single_photon_energy)]

    return nphotons


def scalingConstant(baseWavelengths, baseIntensities, targetWavelengths, targetIntensities, lambdaMax = 467):
    """
    Description
    -----------
    Mainly used for monochromatic stimuli
    scalingConstant() will compute the scaling constant needed for a target spectrum of light to match melanopsin photon absorption of a known base spectrum.
    Ex. when matching the xenon to the idealized 440nm light pulse, the xenon was the target spectrum and the 440nm pulse is the base spectrum

    Parameters
    ----------
    baseWavelengths : np.array shape (1,)
    
    baseIntensities : np.array shape (1,)
    
    targetWavelengths : np.array shape (1,)
    
    targetIntensities : np.array shape (1,)
    
    lambdaMax : scalar
        Wavenlength (in nm) of light stimulus.

    Returns
    -------
    scalingConstant : float
        The scaling constant needed for the target spectrum to match melanopsin photon absorption of the base spectrum
    """
    baseWavelengths = baseWavelengths
    targetWavelengths = targetWavelengths

    baseWlenBin = np.mean(np.diff(baseWavelengths))
    targetWlenBin = np.mean(np.diff(targetWavelengths))

    # R state nomogram
    R_params = {
    'true_lambda_max': 467,
    'lambda_max_param': 494,
    'A': 44.582569987408945,
    'B': 8.614863087478009,
    'C': -11.594413782436712,
    'D': 0.0,
    'a': 0.9324464256209672,
    'b': 0.9901137601105081,
    'c': 1.128289328759767,
    }

    nomoBase = customNomogram(baseWavelengths, R_params)
    nomoTarget = customNomogram(targetWavelengths, R_params)

    # Compute the scaling constant, transposing the Govardovskii nomogram to find the dot product
    scalingConstant = (np.dot(baseIntensities, nomoBase.T)*baseWlenBin)/(np.dot(targetIntensities, nomoTarget.T)*targetWlenBin)
    return(scalingConstant)

def scaleToEquivalentPhotons(specwlen, specint, intensities, wlen_match = 440):
    """
    Description
    -----------
    scaleToEquivalentPhotons() normalizes a given spectrum by its photon count and then produces the scaling constant needed for the spectrum to have the have the same effect on the ground state melanopsin as one 440nm monochromatic photon.

    Parameters
    ----------
    specwlen : array-like, shape (1,)
        Wavelengths of light used in the desired spectrum (x-axis).
    specint : int
        Intensities of wavelengths in the desired spectrum. Must have the same length as 'specwlen'.
    intensities : array, shape (1,)
        Intensities in photons/micron^2 of stimulus light pulses.
    wlen_match : float
        Desired wavelength to match to. Default is 440nm out of convenience.

    Returns
    -------
    scaling_constant : float
        scaling constant needed for the spectrum to have the have the same effect on the ground state melanopsin as one 440nm monochromatic photon.
    """
    
    dwlen = np.mean(np.diff(specwlen))
    photon_count = np.trapezoid(specint, specwlen, dwlen)
    norm_specint = specint / photon_count
    norm_intensity = np.zeros(len(specwlen))
    norm_intensity[np.argmin(np.abs(specwlen - wlen_match))] = 1 / dwlen
    scaling_constant = scalingConstant(specwlen, norm_intensity, specwlen, norm_specint)
    scaled_intensities = np.array([scaling_constant * norm_specint * intensity for intensity in intensities])
    wlens = np.tile(specwlen, (len(intensities), 1))

    return wlens, scaled_intensities

def scaleToPhotonCount(specwlen, specint, intensities):
    """
    Description
    -----------
    scaleToPhotonCount scales an input sepectrum so that it has the desired number of photons specified in intensities

    Parameters
    ----------
    specwlen : array-like, shape (1,)
        Range in the wavelength of light used in the spectrum (x-axis).
    specint : int
        Desired number of data points that "data" will be interpolated up to.
    intensities : array, shape (1,)
        Intensities in photons/micron^2 of stimulus light pulses.

    Returns
    -------
    wlens : array, shape (len(intensities),)
        Wavelengths for each spectral frame.
    scaled_intensities : array, shape (len(intensities),len(wlens[i]))
        Intensities for each spectral frame.
    """
    
    dwlen = np.mean(np.diff(specwlen))
    photon_count = np.trapezoid(specint, specwlen, dwlen)
    norm_specint = specint / photon_count
    scaled_intensities = np.array([norm_specint * intensity for intensity in intensities])
    wlens = np.tile(specwlen, (len(intensities), 1))

    return wlens, scaled_intensities

# Selects a light stimulus based on user input. Returns parameters for the model run and eventual comparison to cellular data (for light stimuli used in the study). Extend this function to describe your own light stimuli.
def setupModelRun(experiment, stimulustype, rate = 0.005, dilation = 1, specwlen=[], specint=[]):
    """
    Description
    -----------
    Creates a light stimulus based on user input.

    Parameters
    ----------
    experiment : string
       The specific type of experiment wanted to run. Current support only for "intensity response".
    stimulustype : string
        The specific type of stimulii to be used in the experiment. Current support is too large to list.
    rate : float
        The timestep of the experiment (in seconds). Seems as though if the timestep is higher than t=0.01 the model blows up, so always use t <=0.01 when running.
    dilation : float
        To scale the time step
    specwlen : array-like, shape (1,)
        The spectrum of wavelengths used in a pulse of light if it is not monochromatic.
    specint : array-like, shape (1, len(specwlen))
        The corresponding intensities for specwlen. Same shape as specwlen.

    Returns
    -------
    outDict
        A dictionary containing the following keys:
        stim_start : array, (1,)
            Stimulus start times
        stim_duration : array, (1,)
            Durations of stimuli.
        timings : Numpy array, shape (1,)
            Array of all start times for each pulse in the experiment. Note: It is only the start time, it does not have an end time.
        interstim : Numpy array, shape (1, len(timings))
            Space between the pulses.
        ti : float
            Ttime interval of the experiment. Or, the total length of the experiment.
        intensities : Numpy array, shape (1, len(timings))
            Array representing the spetral intensities for the light pulse.
        datafilename : string
            Name for the saved data file.
        irdatafilenamet : string
            Name for the intensity response data file name transient.
        irdatafilenamep : string
            Name for the intensity response data file name persistent.
        dataidx : string
            data id name.
        wlen : Numpy array, shape (1, len(timings))
            array containing the wavelengths of each light pulse.
        intensitieslog : array, (1,)
            Conditional parameter. Log10 of each of the intensities
    """

    # Initialize empty data file names
    datafilename = ''
    irdatafilenamet = ''
    irdatafilenamep = ''
    dataidx = 0
    intensities_log=[] #give everything a default intensities_log empty list and then just change it if it applies. Otherwise it doesn't get sent in the output dictionary

    if experiment == "intensity response":

        if stimulustype == "440 440":

            intensities = np.array([10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0, 10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0]) # Retinal irradiance for voltage-clamp 440-440 stimuli
            wlen = np.array((np.ones(int(len(intensities) / 2))*440,np.ones(int(len(intensities) / 2))*440)).reshape(-1)
            datafilename = os.path.join(_RECORDINGS_DIR, "V2308154 10 Hz lowpass 440 440.csv") # Example cell is V2308154
            # irdatafilenamet = os.path.join(_IR_DATA_DIR, "Vclamp 440 440 transient.csv")
            # irdatafilenamep = os.path.join(_IR_DATA_DIR, "Vclamp 440 440 persistent.csv")

            dataidx = 10

        elif stimulustype == "560 560":

            intensities = np.array([10**6.8, 10**7.2, 10**7.4, 10**7.8, 10**8.1, 10**8.4, 10**8.8, 10**9.1, 10**9.4, 10**9.6, 0, 10**6.8, 10**7.2, 10**7.4, 10**7.8, 10**8.1, 10**8.4, 10**8.8, 10**9.1, 10**9.4, 10**9.6, 0]) # Retinal irradiance for voltage-clamp 560-560 stimuli
            wlen = np.array((np.ones(int(len(intensities) / 2))*560,np.ones(int(len(intensities) / 2))*560)).reshape(-1)
            datafilename = os.path.join(_RECORDINGS_DIR, "V2206075 10 Hz lowpass 560 560.csv")  # Example cell is V2206075
            irdatafilenamet = os.path.join(_IR_DATA_DIR, "Vclamp 560 560 transient.csv")
            irdatafilenamep = os.path.join(_IR_DATA_DIR, "Vclamp 560 560 persistent.csv")
            dataidx = 3

        elif stimulustype == "440 560":

            intensities = np.array([10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0, 10**6.9, 10**7.2, 10**7.5, 10**7.8, 10**8.1, 10**8.4, 10**8.8, 10**9.1, 10**9.4, 10**9.7, 0]) # Retinal irradiance for voltage-clamp 440-560 stimuli
            wlen = np.array((np.ones(int(len(intensities) / 2))*440,np.ones(int(len(intensities) / 2))*560)).reshape(-1)
            datafilename = os.path.join(_RECORDINGS_DIR, "V2308043 10 Hz lowpass 440 560.csv") # Example cell is V2308043
            irdatafilenamet = os.path.join(_IR_DATA_DIR, "Vclamp 440 560 transient.csv")
            irdatafilenamep = os.path.join(_IR_DATA_DIR, "Vclamp 440 560 persistent.csv")
            dataidx = 8

        elif stimulustype == "560 440":

            intensities = np.array([10**6.8, 10**7.2, 10**7.4, 10**7.8, 10**8.1, 10**8.4, 10**8.8, 10**9.1, 10**9.4, 10**9.6, 0, 10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**8.0, 10**8.3, 10**8.6, 10**9.0, 10**9.3, 10**9.5, 0]) # Retinal irradiance for voltage-clamp for 560-440 stimuli
            wlen = np.array((np.ones(int(len(intensities) / 2))*560,np.ones(int(len(intensities) / 2))*440)).reshape(-1)
            datafilename = os.path.join(_RECORDINGS_DIR, "V2205264 10 Hz lowpass 560 440.csv") # Example cell is V2205264
            irdatafilenamet = os.path.join(_IR_DATA_DIR, "Vclamp 560 440 transient.csv")
            irdatafilenamep = os.path.join(_IR_DATA_DIR, "Vclamp 560 440 persistent.csv")
            dataidx = 3

        elif stimulustype == "440 440 monochromatic":

            dwlen = np.mean(np.diff(specwlen))
            photon_count = np.trapezoid(specint, specwlen, dwlen)
            norm_specint = specint / photon_count
            norm_intensity_440 = np.zeros(len(specwlen))
            norm_intensity_440[np.argmin(np.abs(specwlen - 440))] = 1 / dwlen
            scale_spec_to440 = scalingConstant(specwlen, norm_intensity_440, specwlen, norm_specint)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 440 nm photons: ', scale_spec_to440)
            intensities440 = [10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0, 10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0] # Retinal irradiance for voltage-clamp 440-440 stimuli
            # intensities440 = [10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0, 10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 10**9.7, 0] # Retinal irradiance for voltage-clamp 440-440 stimuli
            intensities = np.array([scale_spec_to440 * norm_specint * intensity for intensity in intensities440])
            wlen = np.tile(specwlen, (len(intensities), 1))

            plt.plot(wlen[0], intensities[0] / np.max(intensities[0]), color='k')
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Normalized Intensity')
            plt.show()

            dataidx = -1

        elif stimulustype == "560 560 monochromatic":

            dwlen = np.mean(np.diff(specwlen))
            photon_count = np.trapezoid(specint, specwlen, dwlen)
            norm_specint = specint / photon_count
            norm_intensity_440 = np.zeros(len(specwlen))
            norm_intensity_440[np.argmin(np.abs(specwlen - 440))] = 1 / dwlen
            scale_spec_to440 = scalingConstant(specwlen, norm_intensity_440, specwlen, norm_specint)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 440 nm photons: ', scale_spec_to440)
            intensities440 = [10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0, 10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0] # Retinal irradiance for voltage-clamp 440-440 stimuli
            intensities = np.array([scale_spec_to440 * norm_specint * intensity for intensity in intensities440])
            wlen = np.tile(specwlen, (len(intensities), 1))

            plt.plot(wlen[0], intensities[0] / np.max(intensities[0]), color='k')
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Normalized Intensity')
            plt.show()

            dataidx = -1

        elif stimulustype == "440 CC":

            intensities = np.array([10**4.73, 10**5.53, 10**6.15, 10**7.0, 10**7.61, 10**8.23, 10**8.93, 10**9.49])  # Retinal irradiance for 440 nm current clamp stimuli
            wlen = np.array(np.ones(int(len(intensities)))*440).reshape(-1)

        elif stimulustype == "560 CC":

            wlen = np.array(np.ones(int(len(intensities)))*560).reshape(-1)
            intensities = np.array([10**4.92, 10**5.72, 10**6.34, 10**7.18, 10**7.80, 10**8.42, 10**9.12, 10**9.68]) # Retinal irradiance for 560 nm current clamp stimuli

        elif stimulustype == "xenon xenon":

            intensities440 = [10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0, 10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0] # Retinal irradiance for voltage-clamp 440-440 stimuli
            intensities = np.array([cXenon * intensity_xenon/photon_count_xenon*intensity for intensity in intensities440])
            wlen = np.tile(wlenspec, (len(intensities), 1))
            dataidx = -1
            plt.plot(wlen[0], intensities[0] / np.max(intensities[0]), color='k')
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Normalized Intensity')
            plt.show()

        elif (
            stimulustype == "xenon 560 suppression xenon"
            or stimulustype == "xenon 560 suppression xenon first set only"
            or stimulustype == "xenon 560 suppression xenon max sup max"
            or stimulustype == "xenon 560 suppression xenon triple pre sup bright"
            or stimulustype == "xenon 560 suppression xenon eye"
            or stimulustype == "xenon"
            or stimulustype == "xenon eye"
        ):  # Xenon


            if stimulustype == 'xenon':
                intensities_log = [
                                    7.433598706, 7.814542969, 8.007575299, 8.413050505, 8.728291025,
                                    9.162075803, 9.448936069, 9.707742823, 10.12226415, 10.35269065,
                                    7.433598706, 7.814542969, 8.007575299, 8.413050505, 8.728291025,
                                    9.162075803, 9.448936069, 9.707742823, 10.12226415, 10.35269065,
                                ]
            if stimulustype == 'xenon eye':
                intensities_log = [
                                    7.477463625, 7.804534741, 8.091395007, 8.350201762, 8.678404421,
                                    9.046801853, 9.265081359, 9.822578977, 10.1378195, 10.24453316,
                                    7.477463625, 7.804534741, 8.091395007, 8.350201762, 8.678404421,
                                    9.046801853, 9.265081359, 9.822578977, 10.1378195, 10.24453316,
                                ]
            if stimulustype == "xenon 560 suppression xenon":
                intensities_log = [
                                    7.433598706, 7.814542969, 8.007575299, 8.413050505, 8.728291025,
                                    9.162075803, 9.448936069, 9.707742823, 10.12226415, 10.35269065,
                                    7.433598706, 7.814542969, 8.007575299, 8.413050505, 8.728291025,
                                    9.162075803, 9.448936069, 9.707742823, 10.12226415, 10.35269065,
                                    9.810921113, 10.35269065
                                ]
            if stimulustype == "xenon 560 suppression xenon first set only":
                # First xenon train + gap, then 560 suppression and final xenon probe (no second train)
                intensities_log = [
                    7.433598706, 7.814542969, 8.007575299, 8.413050505, 8.728291025,
                    9.162075803, 9.448936069, 9.707742823, 10.12226415, 10.35269065,
                    9.810921113, 10.35269065,
                ]
            if stimulustype == "xenon 560 suppression xenon max sup max":
                # Brightest manuscript xenon step, 560 suppression, brightest xenon again (same ISI base)
                intensities_log = [10.35269065, 9.810921113, 10.35269065]
            if stimulustype == "xenon 560 suppression xenon triple pre sup bright":
                # Three brightest xenon pulses, 560 suppression, one brightest xenon probe (same ISI / timing rules)
                intensities_log = [
                    10.35269065,
                    10.35269065,
                    10.35269065,
                    9.810921113,
                    10.35269065,
                ]
            if stimulustype == "xenon 560 suppression xenon eye":
                intensities_log = [
                                    7.477463625, 7.804534741, 8.091395007, 8.350201762, 8.678404421,
                                    9.046801853, 9.265081359, 9.822578977, 10.1378195, 10.24453316,
                                    7.477463625, 7.804534741, 8.091395007, 8.350201762, 8.678404421,
                                    9.046801853, 9.265081359, 9.822578977, 10.1378195, 10.24453316,
                                    9.810921113, 10.24453316
                                ]

            # Count xenon photons
            photon_count_xenon = np.trapezoid(specint, specwlen, np.mean(np.diff(specwlen)))

            # Set the intensities by scaling the integral-normalized spectrum
            intensities = np.array([specint/photon_count_xenon * 10**intensity for intensity in intensities_log])

            # Insert 560 pulse
            if (
                stimulustype == "xenon 560 suppression xenon"
                or stimulustype == "xenon 560 suppression xenon first set only"
                or stimulustype == "xenon 560 suppression xenon max sup max"
                or stimulustype == "xenon 560 suppression xenon triple pre sup bright"
                or stimulustype == "xenon 560 suppression xenon eye"
            ):

                # Load the spectrum
                spec560 = pd.read_table(
                    os.path.join(_SPECTRA_DIR, "560bp Viet.txt"),
                    names=["wlen", "power"],
                )

                # Extract wavelength
                wlen560 = spec560["wlen"].to_numpy()

                # Convert power to photon flux ('Intensity')
                intensity560 = np.array(powerSpectrumToPhotonFlux(spec560["power"], wlen560))

                # Set 560 spectrum to near-zero close to the edges of the spectrometers range (where there is measurement noise)
                intensity560[wlen560 < 500] = 0.001
                intensity560[wlen560 > 600] = 0.001

                # Create the interpolation function
                f = interp1d(wlen560, intensity560, kind='cubic')

                # Interpolate intensity560 at the xenon spectrum's wavelength values
                intensity560_interp = f(specwlen)

                # Count 560 photons
                photon_count_560 = np.trapezoid(intensity560_interp, specwlen, np.mean(np.diff(specwlen)))

                # Insert the 560 pulse
                idx_560 = -2 # Index of the 560 nm pulse
                intensities[idx_560] = intensity560_interp/photon_count_560 * 10**intensities_log[idx_560] # Insert the 560 spectrum

            wlen = np.tile(specwlen, (len(intensities), 1))
            datafilename = os.path.join(_RECORDINGS_DIR, "V2504041 Xenon 10 Hz lowpass.csv") # TO DO: add example where there's an eye filter
            dataidx = -1

        elif stimulustype == "spectrum matched to 440":

            intensities_log = [6.67, 7.00, 7.29, 7.61, 7.93,8.23, 8.60, 8.93, 9.22, 9.49,
                               6.67, 7.00, 7.29, 7.61, 7.93, 8.23, 8.60, 8.93, 9.22, 9.49]

            dwlen = np.mean(np.diff(specwlen))
            photon_count = np.trapezoid(specint, specwlen, dwlen)
            norm_specint = specint / photon_count
            norm_intensity_440 = np.zeros(len(specwlen))
            norm_intensity_440[np.argmin(np.abs(specwlen - 440))] = 1 / dwlen
            scale_spec_to440 = scalingConstant(specwlen, norm_intensity_440, specwlen, norm_specint)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 440 nm photons: ', scale_spec_to440)
            intensities440 = [ # Retinal irradiance for voltage-clamp 440-440 stimuli
                                  10**6.67, 10**7.00, 10**7.29, 10**7.61, 10**7.93,
                                  10**8.23, 10**8.60, 10**8.93, 10**9.22, 10**9.49,
                                  0,
                                  10**6.67, 10**7.00, 10**7.29, 10**7.61, 10**7.93,
                                  10**8.23, 10**8.60, 10**8.93, 10**9.22, 10**9.49,
                                  0,
                              ]
            intensities = np.array([scale_spec_to440 * norm_specint * intensity for intensity in intensities440])
            wlen = np.tile(specwlen, (len(intensities), 1))
            dataidx = -1

        elif stimulustype == "spectrum matched to 560":

            intensities_log = [6.82, 7.15, 7.45, 7.76, 8.09, 8.38, 8.76, 9.08, 9.38, 9.64, 6.82, 7.15, 7.45, 7.76, 8.09, 8.38, 8.76, 9.08, 9.38, 9.64]

            dwlen = np.mean(np.diff(specwlen))
            photon_count = np.trapezoid(specint, specwlen, dwlen)
            norm_specint = specint / photon_count
            norm_intensity_560 = np.zeros(len(specwlen))
            norm_intensity_560[np.argmin(np.abs(specwlen - 560))] = 1 / dwlen
            scale_spec_to560 = scalingConstant(specwlen, norm_intensity_560, specwlen, norm_specint)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 560 nm photons: ', scale_spec_to560)
            intensities560 = [10**6.82, 10**7.15, 10**7.45, 10**7.76, 10**8.09,
                              10**8.38, 10**8.76, 10**9.08, 10**9.38, 10**9.64,
                              0,
                              10**6.82, 10**7.15, 10**7.45, 10**7.76, 10**8.09,
                              10**8.38, 10**8.76, 10**9.08, 10**9.38, 10**9.64,
                              0,] # Retinal irradiance for voltage-clamp 560-560 stimuli
            intensities = np.array([scale_spec_to560 * norm_specint * intensity for intensity in intensities560])
            wlen = np.tile(specwlen, (len(intensities), 1))

            dataidx = -1

        elif stimulustype == "spectrum matched to 440 560":
            # We are going to have to use a double specwlen for this, one for the 440 spectrum and one for the 560 spectrum
            # IMPORTANT TO LOAD IN SPEC WLEN AS specwlen=[wlen440, wlen560] and similarly for the intensities

            # intensities440 = [10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0, 10**6.7, 10**7.0, 10**7.3, 10**7.6, 10**7.9, 10**8.2, 10**8.6, 10**8.9, 10**9.2, 10**9.5, 0]
            intensities_log = [6.67, 7.00, 7.29, 7.61, 7.93, 8.23, 8.60, 8.93, 9.22, 9.49, 6.67, 7.00, 7.29, 7.61, 7.93, 8.23, 8.60,  8.93, 9.22, 9.49]

            # Matching the 440 spectrum
            dwlen440 = np.mean(np.diff(specwlen[0]))
            photon_count440 = np.trapezoid(specint[0], specwlen[0], dwlen440)
            norm_specint440 = specint[0] / photon_count440
            norm_intensity_440 = np.zeros(len(specwlen[0]))
            norm_intensity_440[np.argmin(np.abs(specwlen[0] - 440))] = 1 / dwlen440
            scale_spec_to440 = scalingConstant(specwlen[0], norm_intensity_440, specwlen[0], norm_specint440)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 440 nm photons: ', scale_spec_to440)
            intensities440 = [10**6.67, 10**7.00, 10**7.29, 10**7.61, 10**7.93, 10**8.23, 10**8.60, 10**8.93, 10**9.22, 10**9.49, 0]
            intensities1 = np.array([scale_spec_to440 * norm_specint440 * intensity for intensity in intensities440])
            wlen1 = np.tile(specwlen[0], (len(intensities1), 1))

            # Matching the 560 spectrum
            dwlen560 = np.mean(np.diff(specwlen[1]))
            photon_count560 = np.trapezoid(specint[1], specwlen[1], dwlen560)
            norm_specint560 = specint[1] / photon_count560
            norm_intensity_560 = np.zeros(len(specwlen[1]))
            norm_intensity_560[np.argmin(np.abs(specwlen[1] - 560))] = 1 / dwlen560
            scale_spec_to560 = scalingConstant(specwlen[1], norm_intensity_560, specwlen[1], norm_specint560)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 560 nm photons: ', scale_spec_to560)
            intensities560 = [10**6.67, 10**7.00, 10**7.29, 10**7.61, 10**7.93, 10**8.23, 10**8.60, 10**8.93, 10**9.22, 10**9.49, 0] # Retinal irradiance for voltage-clamp 440-440 stimuli
            intensities2 = np.array([scale_spec_to560 * norm_specint560 * intensity for intensity in intensities560])
            wlen2 = np.tile(specwlen[1], (len(intensities2), 1))

            # plt.plot(wlen2[0], intensities2[0] / np.max(intensities2[0]), color='k')
            # plt.axvline(x=560, color='r')
            # plt.xlabel('Wavelength (nm)')
            # plt.ylabel('Normalized Intensity')
            # plt.show()

            # plt.plot(wlen[12], intensities[12] / np.max(intensities[12]), color='k')
            # plt.xlabel('Wavelength (nm)')
            # plt.ylabel('Normalized Intensity')
            # plt.show()

            intensities = np.append(intensities1, intensities2, axis=0)
            wlen = np.append(wlen1, wlen2, axis=0)
            dataidx = -1

        elif stimulustype == "spectrum matched to 560 440":
            # We are going to have to use a double specwlen for this, one for the 440 spectrum and one for the 560 spectrum
            # IMPORTANT TO LOAD IN SPEC WLEN AS specwlen=[wlen560, wlen440] and similarly for the intensities
            intensities_log = [6.82, 7.15, 7.45, 7.76, 8.09, 8.38, 8.76, 9.08, 9.38, 9.64, 6.71, 7.03, 7.33, 7.64, 7.97, 8.26, 8.64, 8.96, 9.26, 9.53]

            # Matching the 560 spectrum
            dwlen560 = np.mean(np.diff(specwlen[0]))
            photon_count560 = np.trapezoid(specint[0], specwlen[0], dwlen560)
            norm_specint560 = specint[0] / photon_count560
            norm_intensity_560 = np.zeros(len(specwlen[0]))
            norm_intensity_560[np.argmin(np.abs(specwlen[0] - 560))] = 1 / dwlen560
            scale_spec_to560 = scalingConstant(specwlen[0], norm_intensity_560, specwlen[0], norm_specint560)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 560 nm photons: ', scale_spec_to560)
            intensities560 = [10**6.82, 10**7.15, 10**7.45, 10**7.76, 10**8.09, 10**8.38, 10**8.76, 10**9.08, 10**9.38, 10**9.64, 0] # Retinal irradiance for voltage-clamp 560-440 stimuli
            intensities1 = np.array([scale_spec_to560 * norm_specint560 * intensity for intensity in intensities560])
            wlen2 = np.tile(specwlen[0], (len(intensities1), 1))

            # Matching the 440 spectrum
            dwlen440 = np.mean(np.diff(specwlen[1]))
            photon_count440 = np.trapezoid(specint[1], specwlen[1], dwlen440)
            norm_specint440 = specint[1] / photon_count440
            norm_intensity_440 = np.zeros(len(specwlen[1]))
            norm_intensity_440[np.argmin(np.abs(specwlen[1] - 440))] = 1 / dwlen440
            scale_spec_to440 = scalingConstant(specwlen[1], norm_intensity_440, specwlen[1], norm_specint440)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 440 nm photons: ', scale_spec_to440)
            intensities440 = [10**6.71, 10**7.03, 10**7.33, 10**7.64, 10**7.97, 10**8.26, 10**8.64, 10**8.96, 10**9.26, 10**9.53, 0] # Retinal irradiance for voltage-clamp 560-440 stimuli
            intensities2 = np.array([scale_spec_to440 * norm_specint440 * intensity for intensity in intensities440])
            wlen1 = np.tile(specwlen[1], (len(intensities2), 1))

            intensities = np.append(intensities1, intensities2, axis=0)
            wlen = np.append(wlen1, wlen2, axis=0)

            dataidx = -1

        elif stimulustype == "spectrum matched to 560 test":

            dwlen = np.mean(np.diff(specwlen))
            photon_count = np.trapezoid(specint, specwlen, dwlen)
            norm_specint = specint / photon_count
            norm_intensity_560 = np.zeros(len(specwlen))
            norm_intensity_560[np.argmin(np.abs(specwlen - 560))] = 1 / dwlen
            scale_spec_to560 = scalingConstant(specwlen, norm_intensity_560, specwlen, norm_specint)
            print('Scale factor converting this spectrum to melanopsin (R state) equivalent 560 nm photons: ', scale_spec_to560)
            intensities560 = [2*10**9] # Retinal irradiance for voltage-clamp 440-440 stimuli
            intensities = np.array([scale_spec_to560 * norm_specint * intensity for intensity in intensities560])
            wlen = np.tile(specwlen, (len(intensities), 1))

            plt.plot(wlen[0], intensities[0] / np.max(intensities[0]), color='k')
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Normalized Intensity')
            plt.show()

            dataidx = -1

        # # Parameters shared across stimuli...
        stim_start = 11.0938
        stim_duration = 1
        interstim = 140 # Set this to 140 for voltage clamp, 70 for current clamp
        if stimulustype in ("xenon", "xenon eye"):
            # Old layout had placeholder pulses at indices 10 and 21; dark is implicit between intervals.
            timings = _pulse_starts_skipping_old_indices(stim_start, interstim, 22, (10, 21))
        elif stimulustype in ("xenon 560 suppression xenon", "xenon 560 suppression xenon eye"):
            timings = _pulse_starts_skipping_old_indices(stim_start, interstim, 24, (10, 21))
        elif stimulustype == "xenon 560 suppression xenon first set only":
            timings = _pulse_starts_skipping_old_indices(stim_start, interstim, 13, (10,))
        else:
            timings = np.linspace(stim_start, (len(intensities) - 1)*interstim + stim_start, len(intensities))
        timings_off = timings + stim_duration

        #...except for special cases like xenon and suppression pulses
        if (
            stimulustype == "xenon 560 suppression xenon"
            or stimulustype == "xenon 560 suppression xenon first set only"
            or stimulustype == "xenon 560 suppression xenon max sup max"
            or stimulustype == "xenon 560 suppression xenon triple pre sup bright"
            or stimulustype == "xenon 560 suppression xenon eye"
            or stimulustype == "xenon"
            or stimulustype == "xenon eye"
        ):  # Xenon

            if stimulustype in (
                "xenon 560 suppression xenon",
                "xenon 560 suppression xenon first set only",
                "xenon 560 suppression xenon max sup max",
                "xenon 560 suppression xenon triple pre sup bright",
                "xenon 560 suppression xenon eye",
            ):

                timings[-2] = timings[-3] + 0.5*interstim  # Second to last stimulus starts a little later
                timings_off[-2] = timings[-2] + interstim # ...and lasts one interstim
                timings[-1] = timings_off[-2] + interstim # last stimulus starts one interstim later...
                timings_off[-1] = timings[-1] + stim_duration # ...and lasts the usual duration


        timings = [(timings[n],timings_off[n]) for n in range(len(timings))]
        timings = np.array(timings).reshape(-1)

    # This experiment is accessible in the model by using a blank string for experimentype = ""
    elif stimulustype == "560 conditioning":

        # Alternate between 440 and 480 for 10 pairs, then add 560 and 480, then alternate 440/480 4 more times
        pre_pairs = 10
        post_pairs = 40
        wlen = np.array([440, 480] * pre_pairs + [560] + [440, 480] * post_pairs)

        # 10 values of 10**7, then 2*10**9, then 10 more of 10**7
        intensities = np.array([10**7] * pre_pairs*2 + [2 * 10**9] + [10**7] * post_pairs*2)

        intensities_log = [np.log10(intensity) for intensity in intensities]


        # Total number of stimuli: 2*pre_pairs + 1 (560) + 2*post_pairs
        n_stimuli = 2 * pre_pairs + 1 + 2 * post_pairs

        # Start times spaced 70 s apart starting at 10 s
        start_times = list(range(10, 10 + 70 * n_stimuli, 70))

        durations = 0.05 * np.ones(np.size(start_times))
        durations[wlen == 560] = 30
        end_times = [start_time + duration for start_time, duration in zip(start_times, durations)]
        timings = [val for pair in zip(start_times, end_times) for val in pair]
        dataidx = -1
        interstim = timings[-1] - timings[-3]

    # This experiment is accessible in the model by using a blank string for experimentype = ""
    elif stimulustype == "560 sensitivity recovery":

        # pre and post specify number of probes
        pre = 10
        post = 40
        wlen = np.array([520] * pre + [560] + [520] * post)

        # sensitivity probes, then conditioning, then more probes
        intensities = np.array([10**4] * pre + [3 * 10**9] + [10**7] * post)


        # Total number of stimuli:
        n_stimuli = pre + 1 + post

        # Start times spaced 70 s apart starting at 10 s
        start_times = list(range(10, 10 + 70 * n_stimuli, 70))

        durations = 0.05 * np.ones(np.size(start_times))
        durations[wlen == 560] = 60
        end_times = [start_time + duration for start_time, duration in zip(start_times, durations)]
        timings = [val for pair in zip(start_times, end_times) for val in pair]
        dataidx = -1
        interstim = timings[-1] - timings[-3]

    # Background light to probe equilibrium
    elif stimulustype == "600 nm background":

        wlen = [600]
        intensities = [7*10**8]
        timings = [10,6000]
        interstim = 60

    elif stimulustype == "flashes on 600 nm background":

        # Define wavelength list
        wls = [410.0, 440.0, 480.0, 500.0, 510.0, 520.0, 560.0, 600.0, 640.0, 680.0]

        # Define background and flash intensities
        dark = np.zeros(len(wls))
        #backgroundint = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0*10**8, 0, 0]) # 600 nm background - upper range from Emanuel and Do 2015
        #backgroundint = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.3*10**7, 0, 0]) # 600 nm background - middle range from Emanuel and Do 2015
        backgroundint = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0*10**6,  0, 0]) # 600 nm background - lower range from Emanuel and Do 2015

        #intensities_log = [4.8, 4.6, 4.6, 4.8, 4.8, 5.0, 6.0, 7.2, 8.2, 8.8]
        intensities_log = [4.8, 4.6, 4.6, 4.8, 4.8, 5.0, 7.2, 8.2, 9.2, 9.8] # Govardovskii E

        wavelengths = ["410", "440", "480", "500", "515", "520", "560", "600", "640", "680"]

        flashes = {}

        for i, wl in enumerate(wavelengths):
            row = [0] * len(intensities_log)
            row[i] = 10 ** intensities_log[i]
            flashes[wl] = row


        # Build intensity list alternating dark and flash for each wavelength
        intensities = []
        for key in flashes: # flashes on darkness
            intensities.append(dark)
            intensities.append(np.array(flashes[key]))
        intensities.append(dark) # Rest before background
        intensities.append(backgroundint) # Turn on background
        for key in flashes: # flashes on background
            intensities.append(dark + backgroundint)
            intensities.append(np.array(flashes[key]) + backgroundint) # Make brighter on the background (*10 for 600 nm background, *10**3 for 440 nm)
        intensities.append(dark + backgroundint)  # Final rest period

        # Build wlen: repeats of the same wls for every stimulus intensity
        wlen = [wls] * len(intensities)

        # Generate timings
        start_time = 10
        steady_state_time = 6000 # time to reach steady state after background turns on. 600 is sufficient for brighter background (> ~6 log ph)
        pulse_duration = 0.05
        interstim = 360 # was 90
        timings = [0]

        for i in range(len(wls)):  # flashes on darkness
            stim_on = start_time + i * interstim
            stim_off = stim_on + pulse_duration
            timings.extend([stim_on, stim_on, stim_off, stim_off])

        timings.extend([stim_off+interstim, stim_off+interstim]) # Rest in dark to allow last flash to decay (stim_off) then turn on background (stim_off+interstim)

        background_ss = stim_off + interstim + steady_state_time

        timings.extend([background_ss, background_ss]) # keep background on long enough to reach steady state before delivering flashes

        for i in range(len(wls)):  # flashes on background
            stim_on = background_ss + i * interstim
            stim_off = stim_on + pulse_duration
            timings.extend([stim_on, stim_on, stim_off, stim_off])

        timings.append(stim_off + interstim)  # Final rest period

    elif stimulustype == "check equilibrium":

        intensities_log = [np.log10(10**7.5)]

        dwlen = np.mean(np.diff(specwlen))
        photon_count = np.trapezoid(specint, specwlen, dwlen)
        norm_specint = specint / photon_count
        norm_intensity_440 = np.zeros(len(specwlen))
        norm_intensity_440[np.argmin(np.abs(specwlen - 440))] = 1 / dwlen
        scale_spec_to440 = scalingConstant(specwlen, norm_intensity_440, specwlen, norm_specint)
        print('Scale factor converting this spectrum to melanopsin (R state) equivalent 440 nm photons: ', scale_spec_to440)
        intensities440 = [10**7.5] # Retinal irradiance for voltage-clamp 440-440 stimuli
        intensities = np.array([scale_spec_to440 * norm_specint * intensity for intensity in intensities440])
        wlen = np.tile(specwlen, (len(intensities), 1))
        timings = [1,600]
        interstim = 1
        dataidx = -1


    else: # Enter custom stimuli here.

        wlen = [440, 440, 440, 440, 440, 440, 560, 440, 440, 440, 440, 440, 440, 440]
        intensities = [10**10]*13 + [0]
        timings = [10,11,  70,71,  130,131,  190,191, 250,251, 310,311, 370,400,
                460,461, 520,521, 580,581, 640,641, 700,701, 760,761, 820,831]
        interstim = timings[-1] - timings[-3]

    # Set the simulation time
    ti = np.max(timings) + interstim

    # Throw all the model parameters into a dictionary
    outDict = {}
    #outDict['stim_start'] = stim_start
    #outDict['stim_duration'] = stim_duration
    outDict['timings'] = np.array(timings)
    outDict['interstim'] = np.array(interstim)
    outDict['ti'] = ti*dilation # Lower the simulation time (e.g. ti / 10) for lower runtime during debugging
    outDict['intensities'] = np.array(intensities)
    outDict['datafilename'] = datafilename
    outDict['irdatafilenamet'] = irdatafilenamet
    outDict['irdatafilenamep'] = irdatafilenamep
    outDict['dataidx'] = dataidx
    outDict['wlen'] = np.array(wlen)
    # outDict['intensitieslog'] = intensities_log # Log10 photon counts, for plotting intensity-response relations
    if intensities_log != []:
        outDict['intensitieslog'] = intensities_log # do this to get around the intensities_log issue if not defined for a specific stimulus
    if stimulustype == "spectrum matched to 560 test":
        outDict['ti'] = ti*5

    return outDict


def get_predict_melanopsin_defaults():
    """
    Default parameter dict for ``predictMelanopsin`` (merged with the ``config`` argument).

    Exposed so UIs and scripts can read the same keys and defaults as the model without
    duplicating the table.
    """
    return {

        # Simulation time step, in seconds
        'rate': 0.01, # Euler blow-up at rates coarser than 0.02 sec. 

        # Toggle features
        'activity feedback': True,
        'Mprime feedback': True,

        # Quantum efficiencies
        'Rphi': 0.52,  # measured by Shichida
        'Mphi': 0.22,  # measured by Shichida
        'Ephi': 0.65,  # Not measured by Shichida. Initially guessed an intermediate value of 0.4 but >= 0.6 appears to give more realistic fit, especially to Xenon data

        # Absorption coefficients in units of cm^3/(mol*cm)
        'Reps': 33000 * 1000,
        'Meps': 52600 * 1000,
        'Eeps': 42000 * 1000,

        # Choose nomogram function
        'spectral sensitivity': 'Custom', # 'Govardovskii' or 'Custom'. If Govardvoskii, only the 'true_lambda_max' parameters will be used.

        # R state nomogram
        'R_true_lambda_max': 467,
        'R_lambda_max_param': 494,
        'R_A': 44.582569987408945,
        'R_B': 8.614863087478009,
        'R_C': -11.594413782436712,
        'R_D': 0.0,
        'R_a': 0.9324464256209672,
        'R_b': 0.9901137601105081,
        'R_c': 1.128289328759767,

        # M state nomogram
        'M_true_lambda_max': 476,
        'M_lambda_max_param': 483.0,
        'M_A': 44.402,
        'M_B': 11.768,
        'M_C': -14.088,
        'M_D': 0.2554,
        'M_a': 0.90595,
        'M_b': 0.87117,
        'M_c': 1.1367,

        # E state nomogram
        'E_true_lambda_max': 446,
        'E_lambda_max_param': 446,
        'E_A': 69.7,
        'E_B': 28,
        'E_C': -14.9,
        'E_D': 0.674,
        'E_a': 0.8795 + 0.0459 * np.exp(-(446-300)**2 / 11940),
        'E_b': 0.922,
        'E_c': 1.104,

        # Transition probability of M to E (M to R will be set to 1-[P(M to E)])
        'propME': 0.5,

        ## Rate constants (transitions * s-1) for...
        # ...transitions between high and low gain states
        'konM': 0.012,        # M to M'
        'koffM': 0.0012,      # M' to M
        'konR': 1/3600,       # R to R'
        'koffR': 0.0454545,   # R' to R
        'konE': 1/3600,       # E to E'
        'koffE': 0.0454545,   # E' to E

        # ...leaky integrators for gain control
        'konMG': 1.0,           # onset of M'-driven negative feedback
        'koffMG': 0.355,        # decay of M'-driven negative feedback
        'konGG': 1.0,           # onset of activity-driven negative feedback
        'koffGG': 0.037,        # decay of activity-driven negative feedback

        # ...bleaching
        'kbleachR': 0.75e-2,
        'kbleachRprime': 0.75e-2,
        'kbleachM': 1e-1,
        'kbleachMprime': 1e-4,
        'kbleachE': 5e-4,
        'kbleachEprime': 1e-4,
        'koffO': 0.75e-2,

        # ...thermal transitions of M to ground states
        'kthermM': 0.0, # off by default
        'propMEtherm': 0.5, # given that M has a thermal transtion, propbability it transitions to E

        # Threshold for faster M to M' conversion
        'm deplete thresh': 0.05,
        'm deplete scalar': 2.0,

        # Signaling gain
        'MprimeGain': 0.2,  # M prime
        'Ogain': 10**-6.5,  # Bleached state. Based on Carter Cornwall measurements of salamander rods

        # Parameters for nonlinearity (Weibull function) shaping M'-driven negative feedback
        'k_mg': 1.215,
        'x0_mg': 0.381,

        # Sigmoidal nonlinearity applied to activity-driven negative feedback
        'k': 40.0,
        'L': 0.99999999999999,
        'x0': 0.0008,

    }


class SimulationCancelled(Exception):
    """Raised when ``predictMelanopsin`` exits early because ``cancel_event`` was set."""


def predictMelanopsin(
    stimulus,
    config=None,
    progress_cb=None,
    progress_every=100,
    cancel_event=None,
):
    """
    Description
    -----------
    Simulate opsin equilibrium behavior when given a light stimulus.
    Core function of the model.

    Inputs
    ----------
    stimulus : Python Dictionary
        Input dictionary.
    config: Python Dictionary
        User defined model parameters not accounted for in 'stimulus'.
        Must be model-specific parameters, not absolute custom parameters.
        All specific keys can be found in the appendix of the tutorial notebook.
    progress_cb : callable, optional
        Optional callback for reporting simulation progress. Called as
        ``progress_cb(current_step, total_steps)``.
    progress_every : int
        Callback update interval in loop steps. Lower values update more frequently.
    cancel_event : threading.Event, optional
        If set, the main integration loop raises ``SimulationCancelled`` at the
        next timestep (cooperative cancellation for threaded UIs).

    Returns
    -------
    outDict
        A dictionary containing the following keys:
        xaxis : ndarray
            Time axis in s
        R : ndarray
            Pigment fraction in R state.
        M : ndarray
            Pigment fraction in M state.
        E : ndarray
            Pigment fraction in E state.
        O : ndarray
            Pigment fraction in O state (unavailable e.g., bleached, internalized, etc.).
        Rprime : ndarray
            Pigment fraction in R' state.
        Mprime : ndarray
            Pigment fraction in M' state.
        Eprime : ndarray
            Pigment fraction in E' state.
        gg : float
            global gain
        mg : float
            M gain
        currentOut : ndarray
            1.0 * M + MprimeGain * M'
        currentNoLowGain : ndarray
            1.0 * M + 1.0 * M'
        currentGlobalGain : ndarray
            M and MprimeGains are adjusted by recent current.
        lightIntensity : ndarray
            intensity of light at each sample.
        lightWavelength : ndarray
            The range of wavelengths the model was run over
        R2M, Rprime2Mprime, E2M, Eprime2Mprime : ndarray
            Per-step photoconversion fluxes into M / M' from R, R', E, E' (when light is on).
        lightPulseID : ndarray
            Index of the active stimulus interval per time step (-1 = dark).
    """

    intensity = np.array(stimulus['intensities'])
    intensityChecker(intensity)

    wlen = np.array(stimulus['wlen'])
    wlenChecker(wlen)

    lightTimings = np.array(stimulus['timings'])
    if all(isinstance(i, tuple) for i in lightTimings) == True:
        lightTimings = lightTimings.flatten()
    lightTimingsCheck(lightTimings)

    intWlenLtChecker(intensity, wlen, lightTimings)

    ti = stimulus['ti'] # duration of the simulation (typically end of last light stim + interstimulus interval)
    tiChecker(ti)

    if np.max(lightTimings) > ti:
        raise ValueError("'ti' is smaller than max(lightTimings). 'ti' must be large enough to accomodate for all desired stimulus intervals.")

    # delay the lightTimings by 400 msec to approximate response latency in ipRGCs
    lightTimings = lightTimings + 0.4

    defaults = get_predict_melanopsin_defaults()

    # Merge user-supplied parameters (config) with defaults
    if config is None:
        config = {}
    unknown_keys = set(config.keys()) - set(defaults.keys())
    if unknown_keys:
        warnings.warn(
            f"Config contains key(s) the model does not use (typo or invalid parameter); they have no effect: {sorted(unknown_keys)}",
            UserWarning,
            stacklevel=2,
        )
    params = {**defaults, **config}

    # Adjust rate constants to the simulation rate
    rate_constant_keys = ['konM', 'koffM', 'konR', 'koffR', 'konE', 'koffE',
                          'konMG', 'koffMG',
                          'konGG', 'koffGG',
                          'kbleachR', 'kbleachRprime', 'kbleachM', 'kbleachMprime', 'kbleachE', 'kbleachEprime', 'koffO',
                          'kthermM']
    for key in rate_constant_keys:
        params[key] = params[key] * params['rate']

    # Unpack parameters into variables
    (
    activity_feedback, Mprime_feedback,
    rate, Rphi, Mphi, Ephi, Reps, Meps, Eeps,
    spectral_sensitivity,
    propME,
    konM, koffM, konR, koffR, konE, koffE,
    konMG, koffMG, konGG, koffGG,
    kbleachR, kbleachRprime, kbleachM, kbleachMprime, kbleachE, kbleachEprime,
    koffO, Ogain,
    MprimeGain,
    k_mg, x0_mg,
    k, L, x0,
    mdeplete_thresh, mdeplete_scalar,
    propMEtherm, kthermM

    ) = [params[k] for k in [

    'activity feedback', 'Mprime feedback',
    'rate', 'Rphi', 'Mphi', 'Ephi', 'Reps', 'Meps', 'Eeps',
    'spectral sensitivity',
    'propME',
    'konM', 'koffM', 'konR', 'koffR', 'konE', 'koffE',
    'konMG', 'koffMG', 'konGG', 'koffGG',
    'kbleachR', 'kbleachRprime', 'kbleachM', 'kbleachMprime', 'kbleachE', 'kbleachEprime',
    'koffO', 'Ogain',
    'MprimeGain',
    'k_mg', 'x0_mg',
    'k', 'L', 'x0',
    'm deplete thresh', 'm deplete scalar',
    'propMEtherm', 'kthermM'

    ]]

    # transition probability of M to R
    if propME > 1 or propME < 0:
        raise ValueError()
    propMR = 1 - propME

    # check inputs for validity
    # assert int(len(lightTimings) / len(intensity)) == 2, 'Numbers of pulses don\'t match between intensity & lightTimings'
    # assert len(wlen) / len(intensity) == 2, 'Numbers of pulses don\'t match between wlen & lightTimings'

    # initial global gain and current
    currentGlobalGain = np.zeros(int(ti/rate))
    gg = np.zeros(int(ti/rate))
    mg = np.zeros(int(ti/rate))
    gg_scalar = np.zeros(int(ti/rate))
    mg_scalar = np.zeros(int(ti/rate))

    #gg = global gain, mg = M gain (ends up getting set by M')
    gg[0] = 0
    mg[0] = 0
    currentGlobalGain[0] = 0

    # pigment state arrays
    R = np.zeros(int(ti/rate))
    M = np.zeros(int(ti/rate))
    E = np.zeros(int(ti/rate))
    O = np.zeros(int(ti/rate))
    Rprime = np.zeros(int(ti/rate))
    Mprime = np.zeros(int(ti/rate))
    Eprime = np.zeros(int(ti/rate))

    # initialize photoconversion counters
    R2M = np.zeros(int(ti/rate))
    Rprime2Mprime = np.zeros(int(ti/rate))
    E2M = np.zeros(int(ti/rate))
    Eprime2Mprime = np.zeros(int(ti/rate))
    M2E = np.zeros(int(ti/rate))
    Mprime2Eprime = np.zeros(int(ti/rate))
    M2R = np.zeros(int(ti/rate))
    Mprime2Rprime = np.zeros(int(ti/rate))

    # Dark start: R↔O equilibrium. With default kbleachR = koffO, half of melanopsin
    # is in O; remaining chromophore-bound pigment is in R (M, E, and primed states = 0).
    R_init =  1 / (1 + kbleachR / koffO)
    R[0] = R_init
    M[0] = 0
    E[0] = 0
    O[0] = 1 - R_init
    Rprime[0] = 0
    Mprime[0] = 0
    Eprime[0] = 0
    M_new = np.zeros(int(ti/rate))

    # Handle 1D or 2D intensity array
    if intensity.ndim == 1: # 1D intensity array means the stimulus is defined in monochromatic mode

        # Reshape to a 2D array with one column if it's 1D
        lightIntensity = np.zeros([len(R), 1])
        lightWavelength = np.zeros([len(R), 1])

    else:

        lightIntensity = np.zeros([len(R),len(intensity)]) 
        lightWavelength = np.zeros([len(R),len(intensity)])

    lightPulseID = np.ones(len(R),dtype=np.int32)*-1

    # Extract stimulus intervals (defined by pairs of start,stop times)
    for i in range(0,len(lightTimings),2):
        if intensity.ndim == 1: # monochromatic mode
            lightIntensity[int(lightTimings[i]/rate):int(lightTimings[i+1]/rate)] = intensity[int(i/2)] 
            lightWavelength[int(lightTimings[i]/rate):int(lightTimings[i+1]/rate)] = wlen[int(i/2)] 
        lightPulseID[int(lightTimings[i]/rate):int(lightTimings[i+1]/rate)] = int(i/2)

    # Apply nomogram function (determine spectral sensitivity to all input wavelengths) 
    if spectral_sensitivity == 'Govardovskii':

        nomoR = [gov(n, 467) for n in wlen]
        nomoM = [gov(n, 476) for n in wlen]
        nomoE = [gov(n, 446) for n in wlen]

    elif spectral_sensitivity == 'Custom':

        nomoR = [customNomogram(n, getNomogramParams(params, 'R')) for n in wlen]
        nomoM = [customNomogram(n, getNomogramParams(params, 'M')) for n in wlen]
        nomoE = [customNomogram(n, getNomogramParams(params, 'E')) for n in wlen]

    # Convert intensity to moles per centimeter squared per time slice per nm
    if wlen.ndim > 1:
        scaleFactor = (wlen[0][1] - wlen[0][0]) ## nm per slice
        intensity_conv = [(n/6.022e23)*1e8*scaleFactor for n in intensity] # mol per cm^2 per s per nm. Integrate over the wavelength bin using scaleFactor.
    else:
        intensity_conv = [(n/6.022e23)*1e8 for n in intensity] # mol per cm^2 per s

    intensity_conv = [n*rate for n in intensity_conv] # now mol per cm^2 per time slice
    

    
    # Probability of photon absorption
    crsum = [np.sum(R_[0] * R_[1] * Reps) for R_ in zip(nomoR,intensity_conv)]
    cmsum = [np.sum(M_[0] * M_[1] * Meps) for M_ in zip(nomoM,intensity_conv)]
    cesum = [np.sum(E_[0] * E_[1] * Eeps) for E_ in zip(nomoE,intensity_conv)]

    # Run the simulation
    total_steps = int(ti / rate) - 1
    if progress_cb is not None:
        progress_cb(0, total_steps)
    for i in tqdm(range(total_steps), desc = "The Main Loop"):
        if cancel_event is not None and cancel_event.is_set():
            raise SimulationCancelled

        ## Three stages -- first consider photoconversions
        if lightPulseID[i] != -1 :

            R2M[i+1] = np.log(10)*Rphi*crsum[lightPulseID[i]]*R[i]
            E2M[i+1] = np.log(10)*Ephi*cesum[lightPulseID[i]]*E[i]
            M2E[i+1] = np.log(10)*Mphi*cmsum[lightPulseID[i]]*M[i] * propME
            M2R[i+1] = np.log(10)*Mphi*cmsum[lightPulseID[i]]*M[i] * propMR

            Rprime2Mprime[i+1] = np.log(10)*Rphi*crsum[lightPulseID[i]]*Rprime[i]
            Eprime2Mprime[i+1] = np.log(10)*Ephi*cesum[lightPulseID[i]]*Eprime[i]
            Mprime2Eprime[i+1] = np.log(10)*Mphi*cmsum[lightPulseID[i]]*Mprime[i] * propME
            Mprime2Rprime[i+1] = np.log(10)*Mphi*cmsum[lightPulseID[i]]*Mprime[i] * propMR

            R[i+1] = R[i] - R2M[i+1] + M2R[i+1]
            M[i+1] = M[i] - M2E[i+1] - M2R[i+1] + E2M[i+1] + R2M[i+1]
            E[i+1] = E[i] - E2M[i+1] + M2E[i+1]

            Rprime[i+1] = Rprime[i] - Rprime2Mprime[i+1] + Mprime2Rprime[i+1]
            Mprime[i+1] = Mprime[i] - Mprime2Eprime[i+1] - Mprime2Rprime[i+1] + Rprime2Mprime[i+1] + Eprime2Mprime[i+1]
            Eprime[i+1] = Eprime[i] - Eprime2Mprime[i+1] + Mprime2Eprime[i+1]
        else:
            M_new[i+1] = 0
            R[i+1] = R[i]
            Rprime[i+1] = Rprime[i]
            M[i+1] = M[i]
            Mprime[i+1] = Mprime[i]
            E[i+1] = E[i]
            Eprime[i+1] = Eprime[i]


        ## Second, consider light-independent gain-state conversions

        # When Mprime fraction exceeds a threshold, begin transitioning M to M' faster
        if Mprime[i+1] > mdeplete_thresh:
            mdeplete = mdeplete_scalar
        else:
            mdeplete = 1.0

        O[i+1] = O[i] + M[i+1] * kbleachM + Mprime[i+1] * kbleachMprime + E[i+1] * kbleachE + Eprime[i+1] * kbleachEprime + R[i+1] * kbleachR + Rprime[i+1] * kbleachRprime - O[i] * koffO

        R[i+1] = R[i+1] + Rprime[i+1] * koffR - R[i+1] * konR + O[i] * koffO - R[i+1] * kbleachR + kthermM * M[i+1] * (1-propMEtherm)
        Rprime[i+1] = Rprime[i+1] + R[i+1] * konR - Rprime[i+1] * koffR - Rprime[i+1] * kbleachRprime

        M[i+1] = M[i+1] + Mprime[i+1] * koffM - M[i+1] * mdeplete * konM - M[i+1] * kbleachM - kthermM * M[i+1]
        Mprime[i+1] = Mprime[i+1] + M[i+1] * mdeplete * konM - Mprime[i+1] * koffM - Mprime[i+1] * kbleachMprime

        E[i+1] = E[i+1] + Eprime[i+1] * koffE - E[i+1] * konE - E[i+1] * kbleachE + kthermM * M[i+1] * propMEtherm
        Eprime[i+1] = Eprime[i+1] + E[i+1] * konE - Eprime[i+1] * koffE - Eprime[i+1] * kbleachEprime

        # Renormalization (counters accumulating numerical errors over long simulations)
        total_pigment = R[i+1] + M[i+1] + E[i+1] + O[i+1] + Rprime[i+1] + Mprime[i+1] + Eprime[i+1]
        R[i+1] /= total_pigment
        M[i+1] /= total_pigment
        E[i+1] /= total_pigment
        O[i+1] /= total_pigment
        Rprime[i+1] /= total_pigment
        Mprime[i+1] /= total_pigment
        Eprime[i+1] /= total_pigment

        # Check to make sure state at this timepoint is between [0,1].
        stateBoundsChecker(R, M, E, O, Rprime, Mprime, Eprime, idx=i+1)

        ## Third, consider two forms of gain control.
        # gain for M (mg): M', passed through a leaky integrator and weibull, reduces the gain of M
        # global gain (gg): M, passed through a distinct leaky integrator and sigmoid, reduces the gain of both M and M'

        # Leaky integration of mg
        minGainUpdate_mg = 0.0

        if Mprime[i] > minGainUpdate_mg:

          mg[i+1] = mg[i] + Mprime[i] * konMG - mg[i] * koffMG

        else:

          mg[i+1] = mg[i] - mg[i] * koffMG

        if mg[i+1] < 0:

          mg[i+1] = 0

        # Leaky integration of gg
        if currentGlobalGain[i] > 0: 

          gg[i+1] = gg[i] + currentGlobalGain[i] * konGG - gg[i] * koffGG # GG based on activity

        else:

          gg[i+1] = gg[i] - gg[i] * koffGG

        if gg[i+1] < 0:

          gg[i+1] = 0


        # Apply nonlinearities
        lowestgain = 1e-10 # was 0.01

        # Apply nonlinearity to integrated M' fraction (persistent activity)
        mg_scalar[i+1] = np.max([weib(mg[i+1], k=k_mg, lam=x0_mg), lowestgain])

        # Apply nonlinearity to integrated activity
        gg_scalar[i+1] = np.max([1 - sigmoid(gg[i+1], k, L, x0), lowestgain])

        # Update current
        if activity_feedback and Mprime_feedback:

          currentGlobalGain[i+1] =  gg_scalar[i+1] * (mg_scalar[i+1] * M[i+1] + Mprime[i+1] * MprimeGain + O[i+1] * Ogain)

        elif activity_feedback and not Mprime_feedback:

          currentGlobalGain[i+1] =  gg_scalar[i+1] * (M[i+1] + Mprime[i+1] * MprimeGain + O[i+1] * Ogain)

        elif Mprime_feedback and not activity_feedback:

          currentGlobalGain[i+1] =  mg_scalar[i+1] * M[i+1] + Mprime[i+1] * MprimeGain + O[i+1] * Ogain

        elif not activity_feedback and not Mprime_feedback:

          currentGlobalGain[i+1] =  M[i+1] + Mprime[i+1] * MprimeGain + O[i+1] * Ogain

        if progress_cb is not None and ((i + 1) % max(1, int(progress_every)) == 0 or i + 1 == total_steps):
            progress_cb(i + 1, total_steps)

    ## Set the current
    currentOut = M * 1.0 + Mprime * MprimeGain + O * Ogain # Low gain states

    currentNoLowGain = M * 1.0 + Mprime * 1.0 # No low gain states

    xaxis = np.linspace(0, ti, num=len(currentGlobalGain)) # Timebase

    # Combine the outputs into a dictionary
    outDict = {}
    outDict['xaxis'] = xaxis
    outDict['R'] = R
    outDict['M'] = M
    outDict['E'] = E
    outDict['O'] = O
    outDict['Rprime'] = Rprime
    outDict['Mprime'] = Mprime
    outDict['Eprime'] = Eprime
    outDict['gg'] = gg
    outDict['mg'] = mg
    outDict['gg_scalar'] = gg_scalar
    outDict['mg_scalar'] = mg_scalar
    outDict['currentOut'] = currentOut
    outDict['currentNoLowGain'] = currentNoLowGain
    outDict['currentGlobalGain'] = currentGlobalGain
    outDict['lightIntensity'] = lightIntensity
    outDict['lightWavelength'] = lightWavelength
    outDict['stimMonitor'] = stimMonitor(stimulus, xaxis)
    outDict['rate'] = rate
    outDict['params'] = params
    # Per–time-step photoconversion fluxes (fraction of pigment pool moved in one step; only
    # nonzero at indices i+1 when light was on at step i). Same definitions for Govardovskii / Custom.
    outDict['R2M'] = R2M
    outDict['Rprime2Mprime'] = Rprime2Mprime
    outDict['E2M'] = E2M
    outDict['Eprime2Mprime'] = Eprime2Mprime
    outDict['lightPulseID'] = lightPulseID

    return outDict


def photoconversion_R_vs_E_during_light(model_out):
    """
    Description
    -----------
    Sum R→M and R'→M' vs E→M and E'→M' photoconversion over time steps when light was on.
    Uses the same pulse mask as the simulation (``lightPulseID``), i.e. after any internal
    timing adjustments in :func:`predictMelanopsin`.

    Parameters
    ----------
    model_out : dict
        Output of :func:`predictMelanopsin` (must include R2M, Rprime2Mprime, E2M, Eprime2Mprime, lightPulseID).

    Returns
    -------
    dict
        Keys: ``from_R``, ``from_E`` (sums of step fluxes), ``total``, ``frac_R``, ``frac_E``.
        If ``total`` is 0, fractions are NaN.
    """
    lp = np.asarray(model_out['lightPulseID'])
    R2M = np.asarray(model_out['R2M'])
    Rp = np.asarray(model_out['Rprime2Mprime'])
    E2M = np.asarray(model_out['E2M'])
    Ep = np.asarray(model_out['Eprime2Mprime'])
    on = lp[:-1] != -1
    from_R = float(np.sum(R2M[1:][on]) + np.sum(Rp[1:][on]))
    from_E = float(np.sum(E2M[1:][on]) + np.sum(Ep[1:][on]))
    total = from_R + from_E
    if total <= 0:
        return {
            'from_R': from_R,
            'from_E': from_E,
            'total': total,
            'frac_R': np.nan,
            'frac_E': np.nan,
        }
    return {
        'from_R': from_R,
        'from_E': from_E,
        'total': total,
        'frac_R': from_R / total,
        'frac_E': from_E / total,
    }


def stimMonitor(setup, timebase):
    intensities = setup['intensities']
    timings = setup['timings']
    wlens = setup['wlen']
    monitor = np.ones_like(timebase)


    for i in range(0, len(timings), 2):
        start = timings[i]
        end = timings[i + 1]
        intensity = intensities[i // 2]
        wlen = wlens[i // 2]

        # If the stimulus is defined as spectrum, integrate the spectral intensity to obtain photon count
        if wlen.ndim > 0:
            dwlen = np.median(np.diff(wlen))
            intensity = np.trapezoid(intensity, wlen, dwlen)

        # Find indices in timebase where time is within the pulse window
        idx = np.where((timebase >= start) & (timebase < end))[0]
        monitor[idx] = intensity

    return monitor


# Re-export for backward compatibility; implementation lives in plotting_functions
from .plotting_functions import plotModelRun