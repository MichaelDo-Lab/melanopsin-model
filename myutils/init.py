# Standard Library Imports
import csv
import glob
import math
import os
import re
import subprocess
import time
import warnings
from datetime import datetime
from functools import partial
from pathlib import Path

# Third-Party Imports
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ipywidgets is only needed for the interactive notebook controls. It pulls in a
# large Jupyter dependency tree, so it is optional: the desktop GUI runs without
# it (the values below are simply ``None`` when it is unavailable).
try:
    import ipywidgets as widgets
    from ipywidgets import interact
except Exception:  # pragma: no cover - notebook-only dependency
    widgets = None
    interact = None

from joblib import Parallel, delayed
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfMerger
from scipy import signal
from scipy.interpolate import (
    Akima1DInterpolator,
    CubicSpline,
    PchipInterpolator,
    interp1d,
)
from scipy.optimize import (
    Bounds,
    basinhopping,
    differential_evolution,
    minimize,
)
from skopt import gp_minimize
# tqdm.auto degrades gracefully to a console bar outside of Jupyter.
from tqdm.auto import tqdm

# Custom library imports
from myutils.app_paths import ensure_dir, manuscript_dir, manuscript_spectra_dir
from myutils.model_functions import powerSpectrumToPhotonFlux

def loadInSpectra():
    """
    Description
    -----------
    Loads in experimental spectral data into the notebook environment upon initialization.

    Parameters
    ----------
        None.

    Returns
    -------
        spec440 : array shape (1,)
            440nm experimental spectrum wavelengths.
        spec560 : array shape (1,)
            560nm experimental spectrum wavelengths.
        intensity440 : array shape (1,)
            440nm experimental spectrum intensities. Same length as 'spec440'
        intensity560 : array shape (1,)
            560nm experimental spectrum intensities. Same length as 'spec560'.
        wlenshared : array shape (1,)
            Common wavelength range used by the model.
        xenon : array shape (2,)
            xenon[0] : array shape (1,)
                Wavelengths of xenon spectrum
            xenon[1] : array shape (1,0)
                Intensities of xenon spectrum. Same length as 'xenon[0]'
        xenon_eye : array shape (2,)
            xenon_eye[0] : array shape (1,)
                Wavelengths of xenon spectrum with eye filter.
            xenon_eye[1] : array shape (1,0)
                Intensities of xenon spectrum with eye filter.. Same length as 'xenon_eye[0]'
    """
    data_path = manuscript_spectra_dir()

    # Load in pulses
    spec440 = pd.read_csv(data_path / '440nm spectrum.csv')
    spec560 = pd.read_csv(data_path / '560nm spectrum.csv')

    # Interpolate to a uniform wavelength base
    f440 = interp1d(spec440["wavelength"].to_numpy(), spec440["power"].to_numpy(), kind='linear', bounds_error=False, fill_value='extrapolate')
    f560 = interp1d(spec560["wavelength"].to_numpy(), spec560["power"].to_numpy(), kind='linear', bounds_error=False, fill_value='extrapolate')

    wlenshared = np.linspace(200,800,num=3648)
    spec440["power"] = f440(wlenshared)
    spec560["power"] = f560(wlenshared)

    # Convert power spectra to intensities
    intensity440 = powerSpectrumToPhotonFlux(spec440["power"], wlenshared)
    intensity440 = np.array(intensity440)

    intensity560 = powerSpectrumToPhotonFlux(spec560["power"], wlenshared)
    intensity560 = np.array(intensity560)

    # Trim spectra
    intensity440[wlenshared < 435] = 0.001
    intensity440[wlenshared > 448] = 0.001

    intensity560[wlenshared < 549] = 0.001
    intensity560[wlenshared > 570] = 0.001

    # Load in Xenon and xenon with an eye filter applied
    xenon = pd.read_table(os.path.join(data_path, 'Xenon.txt'),  names=['wavelength', 'power'])
    xenon_eye = pd.read_table(os.path.join(data_path, 'XeEye.txt'),  names=['wavelength', 'power'])

    fxenon = interp1d(xenon['wavelength'].to_numpy(), xenon['power'].to_numpy(), kind='linear', bounds_error=False, fill_value='extrapolate')
    fxenon_eye = interp1d(xenon_eye['wavelength'].to_numpy(), xenon_eye['power'].to_numpy(), kind='linear', bounds_error=False, fill_value='extrapolate')

    xenon['power'] = fxenon(wlenshared)
    xenon_eye['power'] = fxenon_eye(wlenshared)

    # Convert power to photon flux
    xenon['photons'] = powerSpectrumToPhotonFlux(xenon['power'], wlenshared)
    xenon_eye['photons']  = powerSpectrumToPhotonFlux(xenon_eye['power'], wlenshared)

    # Trim noise
    xenon.loc[(wlenshared < 300) | (wlenshared> 770), 'photons'] = 0
    xenon_eye.loc[(wlenshared < 300) | (wlenshared > 770), 'photons'] = 0

    return spec440, spec560, intensity440, intensity560, wlenshared, xenon, xenon_eye

def _run_quiet(cmd: str) -> None:
    """
    Description
    -----------
    Run a shell command with stdout/stderr suppressed; raise on failure.

    Parameters
    ----------
    None.

    Returns
    -------
    Nothing.
    """
    subprocess.run(
        ["bash", "-lc", cmd],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

def _find_arial_path():
    """
    Description
    -----------
    Return a path to an Arial .ttf if available, else None.
    Tries Matplotlib lookup first, then a filesystem search.

    Parameters
    ----------
    None.

    Returns
    -------
    Nothing.
    """
    # 1) Ask Matplotlib (no fallback: either Arial exists or we treat it as missing)
    try:
        p = fm.findfont("Arial", fallback_to_default=False)
        if os.path.exists(p):
            return p
    except Exception:
        pass

    # 2) Filesystem search (covers msttcorefonts and other installs)
    candidates = glob.glob("/usr/share/fonts/**/[Aa]rial*.ttf", recursive=True)
    if candidates:
        # Prefer msttcorefonts if present
        candidates.sort(key=lambda x: ("msttcorefonts" not in x.lower(), len(x)))
        return candidates[0]

    return None

def ensure_arial_font():
    """
    Description
    -----------
    Ensures the arial font is downloaded and active in the system.
    An error tree follows in the font is not installed in the system.

    Parameters
    ----------
    None.

    Returns
    -------
    Nothing.
    """

    # Fast path: already installed
    arial_path = _find_arial_path()
    if arial_path:
        fm.fontManager.addfont(arial_path)
        plt.rcParams["font.family"] = "Arial"
        print(f"Arial already available and active ({arial_path})")
        return

    # Install path: keep everything quiet
    try:
        # Pre-accept EULA
        try:
            _run_quiet(
                'echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula boolean true" | sudo debconf-set-selections'
            )
        except subprocess.CalledProcessError:
            _run_quiet(
                'echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula boolean true" | debconf-set-selections'
            )

        # Update + install + rebuild font cache
        try:
            _run_quiet("sudo apt-get update -qq")
            _run_quiet("sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq ttf-mscorefonts-installer")
            _run_quiet("sudo fc-cache -f")
        except subprocess.CalledProcessError:
            # Some environments don't need sudo
            _run_quiet("apt-get update -qq")
            _run_quiet("DEBIAN_FRONTEND=noninteractive apt-get install -y -qq ttf-mscorefonts-installer")
            _run_quiet("fc-cache -f")

        # Re-check + activate
        arial_path = _find_arial_path()
        if arial_path:
            fm.fontManager.addfont(arial_path)
            plt.rcParams["font.family"] = "Arial"
            print(f"Arial installed and active ({arial_path})")
        else:
            print("Installed mscorefonts, but Arial .ttf was not found on the system.")

    except Exception:
        print("Could not install/activate Arial (see runtime permissions/package manager availability).")

# This function runs the initialization protocol when called from the notebook
def load():
    """
    Description
    -----------
    Initializes the notebook environment with the following:
        - Standard scientific libraries.
        - Custom libraries.
        - Select experimental data.
        - Fonts for Plotting.
        - Data pathing.

    Parameters
    ----------
    None.

    Returns
    -------
    Initialization dictionary : python dict
        A dictionary containing the runtime environment to initalize the notebook.
        Contains values for all the features listed in the description (Writing them all here is redundant. Look at the return disctionary if you want to know whats in it.)
    """

    ensure_arial_font()
    mpl.rcParams['font.family'] = 'LiberationSerif-Regular'
    fonts = fm.findSystemFonts(fontpaths=None, fontext='ttf')


    spec440, spec560, intensity440, intensity560, wlenshared, xenon, xenon_eye = loadInSpectra()
    # wlen_xenon_corneal, intensity_xenon_corneal, other_xenon_spectrum, wlen_xenon, intensity_xenon = loadInActionSpectra()

    #Initialize fonts
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 14,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    # Resolve manuscript spectral assets (works from source and when frozen).
    filepath = ensure_dir(manuscript_dir())

    # Return a manually created python with packages and objects we specifically want in the model that weren't defined earlier.
    return {
        # Define a file path to reach data and save files
        "filepath": filepath,

        # Import all libraries used
        "os": os,
        "re": re,
        "csv": csv,
        "math": math,
        "time": time,
        "warnings": warnings,
        "datetime": datetime,
        "partial": partial,
        "np": np,
        "pd": pd,
        "mpl": mpl,
        "plt": plt,
        "sns": sns,
        # "pyabf": pyabf,
        "tqdm": tqdm,
        "Parallel": Parallel,
        "delayed": delayed,
        "PdfMerger": PdfMerger,
        "PdfPages": PdfPages,
        "signal": signal,
        "Bounds": Bounds,
        "minimize": minimize,
        "differential_evolution": differential_evolution,
        "basinhopping": basinhopping,
        "CubicSpline": CubicSpline,
        "PchipInterpolator": PchipInterpolator,
        "Akima1DInterpolator": Akima1DInterpolator,
        "interp1d": interp1d,
        "gp_minimize": gp_minimize,
        "widgets": widgets,
        "interact": interact,

        # Import stimulus spectra
        "spec440": spec440,
        "spec560": spec560,
        "intensity440": intensity440,
        "intensity560": intensity560,
        "wlenshared": wlenshared,
        "xenon": xenon,
        "xenon_eye":xenon_eye,

        # Fonts from Philippe's additions
        "fonts": fonts
    }