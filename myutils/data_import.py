# Standard library imports
import csv

# Third-party imports
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

filepath = r'.\data\manuscript'
figurepath = r'D:\projects\adaptation manuscript\figures\model fitting'
filepath_optimization_logs = 'Optimization Logs'

def sanityCheck():
    print("data_import.py Sanity Checked")

def resampleSpectrum(power, wlen, new_wlen):
    """
    Description
    -----------
    Interpolates a given power spectrum from original wavelength range wlen to new_wlen. 

    Parameters
    ----------
    power : array-like, shape (1,)
        Array of power values in microwatts.
    wlen : array-like, shape (1,)
        Array of wavelengths in nanometers. Must have the same length as 'power.
    new_wlen : array-like, shape (1,)
        Array of new wavelengths to interpolate to, in nanometers

    Returns
    -------
    imposed_powers : array-like, shape (1,)
        Array of imposed powers in the standard range.
    """
    f = interp1d(wlen, power, kind='linear', bounds_error=False, fill_value='extrapolate')
    imposed_powers = f(new_wlen)

    return imposed_powers

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
        Array of photon counts. Same shape as 'power' and 'wlen'.
    """

    # Set power measurements < 0 to 0
    power[power < 0] = 0

    h = 6.626*10**(-34) # Planck's constant
    c = 2.9979*10**8 # Speed of light
    single_photon_energy = [h*c / (w*10**(-9)) for w in wlen]
    nphotons = np.array([(power_uw * 10**-6) / ev * 10**-8 for power_uw, ev in zip(power, single_photon_energy)])

    return nphotons

def inVivoTransmission(spec_wlen, spec_intensity, trans_wlen, trans_val, tail_nm=100, clamp_min=0.0):
    """
    Apply filter transmission to spectrum(s).

    Parameters
    ----------
    spec_wlen : (N,) array
        Wavelengths of spectrum
    spec_intensity : (N,) or (N, M) array
        Spectrum intensity values. Each entry must have the same length as 'spec_wlen'.
    trans_wlen : (K,) array
        Wavelengths of transmission curve
    trans_val : (K,) array
        Transmission values. Must have the same length as 'trans_wlen'.
    tail_nm : float
        Fit linear extrapolation using last `tail_nm` nm of transmission.
    clamp_min : float or None
        Clamp transmission to be >= clamp_min.

    Returns
    -------
    spec_invivo : (N,) array
        same shape as spec_intensity
    transmission_interp : (N,) array
    """

    spec_wlen = np.asarray(spec_wlen).astype(float)
    spec_intensity = np.asarray(spec_intensity).astype(float)
    trans_wlen = np.asarray(trans_wlen).astype(float)
    trans_val = np.asarray(trans_val).astype(float)

    if spec_wlen.ndim != 1:
        raise ValueError("spec_wlen must be 1D.")
    if trans_wlen.ndim != 1:
        raise ValueError("trans_wlen must be 1D.")
    if trans_val.ndim != 1:
        raise ValueError("trans_val must be 1D.")
    if spec_wlen.shape[0] != spec_intensity.shape[0]:
        raise ValueError("spec_wlen and spec_intensity must have same length.")

    # Ensure transmission curve sorted
    order = np.argsort(trans_wlen)
    trans_wlen = trans_wlen[order]
    trans_val = trans_val[order]

    # --- Linear fit on tail ---
    max_wl = trans_wlen.max()
    tail_mask = trans_wlen >= (max_wl - tail_nm)

    if tail_mask.sum() < 2:
        raise ValueError("Not enough tail points for linear fit.")

    m, b = np.polyfit(trans_wlen[tail_mask], trans_val[tail_mask], 1)

    # --- Interpolate within range ---
    transmission_interp = np.interp(
        spec_wlen,
        trans_wlen,
        trans_val,
        left=np.nan,
        right=np.nan
    )

    # --- Linear extrapolation above range ---
    high_mask = spec_wlen > trans_wlen.max()
    transmission_interp[high_mask] = m * spec_wlen[high_mask] + b

    # Optional clamp
    if clamp_min is not None:
        transmission_interp = np.clip(transmission_interp, clamp_min, None)

    # --- Apply transmission ---
    if spec_intensity.ndim == 1:
        spec_invivo = spec_intensity * transmission_interp
    elif spec_intensity.ndim == 2:
        spec_invivo = spec_intensity * transmission_interp[:, None]
    else:
        raise ValueError("spec_intensity must be 1D or 2D.")

    return spec_invivo, transmission_interp


def _infer_lucas_wlen_col(df):
    """Infer wavelength column from Lucas plate reader DataFrame (case-insensitive)."""
    for name in ['wavelength', 'wl', 'wlen', 'nm', 'lambda', 'wavelength (nm)']:
        for col in df.columns:
            if str(col).strip().lower() == name.replace(' ', ''):
                return col
            if name in str(col).strip().lower():
                return col
    # First column if numeric and looks like wavelength (300-800 range)
    first = df.iloc[:, 0]
    if pd.api.types.is_numeric_dtype(first):
        if first.min() >= 250 and first.max() <= 900:
            return df.columns[0]
    return None


def load_lucas_plate_reader_excel(
    filepath,
    sheet_name=0,
    header_row=11,
    wlen_col=None,
    power_cols=None,
    target_row=1,
):
    """
    Load spectral data from a Lucas lab plate reader Excel file.

    The first 11 rows are metadata (spectrometer measurement method) and are
    skipped. The 12th row (0-based index 11) is the header: 1st column =
    wavelength, remaining columns = power. Power column headers can indicate
    units (e.g. "Ee [W/(sqm*nm)]" for irradiance in W/(m²·nm)). This function
    converts W/(m²·nm) to µW/(µm²·nm) for use with powerSpectrumToPhotonFlux.

    Parameters
    ----------
    filepath : str or path-like
        Path to the .xlsx file (requires openpyxl: pip install openpyxl).
    sheet_name : int or str, optional
        Sheet to read. Default 0 (first sheet).
    header_row : int, optional
        Zero-based row index of the header row. Default 11 (12th row). Rows
        before this are metadata and skipped.
    wlen_col : str or int, optional
        Column name or index for wavelength. If None, first column is used.
    power_cols : list of str or int, optional
        Column names or indices for power spectra. If None, all columns except
        the wavelength column are used (in order).
    target_row : int, optional
        Zero-based row index in the raw spreadsheet containing target
        intensities in log10(photons/(cm²·s)). Default 1 (second row of sheet).

    Returns
    -------
    spectra : list of (wlen, power) tuples
        Each element is (wavelength_1d_array, power_1d_array). Wavelength in nm,
        power in µW/(µm²·nm).
    target_log : ndarray of shape (n_spectra,)
        Target integrated intensity in log10(photons/(cm²·s)) for each spectrum.
    """
    # Read target log intensities from header area (raw spreadsheet second row = index 1)
    df_raw = pd.read_excel(filepath, sheet_name=sheet_name, header=None)
    target_log = pd.to_numeric(df_raw.iloc[target_row, 1:], errors="coerce").values
    if not np.all(np.isfinite(target_log)):
        raise ValueError(
            f"Target row (raw row index {target_row}) must contain numeric "
            "log10(photons/(cm²·s)) for each power column (columns 1 onward)."
        )

    # Read data with header row; spectral data starts at header_row + 1
    df = pd.read_excel(filepath, sheet_name=sheet_name, header=header_row)
    df = df.dropna(how='all').copy()

    if wlen_col is None:
        wlen_col = df.columns[0]
    if power_cols is None:
        power_cols = [c for c in df.columns[1:] if c in df.columns]
    if not power_cols:
        raise ValueError("No power columns found. Pass power_cols= (e.g. list of column names).")
    if len(target_log) != len(power_cols):
        raise ValueError(
            f"Target row has {len(target_log)} values but there are {len(power_cols)} power columns."
        )

    # Determine conversion from power column headers (e.g. "Ee [W/(sqm*nm)]" -> W/(m²·nm))
    # W/(m²·nm) -> µW/(µm²·nm): 1 W/m² = 1e6 µW / 1e12 µm² = 1e-6 µW/µm²
    first_power_header = str(power_cols[0])
    if "W/(sqm*nm)" in first_power_header or "W/(m2*nm)" in first_power_header.lower():
        power_scale = 1e-6
    else:
        raise ValueError(
            f"Unrecognized power unit in header '{first_power_header}'. "
            "Expected e.g. 'Ee [W/(sqm*nm)]' for W/(m²·nm)."
        )

    # Coerce to numeric so any stray non-numeric rows are dropped
    wlen = pd.to_numeric(df[wlen_col], errors="coerce").values
    valid = np.isfinite(wlen)

    # Restrict to rows that are numeric in wavelength and in every power column
    for col in power_cols:
        p = pd.to_numeric(df[col], errors="coerce").values
        valid = valid & np.isfinite(p)

    wlen = wlen[valid]

    spectra = []
    for col in power_cols:
        p = pd.to_numeric(df[col], errors="coerce").values
        p = p[valid]
        p = np.where(np.isfinite(p), p, 0.0)
        p = p * power_scale
        spectra.append((wlen.copy(), p))

    return spectra, target_log


def get_lucas_stimulus_options(
    filepath,
    sheet_name=0,
    wavelength_row=0,
    target_row=1,
):
    """
    Read the Lucas plate reader Excel to get (peak wavelength, target intensity)
    for each stimulus column from the raw spreadsheet top and second row.

    Parameters
    ----------
    filepath : str or path-like
        Path to the .xlsx file.
    sheet_name : int or str, optional
        Sheet to read. Default 0.
    wavelength_row : int, optional
        Zero-based row index in the raw sheet for peak wavelength (nm) per column.
        Default 0 (top row).
    target_row : int, optional
        Zero-based row index for target intensity in log10(photons/(cm²·s)).
        Default 1 (second row).

    Returns
    -------
    options : list of (float, float)
        For each power column (same order as load_lucas / lucas_spectra_to_stimulus_lists),
        the pair (peak_wavelength_nm, target_log_intensity).
    """
    df_raw = pd.read_excel(filepath, sheet_name=sheet_name, header=None)
    # Column 0 = wavelength; columns 1 onward = one spectrum each
    n_cols = len(df_raw.columns) - 1
    if n_cols < 1:
        raise ValueError("Excel must have at least one power column (besides wavelength).")
    # Strip trailing 'nm' from wavelength row (e.g. '365nm' or '365 nm' -> 365)
    raw_wlen = (
        df_raw.iloc[wavelength_row, 1:]
        .astype(str)
        .str.strip()
        .str.replace(r"\s*nm\s*$", "", case=False, regex=True)
        .str.strip()
    )
    wlen_vals = pd.to_numeric(raw_wlen, errors="coerce").values
    target_vals = pd.to_numeric(df_raw.iloc[target_row, 1:], errors="coerce").values
    if not np.all(np.isfinite(wlen_vals)) or not np.all(np.isfinite(target_vals)):
        raise ValueError(
            f"Row {wavelength_row} (wavelength) and row {target_row} (intensity) must "
            "contain numeric values for each power column."
        )
    return list(zip(wlen_vals, target_vals))


def lucas_stimulus_options_table(options):
    """
    Convert the list of (wavelength_nm, target_log) from get_lucas_stimulus_options
    into a readable DataFrame with index, wavelength, and intensity per stimulus.

    Parameters
    ----------
    options : list of (float, float)
        Return value of get_lucas_stimulus_options(filepath, ...).

    Returns
    -------
    pd.DataFrame
        Columns: index, wavelength_nm, target_log. Use print(df) or display(df).
    """
    rows = [
        {"index": i, "wavelength_nm": float(w), "target_log": float(t)}
        for i, (w, t) in enumerate(options)
    ]
    return pd.DataFrame(rows)


def get_lucas_stimulus_index(
    filepath,
    peak_wavelength_nm,
    target_log_intensity,
    sheet_name=0,
    wavelength_row=0,
    target_row=1,
    tol_wlen_nm=2.0,
    tol_log=0.05,
):
    """
    Return the stimulus column index for the given peak wavelength and target intensity.

    Uses the raw Excel top row for wavelength and second row for intensity (same
    as get_lucas_stimulus_options). The first column that matches within tolerance
    is returned.

    Parameters
    ----------
    filepath : str or path-like
        Path to the .xlsx file.
    peak_wavelength_nm : float
        Peak wavelength in nm (matched to wavelength_row).
    target_log_intensity : float
        Target intensity in log10(photons/(cm²·s)) (matched to target_row).
    sheet_name, wavelength_row, target_row
        Passed to get_lucas_stimulus_options.
    tol_wlen_nm : float, optional
        Max absolute difference in nm to consider wavelength a match. Default 2.
    tol_log : float, optional
        Max absolute difference in log units for intensity. Default 0.05.

    Returns
    -------
    index : int
        Zero-based index of the matching stimulus (for use with spec_wlens[i], etc.).
    """
    options = get_lucas_stimulus_options(
        filepath, sheet_name=sheet_name, wavelength_row=wavelength_row, target_row=target_row
    )
    for i, (w, t) in enumerate(options):
        if abs(w - peak_wavelength_nm) <= tol_wlen_nm and abs(t - target_log_intensity) <= tol_log:
            return i
    raise ValueError(
        f"No stimulus found with peak wavelength {peak_wavelength_nm} nm and target log "
        f"{target_log_intensity}. Available: {options[:10]}{'...' if len(options) > 10 else ''}"
    )


def lucas_spectra_to_stimulus_lists(
    filepath,
    specwlen=None,
    sheet_name=0,
    header_row=11,
    target_row=1,
    wlen_col=None,
    power_cols=None,
):
    """
    Convert Lucas lab plate reader Excel spectra into stimulus lists compatible
    with the notebook (Broadband light section): spec_wlens and spec_intensities.

    The raw spreadsheet second row (header area) gives target intensities in
    log10(photons/(cm²·s)); these are converted to photons/(µm²·s) and each
    spectrum is scaled with
    scaleToPhotonCount so its integrated photon count matches the target. No
    pre-receptoral filtering is applied.

    Parameters
    ----------
    filepath : str or path-like
        Path to the Lucas lab .xlsx file.
    specwlen : array-like, optional
        Wavelength grid (nm) for resampling. If None, uses np.linspace(300, 750, 1000).
    sheet_name : int or str, optional
        Excel sheet to read (see load_lucas_plate_reader_excel).
    header_row : int, optional
        Zero-based row index of the header row (default 11 = 12th row). Rows
        before this are metadata and skipped.
    target_row : int, optional
        Raw spreadsheet row index (0-based) containing target log10(photons/(cm²·s)).
        Default 1 (second row of sheet).
    wlen_col : str or int, optional
        Wavelength column (see load_lucas_plate_reader_excel).
    power_cols : list, optional
        Power columns (see load_lucas_plate_reader_excel).

    Returns
    -------
    spec_wlens : list of ndarray
        One wavelength array per stimulus (same grid for all if specwlen given).
    spec_intensities : list of ndarray
        Photon flux density (photons/micron^2/s/nm) for each stimulus, same length as spec_wlens.
    """
    if specwlen is None:
        specwlen = np.linspace(350, 750, 1000)
    specwlen = np.asarray(specwlen, dtype=float)

    spectra, target_log = load_lucas_plate_reader_excel(
        filepath,
        sheet_name=sheet_name,
        header_row=header_row,
        target_row=target_row,
        wlen_col=wlen_col,
        power_cols=power_cols,
    )

    # log10(photons/(cm²·s)) -> photons/(µm²·s): 1 cm² = 1e8 µm²
    target_photons_per_um2_s = np.power(10.0, target_log - 8.0)

    from myutils.model_functions import scaleToPhotonCount

    spec_wlens = []
    spec_intensities = []

    for i, (wlen, power) in enumerate(spectra):
        power = np.asarray(power, dtype=float)
        wlen = np.asarray(wlen, dtype=float)
        power_resampled = resampleSpectrum(power, wlen, specwlen)
        # Copy so powerSpectrumToPhotonFlux does not modify in place
        intensity = np.array(powerSpectrumToPhotonFlux(power_resampled.copy(), specwlen))
        # Scale to target photon count (photons/(µm²·s))
        wlens_scaled, scaled = scaleToPhotonCount(
            specwlen, intensity, [float(target_photons_per_um2_s[i])]
        )
        spec_wlens.append(wlens_scaled[0])
        spec_intensities.append(scaled[0])

    return spec_wlens, spec_intensities


def getDataset(datafilename):
    """
    Description
    -----------
    Put every row from a csv file into a python list.

    Parameters
    ----------
    datafilename : string
        File path to the specified data file.
    Returns
    -------
    data : python list (1,)
        The python list where each index is a row from the data csv file.
    """

    if len(datafilename) > 0:
        data = []
        with open(datafilename) as file:
            reader = csv.reader(file, delimiter = ',')
            for row in reader:
                data.append(float(row[0]))

    return data
