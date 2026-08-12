"""
Manuscript stimulus runner (desktop GUI).

Runs the six intensity-response protocols from the tutorial section
"Simulations from our manuscript" (see Melanopsin_Model_Tutorial.ipynb):
same ``setupModelRun`` / ``predictMelanopsin`` / ``plotModelRun`` pipeline.

**Run controls**

- **Stop** aborts the current integration (cooperative cancel); the previous completed
  run (if any) remains available for export.

**How to run**

- From the repository root (so ``data/`` resolves): ::

    python -m myutils.melanopsin_gui

- Requires the same spectral assets as the tutorial (``data/manuscript/spectra/440nm spectrum.csv``,
  ``560nm spectrum.csv``, ``Xenon.txt``, ``XeEye.txt``).

Internally, stimuli use ``stimulustype`` values such as ``\"spectrum matched to 440\"``
and ``\"xenon\"``; dropdown labels (``440-440``, ``Xenon with ocular filtering``, …)
are display names mapped to those calls.

**File menu**

- **Export Model Run** saves the last completed simulation (run **Run model** first).
  Choose **CSV** for a comma-separated table (``time_s`` plus every 1D time series
  aligned with ``xaxis``), or **NPZ** for a compressed NumPy archive of those arrays
  plus ``stimulus_label``. Nested ``params`` from the model are not included; re-run
  the simulation if you need parameters from the notebook API.
- **Import Spectrum** reads a spectrum file into the user spectrum library.
- **Import Stimulus** reads a v2 stimulus JSON (same schema as
  ``data/user_library/stimuli/``) and registers it in the stimulus library.

**Stimulus Builder** (``Data`` → **Stimulus Builder**)

- **Save** writes the current interval grid into the user stimulus library
  (``data/user_library/stimuli/``) and adds the name to the run dropdown.
- **Load from library…** reloads the user stimulus library from disk and opens
  a picker of saved custom stimuli to edit in the builder.

**Custom Stimulus Builder** (legacy menu placeholder)

- **Load stimulus…** reads a JSON preset and sets the main-window stimulus selection.
  The file must contain either ``\"label\"`` or ``\"stimulus\"`` with one of the same
  strings as the manuscript labels, for example: ``{\"label\": \"440-440\"}``.

**Spectrum Builder** (``Data`` → **Spectrum Builder**)

- **Monochromatic Spectrum** — single-bin impulse on the model ``wlenshared`` grid at the sample
  nearest your target wavelength (readout shows the resolved nm). Saved entries go
  into the custom spectrum library like file imports.
- **Mixed Spectrum** — linear combination of library spectra (built-in + custom).
  Each row can optionally be normalized to its peak (default) or integrated photon
  count before weighting; with neither option selected, raw spectral values are used.
  Preview and **Save to library…** use photon-flux units consistent with the Spectrum
  Library when unnormalized.

**Configure menu**

- **Model parameters…** opens an editor for every key in ``get_predict_melanopsin_defaults()``:
  edit values and click **Update** to apply them to the **next** ``predictMelanopsin`` run.
  **Save Parameters** / **Load Preset** write and read ``.json`` presets under
  ``data/user_library/configurations/``.

**Plot menu**

- **Compare saved predictions…** overlays fields from one or more exported (or current)
  model runs on shared axes.
- **Data Comparator…** overlays experimental data (headered or headerless CSV/Excel;
  time column or user sampling rate) on a model run (current run or saved CSV/NPZ
  export): data are interpolated onto the model timebase with ``syncTimeseries``,
  then plotted with ``plotModelVsData``. The model series defaults to
  ``currentGlobalGain`` but can be changed in the dialog.

**Settings menu**

- **Auto-save options** opens a window with two off-by-default checkboxes: after each
  successful run, optionally write the figure and/or CSV time series to folders of your
  choice (defaults ``./outputs/images`` and ``./outputs/predictions``; **Choose folder…** and
  **Reset to default**). When a checkbox is on, the current save path is shown.

**Help menu**

- **Open README…** opens the local ``README.md`` in the OS default app.
- **Open GitHub…** opens the project repository page in the default browser.
"""

from __future__ import annotations

import csv
import json
import math
import numbers
import os
import re
import subprocess
import sys
import textwrap
import webbrowser
from datetime import datetime
from pathlib import Path
import matplotlib

matplotlib.use("TkAgg")

import threading
import tkinter as tk
import tkinter.font as tkfont
from queue import Empty, Queue
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog
from tkinter import ttk

import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import SpanSelector

from myutils import library_storage
from myutils._version import __version__
from myutils.app_paths import data_dir, outputs_dir, project_root
from myutils.general_functions import syncTimeseries
from myutils.init import load
from myutils.library_storage import (
    BUILTIN_SPECTRUM_NAMES as _BUILTIN_STIMULUS_SPECTRA,
    STIMULUS_SCHEMA_VERSION as _STIMULUS_SCHEMA_VERSION,
)
from myutils.model_functions import (
    SimulationCancelled,
    get_predict_melanopsin_defaults,
    powerSpectrumToPhotonFlux,
    predictMelanopsin,
    setupModelRun,
)
from myutils.plotting_functions import plotModelRun, plotModelVsData
from myutils.startup_splash import StartupSplash

# Same rate as Melanopsin_Model_Tutorial.ipynb manuscript cells
MANUSCRIPT_RATE = 0.005

# Fallback sampling rate (Hz) when a dataset has no time column; used as the
# SamplingRateDialog prefill and the initial editable rate in Data Comparator.
DATA_COMPARE_SAMPLE_RATE_HZ = 100.0
DEFAULT_DATA_COMPARE_SERIES = "currentGlobalGain"

LABELS = (
    "440-440",
    "560-560",
    "440-560",
    "560-440",
    "Xenon",
    "Xenon with ocular filtering",
)
ADD_NEW_STIMULUS_LABEL = "Add stimulus"
NO_STIMULUS_SELECTED_LABEL = "Select stimulus to run..."

# Stimulus Library preview: auto-load scrubber below this interval count.
_STIMULUS_PREVIEW_AUTO_MAX_INTERVALS = 50

# All asset/output locations are anchored to the repository root (resolved via
# myutils.app_paths) so they work identically from a source checkout and from
# the frozen executable, regardless of the current working directory.
DEFAULT_AUTOSAVE_IMAGES_DIR = outputs_dir() / "images"
DEFAULT_AUTOSAVE_DATA_DIR = outputs_dir() / "predictions"

_WAVELENGTH_ALIASES = frozenset({
    "wavelength", "wavelengths", "wlen", "wl", "wave", "lambda",
    "nm", "wav", "wavelen", "wavel",
})
_INTENSITY_ALIASES = frozenset({
    "intensity", "intensities", "power", "counts", "signal", "flux",
    "irradiance", "photons", "int", "value",
})


def _detect_spectrum_columns(df: "pd.DataFrame") -> tuple[str, str] | None:
    """Find wavelength and intensity column names, case/whitespace-insensitive."""
    wlen_col: str | None = None
    intensity_col: str | None = None
    for col in df.columns:
        norm = str(col).strip().lower()
        if norm in _WAVELENGTH_ALIASES and wlen_col is None:
            wlen_col = col
        elif norm in _INTENSITY_ALIASES and intensity_col is None:
            intensity_col = col
    if wlen_col is None or intensity_col is None:
        return None
    return wlen_col, intensity_col


def _columns_look_headerless(columns) -> bool:
    """True if column labels look like data (or Unnamed) rather than names."""
    cols = list(columns)
    if len(cols) < 2:
        return False
    labels = [str(c).strip() for c in cols]
    if all(lab.lower().startswith("unnamed") for lab in labels):
        return True
    try:
        float(labels[0])
        float(labels[1])
        return True
    except ValueError:
        return False


def _infer_spectrum_text_sep(path: str) -> str:
    """Choose a delimiter for CSV/TSV/TXT spectra by peeking at the first line.

    Tab-separated ``.txt`` files (e.g. manuscript Xenon spectra) must not be
    read with the default comma separator: that succeeds with a single column
    whose name still contains an embedded tab.
    """
    lower = path.lower()
    if lower.endswith(".tsv"):
        return "\t"
    if lower.endswith(".csv"):
        return ","
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            line = ""
            for raw in f:
                line = raw.strip()
                if line:
                    break
    except OSError:
        return ","
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    if len(line.split()) >= 2:
        return r"\s+"
    return ","


def _load_spectrum_dataframe(path: str, *, header=0) -> "pd.DataFrame":
    """Load a spectrum spreadsheet; ``header=None`` for headerless files."""
    lower = path.lower()
    if lower.endswith(".xlsx") or lower.endswith(".xls"):
        return pd.read_excel(path, header=header)
    sep = _infer_spectrum_text_sep(path)
    kwargs: dict = {"header": header, "sep": sep}
    if sep == r"\s+":
        kwargs["engine"] = "python"
    return pd.read_csv(path, **kwargs)


_HEADERLESS_SPECTRUM_NOTE = (
    "No column names were found. Assuming the leftmost column is wavelength "
    "and the next column is intensity."
)


def _read_spectrum_file(path: str) -> tuple[np.ndarray, np.ndarray, str | None]:
    """Read a spreadsheet and return (wlen, intensity, note).

    Accepts CSV, TSV, .xlsx/.xls.  Column names are matched
    case-insensitively against common wavelength and intensity aliases.
    If no names are found and the file looks headerless, column 0 is
    wavelength and column 1 is intensity; ``note`` explains that assumption.
    """
    df = _load_spectrum_dataframe(path)
    cols = _detect_spectrum_columns(df)
    note: str | None = None
    if cols is None:
        if not _columns_look_headerless(df.columns):
            raise ValueError(
                "Could not find wavelength and intensity columns.\n"
                f"Found columns: {list(df.columns)}\n\n"
                "Expected a column named like 'wavelength' (or 'wl', 'wlen', …)\n"
                "and one named like 'intensity' (or 'intensities', 'power', …)."
            )
        df = _load_spectrum_dataframe(path, header=None)
        if df.shape[1] < 2:
            raise ValueError(
                "Headerless spectrum file must have at least two columns "
                "(wavelength, intensity)."
            )
        wlen_col, intensity_col = df.columns[0], df.columns[1]
        note = _HEADERLESS_SPECTRUM_NOTE
    else:
        wlen_col, intensity_col = cols
    wlen = np.asarray(df[wlen_col].to_numpy(), dtype=float)
    intensity = np.asarray(df[intensity_col].to_numpy(), dtype=float)
    return wlen, intensity, note




def _power_to_photon_flux(power: np.ndarray, wlen: np.ndarray) -> np.ndarray:
    """Convert a power spectrum to photon flux (negatives clamped to 0).

    Copies the input so the shared call into ``powerSpectrumToPhotonFlux`` (which
    zeroes negatives in place) does not mutate caller-owned arrays.
    """
    power_arr = np.asarray(power, dtype=float).copy()
    wlen_arr = np.asarray(wlen, dtype=float)
    return np.asarray(powerSpectrumToPhotonFlux(power_arr, wlen_arr), dtype=float)


# ---------------------------------------------------------------------------
# Spectrum y-axis unit conversion
# ---------------------------------------------------------------------------
#
# The model consumes spectra in photons/(µm²·nm·s). Imported spectra can arrive
# in a variety of spectrophotometer/spectroradiometer units, so every import is
# converted to that canonical unit before it is stored.
#
# ``kind`` is either:
#   * "power"   - a spectral irradiance (energy) unit. ``scale`` converts it to
#     µW/(cm²·nm); the h·c/λ conversion is then applied via
#     ``powerSpectrumToPhotonFlux`` (which maps µW/(cm²·nm) -> photons/(µm²·nm·s)).
#   * "photon"  - already a photon-flux density. ``scale`` converts it directly
#     to photons/(µm²·nm·s) (no h·c/λ factor).
#
# Avogadro's number (photons per mole) for the µmol option.
_AVOGADRO = 6.02214076e23

# Storage key kept as ``photons_per_um2_nm`` for backward compatibility; the
# physical unit includes /s (see display labels below).
_SPECTRUM_UNIT_INTENSITY = "photons_per_um2_nm"
_SPECTRUM_UNIT_POWER = "uw_per_cm2_nm"

# key -> (label, kind, scale)
_SPECTRUM_UNIT_TABLE: dict[str, tuple[str, str, float]] = {
    # Main options.
    _SPECTRUM_UNIT_INTENSITY: ("Intensity (photons/\u00b5m\u00b2/nm/s)", "photon", 1.0),
    _SPECTRUM_UNIT_POWER: ("Power (\u00b5W/cm\u00b2/nm)", "power", 1.0),
    # "Other" dropdown - power (energy) units, scaled to µW/(cm²·nm).
    "w_per_m2_nm": ("W/m\u00b2/nm", "power", 1e2),
    "mw_per_m2_nm": ("mW/m\u00b2/nm", "power", 1e-1),
    "uw_per_m2_nm": ("\u00b5W/m\u00b2/nm", "power", 1e-4),
    "w_per_cm2_nm": ("W/cm\u00b2/nm", "power", 1e6),
    "mw_per_cm2_nm": ("mW/cm\u00b2/nm", "power", 1e3),
    "erg_per_s_cm2_nm": ("erg/s/cm\u00b2/nm", "power", 1e-1),
    # "Other" dropdown - photon-flux units, scaled to photons/(µm²·nm·s).
    "photons_per_s_cm2_nm": ("photons/s/cm\u00b2/nm", "photon", 1e-8),
    "photons_per_s_m2_nm": ("photons/s/m\u00b2/nm", "photon", 1e-12),
    "umol_per_m2_s_nm": (
        "\u00b5mol photons/m\u00b2/s/nm",
        "photon",
        _AVOGADRO * 1e-6 * 1e-12,
    ),
}

# Units listed in the "Other" dropdown (everything but the two main options),
# preserving insertion order for a stable menu.
_SPECTRUM_OTHER_UNIT_KEYS: list[str] = [
    key
    for key in _SPECTRUM_UNIT_TABLE
    if key not in (_SPECTRUM_UNIT_INTENSITY, _SPECTRUM_UNIT_POWER)
]


def _convert_spectrum_to_photon_flux(
    wlen: np.ndarray, y: np.ndarray, unit_key: str
) -> np.ndarray:
    """Convert an uploaded y-axis array to photons/(µm²·nm·s).

    ``unit_key`` must be a key of ``_SPECTRUM_UNIT_TABLE``. Power-type units are
    rescaled to µW/(cm²·nm) and passed through ``powerSpectrumToPhotonFlux``;
    photon-type units are rescaled directly.
    """
    try:
        _label, kind, scale = _SPECTRUM_UNIT_TABLE[unit_key]
    except KeyError as exc:
        raise ValueError(f"Unknown spectrum unit {unit_key!r}.") from exc
    wlen_arr = np.asarray(wlen, dtype=float)
    y_arr = np.asarray(y, dtype=float).copy()
    if kind == "power":
        return _power_to_photon_flux(y_arr * scale, wlen_arr)
    # Photon-flux density: clamp negatives (measurement noise) then rescale.
    y_arr[y_arr < 0] = 0.0
    return y_arr * scale




class SpectrumUnitsDialog(tk.Toplevel):
    """Modal prompt asking for the intensity units of an uploaded spectrum.

    Wavelength is always assumed to be in nm. Exactly one intensity unit is
    selectable via three mutually exclusive options ("Intensity", "Power",
    "Other"); choosing "Other" enables a dropdown of additional
    spectrophotometer units. The chosen canonical unit key is available in
    ``self.result`` after the dialog closes (``None`` when cancelled).
    """

    _INTENSITY = _SPECTRUM_UNIT_INTENSITY
    _POWER = _SPECTRUM_UNIT_POWER
    _OTHER = "__other__"

    def __init__(self, parent):
        super().__init__(parent)
        self.result: str | None = None
        self.title("Spectrum intensity units")
        self.transient(parent)
        self.resizable(False, False)

        self._choice_var = tk.StringVar(master=self, value=self._INTENSITY)
        self._other_labels = [
            _SPECTRUM_UNIT_TABLE[k][0] for k in _SPECTRUM_OTHER_UNIT_KEYS
        ]
        self._label_to_key = {
            _SPECTRUM_UNIT_TABLE[k][0]: k for k in _SPECTRUM_OTHER_UNIT_KEYS
        }
        self._other_var = tk.StringVar(
            master=self,
            value=self._other_labels[0] if self._other_labels else "",
        )

        frame = ttk.Frame(self, padding=16)
        frame.grid(row=0, column=0, sticky=tk.NSEW)

        ttk.Label(
            frame,
            text=(
                "What are the intensity units of this spectrum?\n"
                "(Wavelength must be in nm.)\n\n"
                "The model runs in photons/\u00b5m\u00b2/nm/s; the data will be "
                "converted for you."
            ),
            justify=tk.LEFT,
        ).grid(row=0, column=0, sticky=tk.W, pady=(0, 12))

        ttk.Radiobutton(
            frame,
            text=_SPECTRUM_UNIT_TABLE[self._INTENSITY][0],
            value=self._INTENSITY,
            variable=self._choice_var,
            command=self._sync_other_state,
        ).grid(row=1, column=0, sticky=tk.W)
        ttk.Radiobutton(
            frame,
            text=_SPECTRUM_UNIT_TABLE[self._POWER][0],
            value=self._POWER,
            variable=self._choice_var,
            command=self._sync_other_state,
        ).grid(row=2, column=0, sticky=tk.W)

        other_row = ttk.Frame(frame)
        other_row.grid(row=3, column=0, sticky=tk.W)
        ttk.Radiobutton(
            other_row,
            text="Other:",
            value=self._OTHER,
            variable=self._choice_var,
            command=self._sync_other_state,
        ).grid(row=0, column=0, sticky=tk.W)
        self._other_combo = ttk.Combobox(
            other_row,
            values=self._other_labels,
            textvariable=self._other_var,
            state="disabled",
            width=28,
        )
        self._other_combo.grid(row=0, column=1, sticky=tk.W, padx=(6, 0))

        button_row = ttk.Frame(frame)
        button_row.grid(row=4, column=0, sticky=tk.E, pady=(16, 0))
        ttk.Button(button_row, text="OK", command=self._on_ok).grid(
            row=0, column=0, padx=(0, 6)
        )
        ttk.Button(button_row, text="Cancel", command=self._on_cancel).grid(
            row=0, column=1
        )

        self._sync_other_state()
        self.protocol("WM_DELETE_WINDOW", self._on_cancel)
        self.bind("<Return>", lambda _e: self._on_ok())
        self.bind("<Escape>", lambda _e: self._on_cancel())
        self.grab_set()
        self.wait_window()

    def _sync_other_state(self) -> None:
        is_other = self._choice_var.get() == self._OTHER
        self._other_combo.configure(state="readonly" if is_other else "disabled")

    def _on_ok(self) -> None:
        choice = self._choice_var.get()
        if choice == self._OTHER:
            label = self._other_var.get()
            key = self._label_to_key.get(label)
            if key is None:
                messagebox.showwarning(
                    "Spectrum intensity units",
                    "Please select a unit from the dropdown.",
                    parent=self,
                )
                return
            self.result = key
        else:
            self.result = choice
        self.destroy()

    def _on_cancel(self) -> None:
        self.result = None
        self.destroy()


def _prompt_spectrum_units(parent) -> str | None:
    """Show the units dialog; return a canonical unit key or None if cancelled."""
    return SpectrumUnitsDialog(parent).result






def _validate_stimulus_spec(data: object) -> dict:
    """Validate and normalize a parsed v2 stimulus spec."""
    return library_storage.validate_stimulus_spec(data)


def _load_stimulus_spec_from_path(path: str) -> dict:
    """Read and validate a stimulus spec from an arbitrary JSON file path."""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return _validate_stimulus_spec(data)


def _serialize_stimulus_spec(
    name: str,
    blocks: list[dict],
    total_duration: float,
    custom_spectra: dict[str, tuple[np.ndarray, np.ndarray]],
) -> dict:
    """Validate ``blocks`` and return a JSON-friendly v2 stimulus spec.

    For interval blocks, any referenced custom spectrum must already exist in
    ``custom_spectra`` -- the spectrum library stores the array bytes and we
    only persist the name reference here.
    """
    name = str(name).strip()
    if not name:
        raise ValueError("Stimulus name is required.")
    if not blocks:
        raise ValueError("Stimulus must have at least one block.")
    for i, block in enumerate(blocks):
        if block.get("type") == "interval":
            spectrum = str(block.get("spectrum_ref", "")).strip()
            if not spectrum:
                raise ValueError(f"Block {i + 1}: spectrum is required.")
            if (
                spectrum not in _BUILTIN_STIMULUS_SPECTRA
                and spectrum not in custom_spectra
            ):
                raise ValueError(
                    f"Block {i + 1}: spectrum {spectrum!r} is not a built-in "
                    "or known custom spectrum; cannot save."
                )
    spec = {
        "version": _STIMULUS_SCHEMA_VERSION,
        "name": name,
        "total_duration": float(total_duration),
        "blocks": blocks,
    }
    return _validate_stimulus_spec(spec)


def _resample_spectrum(
    src_wlen: np.ndarray,
    src_int: np.ndarray,
    target_wlen: np.ndarray,
) -> np.ndarray:
    """Linearly resample ``src_int`` onto ``target_wlen``; zero-fill out of range."""
    src_wlen = np.asarray(src_wlen, dtype=float)
    src_int = np.asarray(src_int, dtype=float)
    target_wlen = np.asarray(target_wlen, dtype=float)
    order = np.argsort(src_wlen)
    return np.interp(target_wlen, src_wlen[order], src_int[order], left=0.0, right=0.0)


def _scale_spectrum_to_integrated_intensity(
    src_wlen: np.ndarray,
    src_intensity: np.ndarray,
    target_intensity: float,
    *,
    target_wlen: np.ndarray | None = None,
    ref_name: str = "spectrum",
) -> np.ndarray:
    """Normalize a spectrum and scale so ∫ I(λ) dλ equals ``target_intensity``.

    ``target_intensity`` is total photons/µm²/s. When ``target_wlen`` is given,
    the spectrum is resampled onto that grid first (model protocol path);
    otherwise scaling stays on the source wavelength grid (builder preview).
    """
    out_wlen = (
        np.asarray(target_wlen, dtype=float)
        if target_wlen is not None
        else np.asarray(src_wlen, dtype=float)
    )
    if float(target_intensity) <= 0:
        return np.zeros_like(out_wlen, dtype=float)

    if target_wlen is not None:
        resampled = _resample_spectrum(src_wlen, src_intensity, out_wlen)
    else:
        resampled = np.asarray(src_intensity, dtype=float)

    dwlen = float(np.mean(np.diff(out_wlen))) if out_wlen.size > 1 else 1.0
    photon_count = float(np.trapezoid(resampled, out_wlen, dwlen))
    if not np.isfinite(photon_count) or photon_count <= 0:
        raise ValueError(
            f"Spectrum {ref_name!r} integrates to zero photons; "
            "cannot scale to a positive intensity."
        )
    return (resampled / photon_count) * float(target_intensity)


def _interval_names_from_spec(spec: dict) -> list[str]:
    """Return one display name per expanded interval in ``spec``."""
    names: list[str] = []
    for block in spec.get("blocks", []):
        if block.get("type") == "interval":
            names.append(str(block.get("spectrum_ref", "Interval")))
    return names


def _library_spectrum_arrays(
    app: "ManuscriptSimApp", name: str
) -> tuple[np.ndarray, np.ndarray] | None:
    """Return ``(wlen, intensity)`` for a library spectrum name, or None if unavailable."""
    if name == "Dark":
        return _DARK_WAVELENGTHS, _DARK_INTENSITY
    env = app._env
    if env is None:
        return None
    wlen = np.asarray(env["wlenshared"], dtype=float)
    if name == "440 nm":
        return wlen, np.asarray(env["intensity440"], dtype=float)
    if name == "560 nm":
        return wlen, np.asarray(env["intensity560"], dtype=float)
    if name == "Xenon":
        return wlen, np.asarray(env["xenon"]["photons"], dtype=float)
    if name == "Xenon (eye)":
        return wlen, np.asarray(env["xenon_eye"]["photons"], dtype=float)
    if name in app._custom_spectra:
        return app._custom_spectra[name]
    return None


def _impulse_monochromatic(
    wlen_grid: np.ndarray, peak_nm: float
) -> tuple[np.ndarray, np.ndarray]:
    """Single-bin impulse on ``wlen_grid`` at the sample nearest ``peak_nm`` (nm).

    Intensity is zero except ``1.0`` at that index. When used in the stimulus
    protocol, total photons scale with the trapezoid integral of this shape times
    the interval intensity the user sets in the Stimulus Builder.
    """
    wlen_grid = np.asarray(wlen_grid, dtype=float)
    lo, hi = float(wlen_grid[0]), float(wlen_grid[-1])
    if not (lo <= peak_nm <= hi):
        raise ValueError(
            f"Target wavelength {peak_nm} nm is outside the model grid "
            f"[{lo:.3f}, {hi:.3f}] nm."
        )
    idx = int(np.argmin(np.abs(wlen_grid - peak_nm)))
    intensity = np.zeros_like(wlen_grid, dtype=float)
    intensity[idx] = 1.0
    return wlen_grid, intensity


_MIX_NORM_NONE = "none"
_MIX_NORM_INTEGRAL = "integral"
_MIX_NORM_PEAK = "peak"


def _normalize_mix_profile(
    profile: np.ndarray,
    wlen: np.ndarray,
    mode: str,
    name: str,
) -> np.ndarray:
    """Normalize a resampled mix-row profile by peak, integral, or leave raw.

    Zero profiles (denominator exactly 0) are returned unchanged. Non-finite or
    negative denominators raise ``ValueError`` naming ``name``.
    """
    if mode == _MIX_NORM_NONE:
        return profile
    if mode == _MIX_NORM_INTEGRAL:
        denom = float(np.trapezoid(profile, wlen))
    elif mode == _MIX_NORM_PEAK:
        denom = float(np.max(profile))
    else:
        raise ValueError(f"Unknown mix normalization mode {mode!r}.")
    if denom == 0.0:
        return profile
    if not np.isfinite(denom) or denom < 0:
        raise ValueError(
            f"Cannot normalize spectrum {name!r}: "
            f"denominator must be finite and non-negative (got {denom})."
        )
    return profile / denom


def _mix_weighted_spectra(
    app: "ManuscriptSimApp",
    rows: list[tuple[str, float]],
    norm_mode: str = _MIX_NORM_PEAK,
) -> tuple[np.ndarray, np.ndarray]:
    """Linearly combine library spectra on ``wlenshared``.

    Each row is resampled, optionally normalized (peak or integrated photon
    count), then multiplied by its weight and summed.
    """
    env = app._env
    if env is None:
        raise ValueError("Spectra not loaded.")
    target_wlen = np.asarray(env["wlenshared"], dtype=float)
    names = [str(n).strip() for n, _ in rows]
    weights = [float(w) for _, w in rows]
    if not rows:
        raise ValueError("Add at least one spectrum row.")
    combined = np.zeros_like(target_wlen, dtype=float)
    for name, w in zip(names, weights):
        if w < 0:
            raise ValueError(f"Negative weight for spectrum {name!r}.")
        data = _library_spectrum_arrays(app, name)
        if data is None:
            raise ValueError(
                f"Spectrum {name!r} could not be resolved (environment not loaded?)."
            )
        src_wlen, src_int = data
        prof = _resample_spectrum(src_wlen, src_int, target_wlen)
        combined += w * _normalize_mix_profile(prof, target_wlen, norm_mode, name)
    return target_wlen, combined


def _spectrum_picker_entries(app: "ManuscriptSimApp") -> list[str]:
    """Names shown in spectrum pickers (built-ins minus hidden, then custom)."""
    entries: list[str] = []
    for name in _SPECTRUM_LIBRARY_ENTRIES:
        if name not in app._hidden_builtin_spectra:
            entries.append(name)
    entries.extend(sorted(app._custom_spectra.keys()))
    return entries


def _builtin_spectrum_for_protocol(name: str, env: dict) -> np.ndarray | None:
    """Return the photon-flux profile on ``env['wlenshared']`` for built-in names."""
    target_wlen = np.asarray(env["wlenshared"], dtype=float)
    if name == "Dark":
        return np.zeros_like(target_wlen, dtype=float)
    if name == "440 nm":
        return np.asarray(env["intensity440"], dtype=float)
    if name == "560 nm":
        return np.asarray(env["intensity560"], dtype=float)
    if name == "Xenon":
        return np.asarray(env["xenon"]["photons"], dtype=float)
    if name == "Xenon (eye)":
        return np.asarray(env["xenon_eye"]["photons"], dtype=float)
    return None


def _build_custom_protocol(
    spec: dict,
    env: dict,
    custom_spectra: dict[str, tuple[np.ndarray, np.ndarray]],
    rate: float = MANUSCRIPT_RATE,
) -> dict:
    """Assemble a ``predictMelanopsin``-ready protocol dict from a v2 spec.

    Walks ``spec["blocks"]`` and expands each interval block into the
    per-interval ``intensities`` / ``wlen`` / ``timings`` arrays the model
    consumes.
    """
    target_wlen = np.asarray(env["wlenshared"], dtype=float)
    blocks = spec["blocks"]

    intensities_rows: list[np.ndarray] = []
    wlen_rows: list[np.ndarray] = []
    timings: list[tuple[float, float]] = []

    t = 0.0
    for block in blocks:
        if block.get("type") == "interval":
            spectrum_name = str(block["spectrum_ref"]).strip()
            duration = float(block["duration"])
            target_intensity = float(block["intensity"])
            builtin = _builtin_spectrum_for_protocol(spectrum_name, env)
            if builtin is not None:
                scaled = _scale_spectrum_to_integrated_intensity(
                    target_wlen,
                    builtin,
                    target_intensity,
                    target_wlen=target_wlen,
                    ref_name=spectrum_name,
                )
            else:
                data = custom_spectra.get(spectrum_name)
                if data is None:
                    raise ValueError(
                        f"Spectrum {spectrum_name!r} is not a built-in or in the "
                        "spectrum library; cannot build protocol."
                    )
                src_wlen, src_int = data
                scaled = _scale_spectrum_to_integrated_intensity(
                    src_wlen,
                    src_int,
                    target_intensity,
                    target_wlen=target_wlen,
                    ref_name=spectrum_name,
                )
            intensities_rows.append(scaled)
            wlen_rows.append(target_wlen.copy())
            timings.append((t, t + duration))
            t += duration
        else:
            raise ValueError(f"Unknown block type: {block.get('type')!r}")

    timings_arr = np.array(timings, dtype=float).reshape(-1)
    intensities_arr = np.array(intensities_rows, dtype=float)
    wlen_arr = np.array(wlen_rows, dtype=float)
    # Use the larger of declared total vs accumulated interval end so trailing
    # dark is honored (binary float summation of many short intervals can
    # undershoot declared ti slightly).
    declared_ti = float(spec.get("total_duration", t))
    ti_val = max(declared_ti, t)
    return {
        "ti": ti_val,
        "timings": timings_arr,
        "intensities": intensities_arr,
        "wlen": wlen_arr,
        "rate": float(rate),
    }





def _format_autosave_path_for_display(p: Path) -> str:
    try:
        return str(p.resolve())
    except OSError:
        return str(p)


def _parse_stimulus_label_from_json(data: object) -> str | None:
    if not isinstance(data, dict):
        return None
    for key in ("label", "stimulus"):
        val = data.get(key)
        if isinstance(val, str):
            return val
    return None


def _iter_exportable_series(result: dict) -> tuple[np.ndarray, list[tuple[str, np.ndarray]]]:
    """Return ``xaxis`` and (name, 1D array) pairs matching ``xaxis`` length."""
    xaxis = np.asarray(result["xaxis"], dtype=float)
    n = int(xaxis.shape[0])
    pairs: list[tuple[str, np.ndarray]] = []
    for key in sorted(result.keys()):
        if key in ("params", "xaxis"):
            continue
        v = result[key]
        if isinstance(v, np.ndarray) and v.ndim == 1 and v.shape[0] == n:
            pairs.append((key, np.asarray(v, dtype=float)))
    return xaxis, pairs


def _export_result_csv(path: str, result: dict) -> None:
    xaxis, series = _iter_exportable_series(result)
    if not series:
        raise ValueError("No time-aligned arrays to export.")
    cols = [xaxis] + [a for _, a in series]
    names = ["time_s"] + [name for name, _ in series]
    data = np.column_stack(cols)
    header = ",".join(names)
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def _export_result_npz(path: str, result: dict, label: str) -> None:
    xaxis, series = _iter_exportable_series(result)
    if not series:
        raise ValueError("No time-aligned arrays to export.")
    payload: dict[str, np.ndarray] = {
        "stimulus_label": np.array(label, dtype=object),
        "time_s": xaxis,
    }
    for name, arr in series:
        payload[name] = arr
    np.savez_compressed(path, **payload)


def _load_saved_prediction_npz(path: str) -> dict:
    """Load exported prediction NPZ into normalized schema."""
    data = np.load(path, allow_pickle=True)
    if "time_s" not in data:
        raise ValueError("NPZ is missing required 'time_s' array.")
    time_s = np.asarray(data["time_s"], dtype=float).reshape(-1)
    if time_s.size == 0:
        raise ValueError("NPZ 'time_s' is empty.")

    label = Path(path).stem
    if "stimulus_label" in data:
        raw = data["stimulus_label"]
        try:
            if np.asarray(raw).shape == ():
                raw = np.asarray(raw).item()
        except Exception:
            pass
        if isinstance(raw, str) and raw.strip():
            label = raw.strip()

    series: dict[str, np.ndarray] = {}
    n = int(time_s.shape[0])
    for key in data.files:
        if key in ("time_s", "stimulus_label"):
            continue
        arr = np.asarray(data[key])
        if arr.ndim != 1 or int(arr.shape[0]) != n:
            continue
        if not np.issubdtype(arr.dtype, np.number):
            continue
        series[key] = np.asarray(arr, dtype=float)
    if not series:
        raise ValueError("NPZ has no numeric 1D time-aligned fields to plot.")
    return {"label": label, "time_s": time_s, "series": series, "source_path": path}


def _load_saved_prediction_csv(path: str) -> dict:
    """Load exported prediction CSV into normalized schema."""
    df = pd.read_csv(path)
    if "time_s" not in df.columns:
        raise ValueError("CSV is missing required 'time_s' column.")
    time_s = pd.to_numeric(df["time_s"], errors="coerce").to_numpy(dtype=float)
    if time_s.size == 0:
        raise ValueError("CSV 'time_s' is empty.")
    if np.any(~np.isfinite(time_s)):
        raise ValueError("CSV 'time_s' contains non-numeric values.")

    series: dict[str, np.ndarray] = {}
    n = int(time_s.shape[0])
    for col in df.columns:
        if col == "time_s":
            continue
        arr = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        if arr.ndim != 1 or int(arr.shape[0]) != n:
            continue
        if np.any(~np.isfinite(arr)):
            continue
        series[str(col)] = arr
    if not series:
        raise ValueError("CSV has no numeric 1D time-aligned fields to plot.")

    return {
        "label": Path(path).stem,
        "time_s": time_s,
        "series": series,
        "source_path": path,
    }


def _load_saved_prediction(path: str) -> dict:
    """Load an exported prediction file (.npz/.csv) into normalized schema."""
    lower = path.lower()
    if lower.endswith(".npz"):
        return _load_saved_prediction_npz(path)
    if lower.endswith(".csv"):
        return _load_saved_prediction_csv(path)
    raise ValueError(f"Unsupported file type for comparison: {Path(path).suffix}")


def _normalize_in_memory_result_for_compare(result: dict, label: str | None) -> dict:
    """Normalize an in-memory model result dict for comparison plotting."""
    xaxis, pairs = _iter_exportable_series(result)
    if not pairs:
        raise ValueError("Current run has no numeric 1D time-aligned fields to plot.")
    run_label = (label or "").strip() or "Current run"
    series = {name: np.asarray(arr, dtype=float) for name, arr in pairs}
    return {
        "label": run_label,
        "time_s": np.asarray(xaxis, dtype=float),
        "series": series,
        "source_path": "<current run>",
    }


def _coerce_model_dict_for_plot_model_vs_data(source: dict) -> dict:
    """Build a ``model`` dict for :func:`plotModelVsData` from live output or saved export."""
    if "time_s" in source and "series" in source:
        time_s = np.asarray(source["time_s"], dtype=float)
        m: dict = {"xaxis": time_s, "params": source.get("params", {})}
        for k, v in source["series"].items():
            m[k] = np.asarray(v, dtype=float)
    elif "xaxis" in source:
        m = dict(source)
    else:
        raise ValueError("Unrecognized run format (need xaxis or time_s+series).")

    if "lightIntensity" not in m:
        if "stimMonitor" in m:
            m["lightIntensity"] = np.asarray(m["stimMonitor"], dtype=float).ravel()
        else:
            raise ValueError(
                "This run has no lightIntensity or stimMonitor, which plotModelVsData "
                "needs for the light strip. Export data from a full simulation run "
                "(CSV/NPZ) that includes one of those series, or use Load current run."
            )
    return m


def _aligned_series_field_names(model: dict) -> list[str]:
    """Names of numeric 1D arrays aligned with ``model['xaxis']`` (for series picker)."""
    _, pairs = _iter_exportable_series(model)
    names = [name for name, _ in pairs]
    return sorted(names)


def _default_data_compare_series_name(model: dict) -> str:
    names = _aligned_series_field_names(model)
    if DEFAULT_DATA_COMPARE_SERIES in names:
        return DEFAULT_DATA_COMPARE_SERIES
    return names[0] if names else ""


# Normalized base names that identify a time column in a headered dataset.
_DATA_TIME_ALIASES = frozenset(
    {
        "time",
        "times",
        "t",
        "time_s",
        "timestamp",
        "sec",
        "secs",
        "second",
        "seconds",
        "elapsed",
        "elapsed time",
    }
)

# Unit token -> multiplier to convert values to seconds.
_TIME_UNIT_SCALES: dict[str, float] = {
    "s": 1.0,
    "sec": 1.0,
    "secs": 1.0,
    "second": 1.0,
    "seconds": 1.0,
    "ms": 1e-3,
    "msec": 1e-3,
    "msecs": 1e-3,
    "millisecond": 1e-3,
    "milliseconds": 1e-3,
    "us": 1e-6,
    "usec": 1e-6,
    "µs": 1e-6,
    "microsecond": 1e-6,
    "microseconds": 1e-6,
    "min": 60.0,
    "mins": 60.0,
    "minute": 60.0,
    "minutes": 60.0,
    "h": 3600.0,
    "hr": 3600.0,
    "hrs": 3600.0,
    "hour": 3600.0,
    "hours": 3600.0,
}


def _parse_time_unit_scale(label: str) -> tuple[float, str]:
    """Infer a seconds-multiplier and display unit from a column header.

    Recognizes trailing parenthesized/bracketed units (``Time (ms)``,
    ``Time [s]``) and underscore suffixes (``time_ms``, ``time_s``). When no
    unit is recognized the values are treated as seconds.
    """
    text = str(label).strip()
    lower = text.lower()

    # Trailing (unit) or [unit]
    m = re.search(r"[\[(]\s*([a-zµμ]+)\s*[\])]\s*$", lower)
    if m:
        token = m.group(1).replace("μ", "µ")
        if token in _TIME_UNIT_SCALES:
            return _TIME_UNIT_SCALES[token], token

    # Underscore / space suffix: time_ms, time_s, elapsed_min
    m = re.search(r"[_\s]([a-zµμ]+)\s*$", lower)
    if m:
        token = m.group(1).replace("μ", "µ")
        if token in _TIME_UNIT_SCALES:
            return _TIME_UNIT_SCALES[token], token

    return 1.0, "s"


def _normalize_time_column_base(label: str) -> str:
    """Strip unit suffixes from a column label for time-alias matching."""
    text = str(label).strip().lower()
    # Drop trailing (unit) / [unit]
    text = re.sub(r"\s*[\[(]\s*[a-zµμ]+\s*[\])]\s*$", "", text)
    # Drop trailing _unit / space-unit when the token is a known time unit
    m = re.search(r"^(.*?)[_\s]([a-zµμ]+)\s*$", text)
    if m and m.group(2).replace("μ", "µ") in _TIME_UNIT_SCALES:
        text = m.group(1).strip()
    return text


def _is_time_column_name(label: str) -> bool:
    """True if ``label`` matches a known time-column alias (unit-tolerant)."""
    base = _normalize_time_column_base(label)
    if base in _DATA_TIME_ALIASES:
        return True
    # Also accept the raw lowercased label (e.g. "time_s" is itself an alias)
    return str(label).strip().lower() in _DATA_TIME_ALIASES


def _trace_columns_look_headerless(columns) -> bool:
    """True if column labels look like data values rather than names.

    Unlike ``_columns_look_headerless``, this accepts single-column files
    (legacy manuscript recordings are one bare numeric column).
    """
    cols = list(columns)
    if len(cols) == 0:
        return False
    labels = [str(c).strip() for c in cols]
    if all(lab.lower().startswith("unnamed") for lab in labels):
        return True
    try:
        for lab in labels:
            float(lab)
        return True
    except ValueError:
        return False


def _load_trace_dataframe(path: str, *, header=0) -> "pd.DataFrame":
    """Load a dataset spreadsheet (CSV/TSV/TXT/Excel) for the Data Comparator."""
    lower = path.lower()
    if lower.endswith(".xlsx") or lower.endswith(".xls"):
        return pd.read_excel(path, header=header)
    sep = _infer_spectrum_text_sep(path)
    kwargs: dict = {"header": header, "sep": sep}
    if sep == r"\s+":
        kwargs["engine"] = "python"
    return pd.read_csv(path, **kwargs)


def _read_data_trace(path: str) -> dict:
    """Load an experimental trace for the Data Comparator.

    Returns a dict with:
      - ``columns``: ordered ``{name: np.ndarray}`` of numeric non-time columns
        (NaNs retained so time/value rows stay aligned until selection time)
      - ``time_s``: time in seconds (zero-based), or ``None`` if no time column
      - ``time_column``: original time column name, or ``None``
      - ``time_unit_label``: display unit token (e.g. ``"ms"``), or ``None``
      - ``time_offset_s``: original ``t[0]`` in seconds before zero-basing (0.0
        when no time column)
      - ``headerless``: whether the file was treated as headerless
    """
    df = _load_trace_dataframe(path, header=0)
    headerless = False
    if _trace_columns_look_headerless(df.columns):
        df = _load_trace_dataframe(path, header=None)
        headerless = True

    if df.shape[1] == 0 or df.shape[0] == 0:
        raise ValueError("Dataset file is empty.")

    time_col = None
    time_s = None
    time_unit_label = None
    time_offset_s = 0.0

    if not headerless:
        for col in df.columns:
            if _is_time_column_name(str(col)):
                time_col = col
                break

    if time_col is not None:
        scale, time_unit_label = _parse_time_unit_scale(str(time_col))
        raw_t = pd.to_numeric(df[time_col], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(raw_t)
        if not np.any(finite):
            raise ValueError(
                f"Time column {str(time_col)!r} has no numeric values."
            )
        t_scaled = raw_t * scale
        # Zero-base using the first finite sample so the trace starts at 0 s
        # (matching the model xaxis convention).
        first_finite = float(t_scaled[np.argmax(finite)])
        time_offset_s = first_finite
        time_s = t_scaled - first_finite

    columns: dict[str, np.ndarray] = {}
    for col in df.columns:
        if time_col is not None and col == time_col:
            continue
        if headerless and df.shape[1] == 1:
            name = "value"
        elif headerless:
            name = f"column_{int(col) + 1}"
        else:
            name = str(col)
        arr = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        columns[name] = arr

    if not columns:
        raise ValueError("No numeric data columns found (only a time column?).")

    # Ensure at least one column has some finite values
    if not any(np.any(np.isfinite(v)) for v in columns.values()):
        raise ValueError("No numeric values found in data columns.")

    return {
        "columns": columns,
        "time_s": time_s,
        "time_column": str(time_col) if time_col is not None else None,
        "time_unit_label": time_unit_label,
        "time_offset_s": float(time_offset_s),
        "headerless": headerless,
    }


def _build_stimulus(label: str, env: dict):
    """Return protocol dict for ``predictMelanopsin`` (matches tutorial)."""
    wlenshared = env["wlenshared"]
    intensity440 = env["intensity440"]
    intensity560 = env["intensity560"]
    xenon = env["xenon"]
    xenon_eye = env["xenon_eye"]
    xph = np.asarray(xenon["photons"], dtype=float)
    xeph = np.asarray(xenon_eye["photons"], dtype=float)

    if label == "440-440":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 440",
            rate=MANUSCRIPT_RATE,
            specwlen=wlenshared,
            specint=intensity440,
        )
    if label == "560-560":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 560",
            rate=MANUSCRIPT_RATE,
            specwlen=wlenshared,
            specint=intensity560,
        )
    if label == "440-560":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 440 560",
            rate=MANUSCRIPT_RATE,
            specwlen=[wlenshared, wlenshared],
            specint=[intensity440, intensity560],
        )
    if label == "560-440":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 560 440",
            rate=MANUSCRIPT_RATE,
            specwlen=[wlenshared, wlenshared],
            specint=[intensity560, intensity440],
        )
    if label == "Xenon":
        return setupModelRun(
            "intensity response",
            "xenon",
            rate=MANUSCRIPT_RATE,
            specwlen=wlenshared,
            specint=xph,
        )
    if label == "Xenon with ocular filtering":
        return setupModelRun(
            "intensity response",
            "xenon eye",
            rate=MANUSCRIPT_RATE,
            specwlen=wlenshared,
            specint=xeph,
        )
    raise ValueError(f"Unknown stimulus label: {label!r}")


SPECTRAL_SENSITIVITY_VALUES = ("Custom", "Govardovskii")
# Left-column labels in the configure dialog (key is still the model dict key)
_PARAM_DISPLAY_LABELS = {
    "rate": "rate (s)",
}
CONFIG_SECTION_STARTS = {
    "activity feedback": "Toggle features",
    "Rphi": "Quantum efficiencies",
    "Reps": "Absorption coefficients",
    "spectral sensitivity": "Nomogram function",
    "R_true_lambda_max": "R state nomogram",
    "M_true_lambda_max": "M state nomogram",
    "E_true_lambda_max": "E state nomogram",
    "propME": "Transition probability",
    "konM": "Rate constants",
    "konMG": "Leaky integrators for gain control",
    "kbleachR": "Bleaching",
    "kthermM": "Thermal transitions of M",
    "m deplete thresh": "M to M' threshold behavior",
    "MprimeGain": "Signaling gain",
    "k_mg": "M'-driven feedback nonlinearity",
    "k": "Activity-driven feedback nonlinearity",
}


def _format_param_for_entry(value: object) -> str:
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, numbers.Real):
        return str(float(value))
    return str(value)


def _parse_param_value(key: str, text: str, template: object) -> object:
    s = text.strip()
    if s.startswith("np.") and s.endswith(")") and "(" in s:
        s = s[s.find("(") + 1 : -1].strip()
    if isinstance(template, bool):
        low = s.lower()
        if low in ("true", "1", "yes"):
            return True
        if low in ("false", "0", "no"):
            return False
        raise ValueError(f'expected True or False, got {text!r}')
    if isinstance(template, str):
        return s
    if isinstance(template, numbers.Integral):
        return int(float(s))
    if isinstance(template, numbers.Real):
        return float(s)
    raise TypeError(f"unsupported template type for {key!r}: {type(template)!r}")


class ModelConfigDialog(tk.Toplevel):
    """Scrollable editor for ``predictMelanopsin`` default parameters."""

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Model parameters")
        self._parent = parent
        self.transient(parent)
        self._wheel_bindtags_active = False
        self._template = get_predict_melanopsin_defaults()
        self._param_templates: list[tuple[str, object]] = [
            (k, v) for k, v in self._template.items() if k != "ti"
        ]
        self._entries: dict[str, ttk.Entry | ttk.Combobox] = {}
        self._section_font = tkfont.nametofont("TkDefaultFont", root=self).copy()
        self._section_font.configure(weight="bold")

        outer = ttk.Frame(self, padding=8)

        canvas = tk.Canvas(outer, highlightthickness=0)
        vsb = ttk.Scrollbar(outer, orient=tk.VERTICAL, command=canvas.yview)
        inner = ttk.Frame(canvas)

        def _on_inner_configure(_event=None) -> None:
            canvas.configure(scrollregion=canvas.bbox("all"))

        inner.bind("<Configure>", _on_inner_configure)

        win_id = canvas.create_window((0, 0), window=inner, anchor=tk.NW)

        def _on_canvas_configure(event: tk.Event) -> None:
            canvas.itemconfigure(win_id, width=event.width)

        canvas.bind("<Configure>", _on_canvas_configure)

        canvas.configure(yscrollcommand=vsb.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)

        for key, tmpl in self._param_templates:
            section = CONFIG_SECTION_STARTS.get(key)
            if section is not None:
                ttk.Separator(inner, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=(8, 4))
                ttk.Label(inner, text=section, font=self._section_font).pack(
                    fill=tk.X, pady=(2, 2), anchor=tk.W
                )
            cur = parent._model_config.get(key, tmpl)
            row = ttk.Frame(inner)
            row.pack(fill=tk.X, pady=2)
            ttk.Label(
                row, text=_PARAM_DISPLAY_LABELS.get(key, key), width=28
            ).pack(side=tk.LEFT, padx=(0, 8))
            if key == "spectral sensitivity":
                cb = ttk.Combobox(
                    row,
                    values=SPECTRAL_SENSITIVITY_VALUES,
                    state="readonly",
                    width=36,
                )
                val = str(cur)
                if val not in SPECTRAL_SENSITIVITY_VALUES:
                    val = SPECTRAL_SENSITIVITY_VALUES[0]
                cb.set(val)
                cb.bind("<MouseWheel>", lambda _e: "break")
                cb.bind("<Button-4>", lambda _e: "break")
                cb.bind("<Button-5>", lambda _e: "break")
                cb.pack(side=tk.LEFT, fill=tk.X, expand=True)
                self._entries[key] = cb
            else:
                ent = ttk.Entry(row, width=40)
                ent.insert(0, _format_param_for_entry(cur))
                ent.pack(side=tk.LEFT, fill=tk.X, expand=True)
                self._entries[key] = ent

        def _wheel(event: tk.Event) -> str:
            if getattr(event, "num", None) == 4:
                canvas.yview_scroll(-1, "units")
                return "break"
            if getattr(event, "num", None) == 5:
                canvas.yview_scroll(1, "units")
                return "break"
            if event.delta:
                canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
            return "break"

        self._wheel_handler = _wheel
        self._bind_global_wheel()

        btn_row = ttk.Frame(self, padding=(8, 0, 8, 8))
        ttk.Button(
            btn_row,
            text="Save Parameters",
            command=self._on_save_parameters,
        ).pack(side=tk.LEFT)
        ttk.Button(
            btn_row,
            text="Load Preset",
            command=self._on_load_preset,
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(btn_row, text="Update", command=self._on_save).pack(side=tk.RIGHT)
        ttk.Button(
            btn_row,
            text="Reset to defaults",
            command=self._on_reset_defaults,
        ).pack(side=tk.RIGHT, padx=(0, 8))
        btn_row.pack(side=tk.BOTTOM, fill=tk.X)
        outer.pack(fill=tk.BOTH, expand=True)

        self.minsize(520, 400)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _bind_global_wheel(self, _event=None) -> None:
        if self._wheel_bindtags_active:
            return
        self.bind_all("<MouseWheel>", self._wheel_handler, add="+")
        self.bind_all("<Button-4>", self._wheel_handler, add="+")
        self.bind_all("<Button-5>", self._wheel_handler, add="+")
        self._wheel_bindtags_active = True

    def _unbind_global_wheel(self, _event=None) -> None:
        if not self._wheel_bindtags_active:
            return
        self.unbind_all("<MouseWheel>")
        self.unbind_all("<Button-4>")
        self.unbind_all("<Button-5>")
        self._wheel_bindtags_active = False

    def _on_close(self) -> None:
        self._unbind_global_wheel()
        self._parent._config_dialog = None
        self.destroy()

    def _on_reset_defaults(self) -> None:
        defaults = get_predict_melanopsin_defaults()
        self._apply_params_to_entries(defaults)

    def _collect_params_from_entries(self) -> dict | None:
        """Parse the form into a model-config dict, or ``None`` on error."""
        out: dict = {}
        for key, tmpl in self._param_templates:
            w = self._entries[key]
            if key == "spectral sensitivity":
                raw = w.get()
                if raw not in SPECTRAL_SENSITIVITY_VALUES:
                    messagebox.showerror(
                        "Model parameters",
                        f'spectral sensitivity must be one of {SPECTRAL_SENSITIVITY_VALUES!r}.',
                        parent=self,
                    )
                    return None
                out[key] = raw
            else:
                assert isinstance(w, ttk.Entry)
                raw = w.get()
                try:
                    parsed = _parse_param_value(key, raw, tmpl)
                except (ValueError, TypeError, OverflowError) as exc:
                    messagebox.showerror(
                        "Model parameters",
                        f"Could not parse {key!r} = {raw!r}:\n{exc}",
                        parent=self,
                    )
                    return None
                out[key] = parsed
        return out

    def _apply_params_to_entries(self, params: dict) -> bool:
        """Fill widgets from a parameter dict (unknown keys ignored).

        Returns ``True`` on success, ``False`` if a value could not be applied.
        """
        for key, tmpl in self._param_templates:
            if key not in params:
                continue
            w = self._entries[key]
            value = params[key]
            if key == "spectral sensitivity":
                v = str(value)
                if v not in SPECTRAL_SENSITIVITY_VALUES:
                    messagebox.showerror(
                        "Model parameters",
                        f'spectral sensitivity must be one of {SPECTRAL_SENSITIVITY_VALUES!r} '
                        f"(got {v!r}).",
                        parent=self,
                    )
                    return False
                w.set(v)
                continue
            assert isinstance(w, ttk.Entry)
            try:
                parsed = _parse_param_value(
                    key, _format_param_for_entry(value), tmpl
                )
            except (ValueError, TypeError, OverflowError) as exc:
                messagebox.showerror(
                    "Model parameters",
                    f"Could not apply {key!r} = {value!r}:\n{exc}",
                    parent=self,
                )
                return False
            w.delete(0, tk.END)
            w.insert(0, _format_param_for_entry(parsed))
        return True

    def _on_save_parameters(self) -> None:
        out = self._collect_params_from_entries()
        if out is None:
            return
        initial_dir = library_storage.configurations_dir()
        try:
            initial_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        path = filedialog.asksaveasfilename(
            parent=self,
            title="Save model parameters",
            defaultextension=".json",
            initialdir=str(initial_dir),
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(out, f, indent=2, sort_keys=False)
                f.write("\n")
        except OSError as exc:
            messagebox.showerror(
                "Model parameters",
                f"Could not write parameters:\n{exc}",
                parent=self,
            )
            return
        self._parent._status.config(
            text=f"Model parameters written to {path}."
        )
        messagebox.showinfo(
            "Model parameters",
            f"Parameters saved to:\n{path}",
            parent=self,
        )

    def _on_load_preset(self) -> None:
        initial_dir = library_storage.configurations_dir()
        try:
            initial_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
        path = filedialog.askopenfilename(
            parent=self,
            title="Load model parameter preset",
            initialdir=str(initial_dir) if initial_dir.is_dir() else None,
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            messagebox.showerror(
                "Model parameters",
                f"Could not read preset:\n{exc}",
                parent=self,
            )
            return
        if not isinstance(data, dict):
            messagebox.showerror(
                "Model parameters",
                "Preset JSON must be an object mapping parameter names to values.",
                parent=self,
            )
            return
        known = {key for key, _ in self._param_templates}
        usable = {k: v for k, v in data.items() if k in known}
        if not usable:
            messagebox.showerror(
                "Model parameters",
                "Preset did not contain any recognized model parameters.",
                parent=self,
            )
            return
        if not self._apply_params_to_entries(usable):
            return
        self._parent._status.config(
            text=f"Loaded model parameter preset from {path}."
        )
        messagebox.showinfo(
            "Model parameters",
            (
                f"Loaded {len(usable)} parameter(s) from:\n{path}\n\n"
                "Click Update to apply them to the next model run."
            ),
            parent=self,
        )

    def _on_save(self) -> None:
        out = self._collect_params_from_entries()
        if out is None:
            return
        self._parent._model_config = out
        self._parent._status.config(text="Model parameters updated.")
        messagebox.showinfo(
            "Model parameters",
            "Updated. The next run will use these values.",
            parent=self,
        )
        self._on_close()


class CustomStimulusDialog(tk.Toplevel):
    """Placeholder window for future custom stimulus creation."""

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Custom Stimulus Builder")
        self._parent = parent
        self.transient(parent)

        outer = ttk.Frame(self, padding=12)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(
            outer,
            text="Custom stimulus creation will be added here.",
        ).pack(anchor=tk.W)

        btn_row = ttk.Frame(outer, padding=(0, 10, 0, 0))
        btn_row.pack(fill=tk.X)
        ttk.Button(
            btn_row, text="Load stimulus…", command=self._on_load_stimulus
        ).pack(side=tk.LEFT)
        ttk.Button(btn_row, text="Close", command=self._on_close).pack(side=tk.RIGHT)

        self.minsize(420, 180)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _on_close(self) -> None:
        self._parent._custom_stim_dialog = None
        self.destroy()

    def _on_load_stimulus(self) -> None:
        self._parent._on_file_load_stimulus()


class _SpectrumPickerPopup(tk.Toplevel):
    """Small modal popup listing all available spectra for selection."""

    def __init__(
        self,
        parent: tk.Toplevel,
        app: "ManuscriptSimApp",
        current: str,
        on_select,
    ) -> None:
        super().__init__(parent)
        self.title("Select Spectrum")
        self.transient(parent)
        self.resizable(False, False)
        self._on_select = on_select

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        list_frame = ttk.Frame(outer)
        list_frame.pack(fill=tk.BOTH, expand=True)
        self._listbox = tk.Listbox(list_frame, height=10, exportselection=False, width=30)
        self._listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self._listbox.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self._listbox.configure(yscrollcommand=sb.set)

        entries: list[str] = []
        for name in _SPECTRUM_LIBRARY_ENTRIES:
            if name not in app._hidden_builtin_spectra:
                entries.append(name)
        for name in app._custom_spectra:
            entries.append(name)
        for name in entries:
            self._listbox.insert(tk.END, name)
        if current in entries:
            idx = entries.index(current)
            self._listbox.selection_set(idx)
            self._listbox.see(idx)

        self._listbox.bind("<Double-Button-1>", self._on_choose)
        self._listbox.bind("<Return>", self._on_choose)

        btn_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        btn_row.pack(fill=tk.X)
        ttk.Button(btn_row, text="Select", command=self._on_choose).pack(side=tk.LEFT)
        ttk.Button(btn_row, text="Cancel", command=self.destroy).pack(side=tk.RIGHT)

        self.grab_set()

    def _on_choose(self, _event=None) -> None:
        sel = self._listbox.curselection()
        if not sel:
            return
        self._on_select(str(self._listbox.get(sel[0])))
        self.destroy()


class StimulusBuilderDialog(tk.Toplevel):
    """Spreadsheet-style interval editor for building custom stimuli."""

    _NONE_LABEL = "-- none --"
    _SELECT_BG = "#9ec9f5"
    # Skip per-interval monitor/grid refresh when loading more than this many intervals.
    _DEFER_REFRESH_INTERVAL_THRESHOLD = 200
    # Win32 Tk canvas coords overflow near 32767 px; keep each embedded frame narrower.
    _MAX_CANVAS_CHUNK_WIDTH = 28000

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Stimulus Builder")
        self._parent = parent
        self.resizable(True, True)
        self.minsize(700, 340)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        # list of dicts, one per interval column
        self._intervals: list[dict] = []
        self._selected_col: int | None = None
        self._intensity_log_mode = tk.BooleanVar(master=self, value=True)
        self._compact_summary_window: dict | None = None
        self._layout_syncing = False

        # read the platform default widget background once
        _tmp = tk.Label(self)
        self._default_cell_bg: str = str(_tmp.cget("bg"))
        self._default_cell_fg: str = str(_tmp.cget("fg"))
        _tmp.destroy()

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        # ── top button bar ───────────────────────────────────────────
        top_bar = ttk.Frame(outer)
        top_bar.pack(fill=tk.X, pady=(0, 8))
        ttk.Button(
            top_bar, text="Add...", command=self._add_interval
        ).pack(side=tk.LEFT)
        ttk.Button(
            top_bar, text="Delete...", command=self._delete_interval
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(
            top_bar, text="Save", command=self._on_save_stimulus
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(
            top_bar, text="Load from library...", command=self._on_load_from_library
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(
            top_bar, text="Load manuscript...", command=self._on_load_manuscript_stimulus
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Checkbutton(
            top_bar,
            text="Log10 intensity input",
            variable=self._intensity_log_mode,
            command=self._on_intensity_mode_toggle,
        ).pack(side=tk.LEFT, padx=(12, 0))
        ttk.Label(top_bar, text="Duration (s):").pack(side=tk.LEFT, padx=(16, 4))
        self._duration_var = tk.StringVar(master=self, value="")
        self._duration_entry = ttk.Entry(
            top_bar, textvariable=self._duration_var, width=12
        )
        self._duration_entry.pack(side=tk.LEFT)
        self._duration_entry.bind("<Return>", self._on_declared_duration_commit)
        self._duration_entry.bind("<FocusOut>", self._on_declared_duration_commit)

        # ── main content: full-width grid on top, monitor+preview row below ────
        content = ttk.Frame(outer)
        content.pack(fill=tk.BOTH, expand=True)

        # top — horizontally scrollable canvas holding the interval grid
        # expand=False so leftover height goes to the monitor/preview row.
        grid_container = ttk.Frame(content)
        grid_container.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        self._h_scroll = ttk.Scrollbar(grid_container, orient=tk.HORIZONTAL)
        self._h_scroll.pack(side=tk.BOTTOM, fill=tk.X)
        self._v_scroll = ttk.Scrollbar(grid_container, orient=tk.VERTICAL)
        self._v_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        grid_pane = ttk.Frame(grid_container)
        grid_pane.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self._row_header_frame = tk.Frame(
            grid_pane, bg=self._default_cell_bg
        )
        self._row_header_frame.pack(side=tk.LEFT, fill=tk.Y)
        self._row_header_frame.pack_propagate(False)

        self._canvas = tk.Canvas(
            grid_pane,
            highlightthickness=0,
            bg=self._default_cell_bg,
            xscrollcommand=self._h_scroll.set,
            yscrollcommand=self._v_scroll.set,
        )
        self._canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._h_scroll.config(command=self._canvas.xview)
        self._v_scroll.config(command=self._canvas.yview)

        self._grid_chunks: list[dict] = []
        self._measured_col_width: int | None = None
        self._measured_col_pitch: int | None = None
        self._canvas.bind("<Configure>", self._on_canvas_configure)

        # separator
        ttk.Separator(content, orient=tk.HORIZONTAL).pack(
            side=tk.TOP, fill=tk.X, pady=(4, 4)
        )

        # bottom row — light monitor (left, 3/4) + spectrum preview (right, 1/4)
        bottom_row = ttk.Frame(content, height=220)
        bottom_row.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        bottom_row.pack_propagate(True)
        self._bottom_row = bottom_row

        monitor_frame = ttk.Frame(bottom_row, padding=(0, 0, 0, 0))
        monitor_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        ttk.Separator(bottom_row, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=(4, 4)
        )
        preview_frame = ttk.Frame(bottom_row, padding=(0, 0, 0, 0), width=320)
        preview_frame.pack(side=tk.LEFT, fill=tk.BOTH)
        preview_frame.pack_propagate(False)

        self._monitor_fig = plt.Figure(figsize=(9, 3))
        self._monitor_ax = self._monitor_fig.add_subplot(111)
        self._monitor_canvas = FigureCanvasTkAgg(self._monitor_fig, master=monitor_frame)
        self._monitor_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._monitor_interval_bounds: list[tuple[float, float, int]] = []
        self._monitor_xlim_full: tuple[float, float] | None = None
        self._monitor_last_hover_t: float | None = None
        self._monitor_selection_patch = None
        self._monitor_hover_col: int | None = None
        self._monitor_motion_cid = self._monitor_canvas.mpl_connect(
            "motion_notify_event", self._on_monitor_motion
        )
        self._monitor_leave_cid = self._monitor_canvas.mpl_connect(
            "axes_leave_event", self._on_monitor_leave
        )
        self._monitor_scroll_cid = self._monitor_canvas.mpl_connect(
            "scroll_event", self._on_monitor_scroll
        )
        self._monitor_click_cid = self._monitor_canvas.mpl_connect(
            "button_press_event", self._on_monitor_click
        )
        monitor_widget = self._monitor_canvas.get_tk_widget()
        monitor_widget.bind("<MouseWheel>", self._on_monitor_tk_wheel, add="+")
        monitor_widget.bind("<Button-4>", self._on_monitor_tk_wheel, add="+")
        monitor_widget.bind("<Button-5>", self._on_monitor_tk_wheel, add="+")
        monitor_widget.bind("<Double-Button-1>", self._on_monitor_reset_zoom, add="+")

        self._preview_fig = plt.Figure(figsize=(3, 3), tight_layout=True)
        self._preview_ax = self._preview_fig.add_subplot(111)
        self._preview_canvas = FigureCanvasTkAgg(self._preview_fig, master=preview_frame)
        self._preview_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._render_light_monitor()
        self._render_preview(None)

        # row-label column (left side of the grid)
        self._build_row_labels()

        # Left/Right navigate intervals when not editing a text field.
        self.bind("<Left>", self._on_interval_nav_left)
        self.bind("<Right>", self._on_interval_nav_right)

        # seed with the first interval
        self._stimulus_dirty = False
        self._tracking_edits = False
        self._add_interval()
        self._sync_duration_to_parent()
        self._loaded_source_name: str | None = None
        self._tracking_edits = True

    def _mark_stimulus_dirty(self) -> None:
        """Flag the builder as having unsaved edits (ignored during load/init)."""
        if self._tracking_edits:
            self._stimulus_dirty = True

    def _mark_stimulus_clean(self) -> None:
        """Clear the unsaved-edits flag after save or replace-load."""
        self._stimulus_dirty = False

    def _has_unsaved_stimulus_changes(self) -> bool:
        return bool(self._stimulus_dirty)

    # ── grid layout helpers ──────────────────────────────────────────

    def _column_width_px(self) -> int:
        if self._measured_col_width and self._measured_col_width > 0:
            return self._measured_col_width
        return 140

    def _column_pitch_px(self) -> int:
        """Horizontal space per column including ``padx`` between cells."""
        if self._measured_col_pitch and self._measured_col_pitch > 0:
            return self._measured_col_pitch
        return self._column_width_px() + 1

    def _max_cols_per_chunk(self) -> int:
        return max(1, self._MAX_CANVAS_CHUNK_WIDTH // self._column_width_px())

    def _chunk_index_for_interval(self, interval_index: int) -> tuple[int, int]:
        """Return ``(chunk_index, local_column)`` for a 0-based interval index."""
        per_chunk = self._max_cols_per_chunk()
        return interval_index // per_chunk, (interval_index % per_chunk) + 1

    def _chunk_x_offset(self, chunk_index: int) -> int:
        """Pixel x-offset for a chunk (do not use ``winfo_width``; it is stale until laid out)."""
        return chunk_index * self._max_cols_per_chunk() * self._column_pitch_px()

    def _ensure_grid_chunk(self, chunk_index: int) -> dict:
        while len(self._grid_chunks) <= chunk_index:
            idx = len(self._grid_chunks)
            x_offset = self._chunk_x_offset(idx)
            frame = tk.Frame(self._canvas, bg=self._default_cell_bg)
            window_id = self._canvas.create_window(
                x_offset, 0, window=frame, anchor=tk.NW
            )
            chunk = {"frame": frame, "window_id": window_id, "x_offset": x_offset}
            frame.bind("<Configure>", self._on_grid_configure)
            self._grid_chunks.append(chunk)
            self._apply_grid_vertical_size()
        return self._grid_chunks[chunk_index]

    def _teardown_interval_widgets(self) -> None:
        for entry in self._intervals:
            for w in entry.get("widgets", []):
                w.destroy()
            header = entry.get("header_label")
            if header is not None:
                header.destroy()
            entry["widgets"] = []
            entry["header_label"] = None
        for chunk in self._grid_chunks:
            self._canvas.delete(chunk["window_id"])
        self._grid_chunks.clear()
        self._hide_compact_grid_summary()

    def _interval_is_compact(self, interval: dict) -> bool:
        return bool(interval.get("compact"))

    def _grid_shows_compact_summary(self) -> bool:
        return bool(self._intervals) and all(
            self._interval_is_compact(iv) for iv in self._intervals
        )

    def _show_compact_grid_summary(self, n_intervals: int, source_name: str | None) -> None:
        """Replace the scrollable column grid with a compact-load summary panel."""
        self._hide_compact_grid_summary()
        frame = tk.Frame(self._canvas, bg=self._default_cell_bg)

        manual_compact = [
            iv for iv in self._intervals if self._interval_is_compact(iv)
        ]

        lines: list[str] = []
        if source_name:
            lines.append(f"Stimulus: {source_name}")
        if manual_compact:
            lines.append(
                f"{len(manual_compact):,} manual intervals loaded in compact view."
            )
            lines.append(
                "Load again and choose \"Show interval columns\" to edit each column."
            )
        lines.append("Light monitor and Save use every interval.")
        tk.Label(
            frame,
            text="\n".join(lines),
            justify=tk.LEFT,
            anchor=tk.NW,
            padx=12,
            pady=12,
            bg=self._default_cell_bg,
            fg=self._default_cell_fg,
        ).pack(anchor=tk.NW)
        window_id = self._canvas.create_window(0, 0, window=frame, anchor=tk.NW)
        self._compact_summary_window = {"frame": frame, "window_id": window_id}
        self._on_grid_configure()

    def _hide_compact_grid_summary(self) -> None:
        if self._compact_summary_window is None:
            return
        self._canvas.delete(self._compact_summary_window["window_id"])
        self._compact_summary_window["frame"].destroy()
        self._compact_summary_window = None

    def _update_compact_grid_display(self) -> None:
        if self._grid_shows_compact_summary():
            self._show_compact_grid_summary(len(self._intervals), self._loaded_source_name)
        else:
            self._hide_compact_grid_summary()

    def _sync_grid_scroll_background(self) -> None:
        """Size the scroll region to span all grid chunk windows."""
        if self._compact_summary_window is not None:
            self._canvas.update_idletasks()
            bbox = self._canvas.bbox(self._compact_summary_window["window_id"])
            canvas_h = self._canvas.winfo_height()
            if bbox:
                x2, y2 = bbox[2], max(bbox[3], canvas_h if canvas_h > 0 else 0)
            elif canvas_h > 0:
                x2, y2 = 400, canvas_h
            else:
                x2, y2 = 400, 200
            self._canvas.configure(scrollregion=(0, 0, x2, y2))
            return
        if not self._grid_chunks:
            return
        per = self._max_cols_per_chunk()
        pitch = self._column_pitch_px()
        n_intervals = len(self._intervals)
        n_chunks = len(self._grid_chunks)
        cols_in_last = n_intervals - (n_chunks - 1) * per if n_intervals else 0
        if cols_in_last <= 0:
            cols_in_last = per
        x2 = (n_chunks - 1) * per * pitch + cols_in_last * pitch

        for i, chunk in enumerate(self._grid_chunks):
            x = self._chunk_x_offset(i)
            if chunk["x_offset"] != x:
                chunk["x_offset"] = x
                self._canvas.coords(chunk["window_id"], x, 0)

        self._canvas.update_idletasks()
        y2 = 0
        for chunk in self._grid_chunks:
            bbox = self._canvas.bbox(chunk["window_id"])
            if bbox:
                x2 = max(x2, bbox[2])
                y2 = max(y2, bbox[3])
        canvas_h = self._canvas.winfo_height()
        if canvas_h > 0 and y2 < canvas_h:
            y2 = canvas_h
        self._canvas.configure(scrollregion=(0, 0, x2, y2))

    def _apply_grid_vertical_size(self, height: int | None = None) -> None:
        """Size label column and canvas chunk windows to natural content height.

        Rows are not stretched to fill the canvas; leftover window height goes to
        the light monitor / spectrum preview. ``height`` is ignored except for
        compact-summary mode, where the summary panel may fill the canvas.

        Callers that need reentrancy protection should set ``_layout_syncing``
        around this call. Updates are skipped when sizes are already correct.
        """
        self.update_idletasks()
        content_h = 0
        for row in range(4):
            info = self._row_header_frame.grid_rowconfigure(row)
            content_h += int(info.get("minsize") or 0)
        if content_h <= 1:
            for child in list(self._row_header_frame.winfo_children())[:4]:
                content_h += max(int(child.winfo_reqheight()), 1)
        content_h = max(content_h, 1)

        if self._compact_summary_window is not None:
            canvas_h = max(int(self._canvas.winfo_height()), 1)
            h = max(
                int(height if height is not None else canvas_h),
                content_h,
                120,
            )
            cur = self._canvas.itemcget(
                self._compact_summary_window["window_id"], "height"
            )
            if str(cur) != str(h):
                self._canvas.itemconfigure(
                    self._compact_summary_window["window_id"], height=h
                )
            if int(self._row_header_frame.cget("height") or 0) != h:
                self._row_header_frame.configure(height=h)
            return

        h = content_h
        for chunk in self._grid_chunks:
            cur = self._canvas.itemcget(chunk["window_id"], "height")
            if str(cur) != str(h):
                self._canvas.itemconfigure(chunk["window_id"], height=h)
        if int(self._row_header_frame.cget("height") or 0) != h:
            self._row_header_frame.configure(height=h)
        # Keep the packed canvas at content height so expand=False grid stays compact.
        # Only set when changed to avoid a Configure feedback loop.
        if int(self._canvas.cget("height") or 0) != h:
            self._canvas.configure(height=h)

    def _on_grid_configure(self, _event=None) -> None:
        if self._layout_syncing:
            return
        self._layout_syncing = True
        try:
            self._sync_row_heights()
            self._apply_grid_vertical_size()
            self._sync_grid_scroll_background()
        finally:
            self._layout_syncing = False

    def _on_canvas_configure(self, event: tk.Event) -> None:
        """React to canvas size changes without re-entering a height feedback loop.

        Normal mode only refreshes the scrollregion. Content height is owned by
        ``_on_grid_configure`` / mount / refresh paths. Compact-summary mode may
        still size the summary panel to the canvas height.
        """
        if self._layout_syncing:
            return
        if self._compact_summary_window is not None:
            self._layout_syncing = True
            try:
                self._apply_grid_vertical_size(max(int(event.height), 1))
            finally:
                self._layout_syncing = False
        self._sync_grid_scroll_background()

    def _sync_row_heights(self) -> None:
        """Keep the fixed label column row heights aligned with interval columns."""
        if self._grid_shows_compact_summary():
            return
        if not self._intervals:
            return
        ref = self._intervals[0]
        header = ref.get("header_label")
        widgets = ref.get("widgets") or []
        if header is None or len(widgets) < 3:
            return
        self.update_idletasks()
        label_children = [
            w
            for w in self._row_header_frame.winfo_children()
            if isinstance(w, (tk.Label, ttk.Label))
        ][:4]
        row_widgets = [header, widgets[0], widgets[1], widgets[2]]
        interval_heights: list[int] = []
        for widget in row_widgets:
            # Use requested (natural) height so maximize cannot lock in stretched minsizes.
            height = int(widget.winfo_reqheight())
            if height <= 1:
                return
            interval_heights.append(height)
        # Prefer the leftmost label column as the row-height authority.
        if len(label_children) == 4:
            label_heights = [int(w.winfo_reqheight()) for w in label_children]
            if any(h <= 1 for h in label_heights):
                heights = interval_heights
            else:
                heights = [
                    max(lh, ih) for lh, ih in zip(label_heights, interval_heights)
                ]
        else:
            heights = interval_heights
        # weight=0: rows stay at natural height; do not stretch to fill the canvas.
        for row, minsize in enumerate(heights):
            info = self._row_header_frame.grid_rowconfigure(row)
            if (
                int(info.get("minsize") or 0) == minsize
                and int(info.get("weight") or 0) == 0
            ):
                # Still ensure chunks match in case a new chunk was added.
                for chunk in self._grid_chunks:
                    cinfo = chunk["frame"].grid_rowconfigure(row)
                    if (
                        int(cinfo.get("minsize") or 0) != minsize
                        or int(cinfo.get("weight") or 0) != 0
                    ):
                        chunk["frame"].grid_rowconfigure(
                            row, weight=0, minsize=minsize, uniform=""
                        )
                continue
            opts = {"weight": 0, "minsize": minsize, "uniform": ""}
            self._row_header_frame.grid_rowconfigure(row, **opts)
            for chunk in self._grid_chunks:
                chunk["frame"].grid_rowconfigure(row, **opts)

    def _intensity_row_label_text(self) -> str:
        """Row-header text for intensity, reflecting log10 vs linear input mode."""
        units = "ph/\u00b5m\u00b2/s"
        if self._intensity_log_mode.get():
            return f"Log10 Intensity ({units})"
        return f"Intensity ({units})"

    def _build_row_labels(self) -> None:
        """Draw the fixed row-header column (col 0)."""
        header_font = tkfont.nametofont("TkDefaultFont", root=self).copy()
        header_font.configure(weight="bold")
        label_width = 28
        label_kwargs = dict(
            width=label_width,
            font=header_font,
            relief="ridge",
            padx=4,
            pady=4,
            bg=self._default_cell_bg,
            fg=self._default_cell_fg,
        )

        tk.Label(self._row_header_frame, text="", **label_kwargs).grid(
            row=0, column=0, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1)
        )
        tk.Label(
            self._row_header_frame, text="Spectrum", anchor=tk.W, **label_kwargs
        ).grid(row=1, column=0, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1))

        self._intensity_row_label = tk.Label(
            self._row_header_frame,
            text=self._intensity_row_label_text(),
            anchor=tk.W,
            **label_kwargs,
        )
        self._intensity_row_label.grid(
            row=2, column=0, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1)
        )

        tk.Label(
            self._row_header_frame, text="Duration (s)", anchor=tk.W, **label_kwargs
        ).grid(row=3, column=0, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1))
        self.update_idletasks()
        self._row_header_frame.configure(
            width=max(int(self._intensity_row_label.winfo_reqwidth()), 1)
        )

    # ── interval management ──────────────────────────────────────────

    def _wire_interval_widgets(self, interval_entry: dict) -> None:
        spectrum_var = interval_entry["spectrum_var"]
        intensity_var = interval_entry["intensity_var"]
        duration_var = interval_entry["duration_var"]
        header_label = interval_entry["header_label"]
        spec_btn, intensity_entry, duration_entry = interval_entry["widgets"]

        spec_btn.config(
            command=lambda e=interval_entry: self._open_spectrum_picker(
                e["spectrum_var"], self._intervals.index(e)
            )
        )

        def _on_click(_event, e=interval_entry) -> None:
            try:
                self._select_column(self._intervals.index(e))
            except ValueError:
                pass

        for w in (header_label, spec_btn, intensity_entry, duration_entry):
            w.bind("<Button-1>", _on_click, add="+")
        intensity_entry.bind(
            "<Return>",
            lambda _e, v=intensity_var: self._on_interval_numeric_change(
                v, interval=self._interval_for_var(v)
            ),
        )
        duration_entry.bind(
            "<Return>", lambda _e: self._on_interval_numeric_change(duration_var)
        )

    def _mount_interval_widgets(
        self,
        interval_entry: dict,
        interval_index: int,
        *,
        defer_ui_sync: bool = False,
    ) -> None:
        """Create and grid widgets for one interval inside the correct canvas chunk."""
        chunk_idx, grid_col = self._chunk_index_for_interval(interval_index)
        chunk = self._ensure_grid_chunk(chunk_idx)
        parent = chunk["frame"]
        col_idx = interval_index + 1

        header_font = tkfont.nametofont("TkDefaultFont", root=self).copy()
        header_font.configure(weight="bold")

        header_label = tk.Label(
            parent,
            text=f"{col_idx}",
            font=header_font,
            relief="ridge",
            padx=4,
            pady=4,
            anchor=tk.CENTER,
            width=18,
            bg=self._default_cell_bg,
        )
        header_label.grid(row=0, column=grid_col, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1))

        spec_btn = tk.Button(
            parent,
            textvariable=interval_entry["spectrum_var"],
            width=18,
            relief="raised",
            bg=self._default_cell_bg,
        )
        spec_btn.grid(row=1, column=grid_col, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1))

        intensity_entry = tk.Entry(
            parent,
            textvariable=interval_entry["intensity_var"],
            width=18,
            bg=self._default_cell_bg,
        )
        intensity_entry.grid(
            row=2, column=grid_col, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1)
        )
        intensity_entry.bind(
            "<FocusOut>",
            lambda _e, v=interval_entry["intensity_var"]: self._on_interval_numeric_change(
                v, interval=interval_entry
            ),
        )

        duration_entry = tk.Entry(
            parent,
            textvariable=interval_entry["duration_var"],
            width=18,
            bg=self._default_cell_bg,
        )
        duration_entry.grid(
            row=3, column=grid_col, sticky=tk.NSEW, padx=(0, 1), pady=(0, 1)
        )
        duration_entry.bind(
            "<FocusOut>",
            lambda _e, v=interval_entry["duration_var"]: self._on_interval_numeric_change(v),
        )

        interval_entry["header_label"] = header_label
        interval_entry["widgets"] = [spec_btn, intensity_entry, duration_entry]
        self._wire_interval_widgets(interval_entry)

        if self._measured_col_width is None:
            self._canvas.update_idletasks()
            measured = header_label.winfo_width()
            if measured > 0:
                self._measured_col_width = measured
                self._measured_col_pitch = measured + 1
                self._canvas.configure(xscrollincrement=self._measured_col_pitch)
        if not defer_ui_sync:
            self._sync_row_heights()

    def _refresh_interval_ui(self) -> None:
        """Update grid layout, light monitor, duration, and row heights."""
        self._on_grid_configure()
        self._render_light_monitor()
        self._sync_duration_to_parent()

    def _rebuild_all_interval_widgets(self) -> None:
        self._teardown_interval_widgets()
        for idx, entry in enumerate(self._intervals):
            if self._interval_is_compact(entry):
                continue
            self._mount_interval_widgets(entry, idx)
        self._update_compact_grid_display()
        self._on_grid_configure()

    def _add_interval(self, *, defer_ui_refresh: bool = False) -> None:
        interval_entry: dict = {
            "spectrum_var": tk.StringVar(master=self, value=self._NONE_LABEL),
            "intensity_var": tk.StringVar(master=self, value=""),
            "duration_var": tk.StringVar(master=self, value=""),
            "widgets": [],
            "header_label": None,
            "block_kind": "interval",
        }
        self._intervals.append(interval_entry)
        self._mount_interval_widgets(
            interval_entry,
            len(self._intervals) - 1,
            defer_ui_sync=defer_ui_refresh,
        )
        if not defer_ui_refresh:
            self._refresh_interval_ui()
        self._mark_stimulus_dirty()

    def _add_interval_data(self, spectrum: str, intensity: float, duration: float) -> dict:
        """Append one interval without mounting per-column grid widgets."""
        spectrum = str(spectrum).strip()
        interval_entry: dict = {
            "spectrum_var": tk.StringVar(master=self, value=spectrum),
            "intensity_var": tk.StringVar(
                master=self,
                value="0"
                if spectrum == "Dark"
                else self._format_intensity_for_display(float(intensity)),
            ),
            "duration_var": tk.StringVar(master=self, value=f"{float(duration):g}"),
            "widgets": [],
            "header_label": None,
            "compact": True,
            "block_kind": "interval",
            "dark_locked": spectrum == "Dark",
        }
        self._intervals.append(interval_entry)
        return interval_entry

    def _delete_interval(self) -> None:
        if len(self._intervals) <= 1:
            messagebox.showinfo(
                "Delete Interval",
                "At least one interval must remain.",
                parent=self,
            )
            return

        # Ask which interval to delete
        choice: dict[str, int | None] = {"idx": None}
        picker = tk.Toplevel(self)
        picker.title("Delete Interval")
        picker.transient(self)
        picker.resizable(False, False)

        outer = ttk.Frame(picker, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(outer, text="Select an interval to delete:").pack(anchor=tk.W, pady=(0, 6))

        n = len(self._intervals)
        list_frame = ttk.Frame(outer)
        list_frame.pack(fill=tk.BOTH, expand=True)
        listbox = tk.Listbox(list_frame, height=min(10, n),
                             exportselection=False, width=24)
        listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        if n > 10:
            sb = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=listbox.yview)
            sb.pack(side=tk.RIGHT, fill=tk.Y)
            listbox.configure(yscrollcommand=sb.set)
        for i in range(n):
            listbox.insert(tk.END, f"Interval {i + 1}")
        listbox.selection_set(0)

        def _confirm(_event=None) -> None:
            sel = listbox.curselection()
            if sel:
                choice["idx"] = int(sel[0])
            picker.destroy()

        listbox.bind("<Double-Button-1>", _confirm)

        btn_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        btn_row.pack(fill=tk.X)
        ttk.Button(btn_row, text="Delete", command=_confirm).pack(side=tk.LEFT)
        ttk.Button(btn_row, text="Cancel", command=picker.destroy).pack(side=tk.RIGHT)

        picker.protocol("WM_DELETE_WINDOW", picker.destroy)
        picker.grab_set()
        picker.wait_window()

        idx = choice["idx"]
        if idx is None:
            return

        del self._intervals[idx]
        self._selected_col = None
        self._rebuild_all_interval_widgets()
        self._render_light_monitor()
        self._render_preview(None)
        self._sync_duration_to_parent()
        self._mark_stimulus_dirty()

    # ── column selection & preview ───────────────────────────────────

    def _focus_is_text_entry(self) -> bool:
        """True when focus is in a text field where Left/Right move the caret."""
        w = self.focus_get()
        return isinstance(w, (tk.Entry, ttk.Entry, tk.Text, tk.Spinbox))

    def _ensure_column_visible(self, col_idx: int) -> None:
        """Scroll the interval grid so ``col_idx`` is within the visible x-range."""
        if self._grid_shows_compact_summary():
            return
        if col_idx < 0 or col_idx >= len(self._intervals):
            return
        pitch = self._column_pitch_px()
        if pitch <= 0:
            return
        self._canvas.update_idletasks()
        try:
            region = [float(v) for v in str(self._canvas.cget("scrollregion")).split()]
        except (TypeError, ValueError):
            return
        if len(region) != 4:
            return
        total_w = region[2] - region[0]
        if total_w <= 0:
            return
        canvas_w = max(int(self._canvas.winfo_width()), 1)
        x0 = float(col_idx * pitch)
        x1 = x0 + float(pitch)
        left_frac, right_frac = self._canvas.xview()
        left_px = left_frac * total_w
        right_px = right_frac * total_w
        if x0 < left_px:
            self._canvas.xview_moveto(max(0.0, x0 / total_w))
        elif x1 > right_px:
            self._canvas.xview_moveto(max(0.0, (x1 - canvas_w) / total_w))

    def _navigate_selected_interval(self, delta: int) -> str | None:
        """Move selection by ``delta`` columns. Returns ``break`` when handled."""
        if self._focus_is_text_entry():
            return None
        n = len(self._intervals)
        if n == 0:
            return "break"
        if self._selected_col is None:
            idx = 0 if delta > 0 else n - 1
        else:
            idx = max(0, min(n - 1, self._selected_col + delta))
        if idx != self._selected_col:
            self._select_column(idx)
        return "break"

    def _on_interval_nav_left(self, _event=None):
        return self._navigate_selected_interval(-1)

    def _on_interval_nav_right(self, _event=None):
        return self._navigate_selected_interval(1)

    def _select_column(self, col_idx: int) -> None:
        """Highlight col_idx, sync the light-monitor band, and render its spectrum."""
        if col_idx < 0 or col_idx >= len(self._intervals):
            return
        entry = self._intervals[col_idx]
        if self._interval_is_compact(entry):
            self._selected_col = col_idx
            self._render_interval_preview(entry)
            self._update_monitor_selection_band(col_idx)
            self._ensure_column_visible(col_idx)
            return
        # deselect previous
        if self._selected_col is not None and self._selected_col < len(self._intervals):
            prev = self._intervals[self._selected_col]
            if not self._interval_is_compact(prev):
                prev["header_label"].config(bg=self._default_cell_bg)
                for w in prev["widgets"]:
                    w.config(bg=self._default_cell_bg)

        self._selected_col = col_idx

        # highlight new selection
        entry["header_label"].config(bg=self._SELECT_BG)
        for w in entry["widgets"]:
            w.config(bg=self._SELECT_BG)

        self._render_interval_preview(entry)
        self._update_monitor_selection_band(col_idx)
        self._ensure_column_visible(col_idx)

    def _spectrum_for_name(
        self, name: str
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Return (wlen, intensity) for a named spectrum, or None if env not loaded."""
        return _library_spectrum_arrays(self._parent, name)

    def _interval_target_intensity(self, interval: dict) -> float:
        """Linear integrated intensity for preview scaling; invalid/blank → 0."""
        intensity = self._parse_interval_intensity_for_monitor(interval)
        if intensity is None:
            return 0.0
        return float(intensity)

    def _render_interval_preview(self, interval: dict) -> None:
        """Render the spectrum preview scaled to this interval's intensity."""
        self._render_preview(
            interval["spectrum_var"].get(),
            self._interval_target_intensity(interval),
        )

    def _render_preview(
        self, name: str | None, target_intensity: float | None = None
    ) -> None:
        ax = self._preview_ax
        ax.clear()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Intensity (ph/\u00b5m\u00b2/nm/s)")
        ax.set_xlim(100, 1000)
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0), useMathText=True)

        if name is None or name == self._NONE_LABEL:
            ax.set_title("Select a spectrum")
        else:
            result = self._spectrum_for_name(name)
            if result is None:
                ax.set_title(f"{name}")
                ax.text(0.5, 0.5, "Spectra not loaded yet.",
                        ha="center", va="center", transform=ax.transAxes)
            else:
                wlen, intensity = result
                plot_intensity = np.asarray(intensity, dtype=float)
                if target_intensity is not None:
                    try:
                        plot_intensity = _scale_spectrum_to_integrated_intensity(
                            wlen,
                            intensity,
                            float(target_intensity),
                            ref_name=str(name),
                        )
                    except ValueError:
                        ax.set_title(f"{name}")
                        ax.text(
                            0.5,
                            0.5,
                            "Cannot scale spectrum to intensity.",
                            ha="center",
                            va="center",
                            transform=ax.transAxes,
                        )
                        self._preview_canvas.draw_idle()
                        return
                ax.plot(wlen, plot_intensity, color="k")
                ax.set_title(f"{name}")
                # Flat/zero spectra (e.g. Dark) otherwise autoscaled symmetrically
                # around zero; keep zero near the bottom of the pane.
                imax = float(np.max(plot_intensity)) if plot_intensity.size else 0.0
                if not np.isfinite(imax) or imax <= 0.0:
                    ax.set_ylim(-1e-3, 1e-2)
                else:
                    ax.set_ylim(-0.02 * imax, imax * 1.05)

        self._preview_canvas.draw_idle()

    def _parse_interval_intensity_for_monitor(self, interval: dict) -> float | None:
        """Linear intensity for monitor plot/hover; blank entry is treated as zero."""
        raw_i = interval["intensity_var"].get().strip()
        if raw_i == "":
            return 0.0
        i = self._parse_intensity_text(raw_i, interval)
        if i is None:
            return None
        return max(0.0, float(i))

    def _iter_monitor_interval_segments(
        self,
    ) -> list[tuple[float, float, int, float]]:
        """Return ``(t_start, t_end, col_idx, intensity)`` segments for the monitor."""
        segments: list[tuple[float, float, int, float]] = []
        t = 0.0
        for col_idx, interval in enumerate(self._intervals):
            raw_d = interval["duration_var"].get().strip()
            try:
                d = float(raw_d)
            except ValueError:
                continue
            if d <= 0:
                continue
            i = self._parse_interval_intensity_for_monitor(interval)
            if i is None:
                continue
            t_end = t + d
            segments.append((t, t_end, col_idx, i))
            t = t_end
        return segments


    def _build_light_monitor_interval_bounds(
        self,
    ) -> list[tuple[float, float, int]]:
        """Return ``(t_start, t_end, col_idx)`` for each interval on the monitor timeline."""
        return [
            (t_start, t_end, col_idx)
            for t_start, t_end, col_idx, _ in self._iter_monitor_interval_segments()
        ]

    def _interval_index_at_monitor_time(self, t: float) -> int | None:
        """Map a monitor x-coordinate (seconds) to a grid column index, if any."""
        bounds = self._monitor_interval_bounds
        if not bounds:
            return None
        # Inclusive end times; at shared boundaries prefer the later interval.
        found: int | None = None
        for t_start, t_end, col_idx in bounds:
            if t_start <= t <= t_end:
                found = col_idx
        return found

    def _update_monitor_selection_band(self, col_idx: int | None) -> None:
        """Draw or clear the sticky vertical highlight for the selected interval."""
        if self._monitor_selection_patch is not None:
            try:
                self._monitor_selection_patch.remove()
            except ValueError:
                pass
            self._monitor_selection_patch = None
        if col_idx is not None:
            spans = [
                (t_start, t_end)
                for t_start, t_end, c in self._monitor_interval_bounds
                if c == col_idx
            ]
            if spans:
                t_start = min(s for s, _ in spans)
                t_end = max(e for _, e in spans)
                self._monitor_selection_patch = self._monitor_ax.axvspan(
                    t_start,
                    t_end,
                    color=self._SELECT_BG,
                    alpha=0.65,
                    zorder=0,
                )
        self._monitor_canvas.draw_idle()

    def _monitor_time_from_event(self, event) -> float | None:
        """Return time (s) under the cursor anywhere in the monitor axes, or None."""
        if event.x is None or event.y is None:
            return None
        if event.inaxes == self._monitor_ax and event.xdata is not None:
            return float(event.xdata)
        try:
            xdata, _yd = self._monitor_ax.transData.inverted().transform(
                (event.x, event.y)
            )
        except Exception:
            return None
        xmin, xmax = self._monitor_ax.get_xlim()
        if xmin <= xmax:
            if xmin <= xdata <= xmax:
                return float(xdata)
        elif xmax <= xdata <= xmin:
            return float(xdata)
        return None

    def _clamp_monitor_xlim(self, xmin: float, xmax: float) -> tuple[float, float]:
        """Keep the visible x-range inside the full stimulus timeline."""
        if self._monitor_xlim_full is None:
            return xmin, xmax
        lo, hi = self._monitor_xlim_full
        span = xmax - xmin
        full_span = hi - lo
        if span >= full_span:
            return lo, hi
        if xmin < lo:
            xmax = lo + span
            xmin = lo
        if xmax > hi:
            xmin = hi - span
            xmax = hi
        return xmin, xmax

    def _monitor_zoom_center(self, t: float | None) -> float:
        """Pick a time (s) to anchor zoom/pan, preferring the cursor position."""
        ax = self._monitor_ax
        xmin, xmax = ax.get_xlim()
        if t is not None and xmin <= t <= xmax:
            return float(t)
        if self._monitor_last_hover_t is not None and xmin <= self._monitor_last_hover_t <= xmax:
            return float(self._monitor_last_hover_t)
        return float(0.5 * (xmin + xmax))

    def _apply_monitor_x_zoom(self, step: int, t_center: float | None) -> None:
        """Zoom the monitor time axis in (step > 0) or out (step < 0) around ``t_center``."""
        if self._monitor_xlim_full is None or step == 0:
            return
        ax = self._monitor_ax
        lo, hi = self._monitor_xlim_full
        full_span = hi - lo
        if full_span <= 0:
            return
        xmin, xmax = ax.get_xlim()
        span = xmax - xmin
        if span <= 0:
            return
        t_anchor = self._monitor_zoom_center(t_center)
        zoom_factor = 0.82 if step > 0 else 1.0 / 0.82
        new_span = span * zoom_factor
        min_span = max(full_span / 5000.0, 1e-6)
        new_span = max(min_span, min(new_span, full_span))
        frac = (t_anchor - xmin) / span
        new_xmin = t_anchor - frac * new_span
        new_xmax = new_xmin + new_span
        new_xmin, new_xmax = self._clamp_monitor_xlim(new_xmin, new_xmax)
        ax.set_xlim(new_xmin, new_xmax)
        if self._selected_col is not None:
            self._update_monitor_selection_band(self._selected_col)
        else:
            self._monitor_canvas.draw_idle()

    def _apply_monitor_x_pan(self, step: int, t_center: float | None) -> None:
        """Pan the monitor time axis (Shift + scroll)."""
        if self._monitor_xlim_full is None or step == 0:
            return
        ax = self._monitor_ax
        xmin, xmax = ax.get_xlim()
        span = xmax - xmin
        if span <= 0:
            return
        # Scroll up moves the view toward earlier times.
        dx = -0.12 * span * step
        new_xmin, new_xmax = self._clamp_monitor_xlim(xmin + dx, xmax + dx)
        ax.set_xlim(new_xmin, new_xmax)
        if self._selected_col is not None:
            self._update_monitor_selection_band(self._selected_col)
        else:
            self._monitor_canvas.draw_idle()

    def _reset_monitor_xlim(self) -> None:
        """Restore the full timeline view on the light monitor."""
        if self._monitor_xlim_full is None:
            return
        lo, hi = self._monitor_xlim_full
        self._monitor_ax.set_xlim(lo, hi)
        if self._selected_col is not None:
            self._update_monitor_selection_band(self._selected_col)
        else:
            self._monitor_canvas.draw_idle()

    def _monitor_time_from_tk_wheel(self, event: tk.Event) -> float | None:
        """Map a Tk wheel event position to monitor time (s)."""
        try:
            inv = self._monitor_ax.transData.inverted()
            widget = self._monitor_canvas.get_tk_widget()
            height = max(widget.winfo_height(), 1)
            xdata, _ydata = inv.transform((float(event.x), float(height - event.y)))
        except Exception:
            return None
        xmin, xmax = self._monitor_ax.get_xlim()
        if xmin <= xmax:
            if xmin <= xdata <= xmax:
                return float(xdata)
        elif xmax <= xdata <= xmin:
            return float(xdata)
        return None

    def _wheel_step_from_tk_event(self, event: tk.Event) -> int:
        if getattr(event, "num", None) == 4:
            return 1
        if getattr(event, "num", None) == 5:
            return -1
        if event.delta:
            return 1 if event.delta > 0 else -1
        return 0

    def _on_monitor_scroll(self, event) -> None:
        """Matplotlib scroll: zoom time axis; Shift+scroll pans."""
        if event.inaxes != self._monitor_ax:
            return
        raw_step = getattr(event, "step", None)
        if raw_step is None or raw_step == 0:
            btn = getattr(event, "button", None)
            if btn == "up":
                step = 1
            elif btn == "down":
                step = -1
            else:
                return
        else:
            step = int(raw_step)
        t = self._monitor_time_from_event(event)
        if bool(getattr(event, "key", None)) and "shift" in str(event.key).lower():
            self._apply_monitor_x_pan(step, t)
        else:
            self._apply_monitor_x_zoom(step, t)

    def _on_monitor_tk_wheel(self, event: tk.Event) -> str:
        """Tk wheel binding (Windows/Linux) for monitor zoom and pan."""
        step = self._wheel_step_from_tk_event(event)
        if step == 0:
            return "break"
        t = self._monitor_time_from_tk_wheel(event)
        shift = bool(getattr(event, "state", 0) & 0x1)
        if shift:
            self._apply_monitor_x_pan(step, t)
        else:
            self._apply_monitor_x_zoom(step, t)
        return "break"

    def _on_monitor_reset_zoom(self, _event: tk.Event | None = None) -> None:
        self._reset_monitor_xlim()
        return "break"

    def _on_monitor_click(self, event) -> None:
        """Click a monitor segment to select the matching interval column."""
        if getattr(event, "dblclick", False):
            return
        if getattr(event, "button", None) != 1:
            return
        t = self._monitor_time_from_event(event)
        if t is None:
            return
        col = self._interval_index_at_monitor_time(t)
        if col is None:
            return
        self._monitor_hover_col = col
        self._select_column(col)

    def _on_monitor_motion(self, event) -> None:
        """Hover over the light monitor to select an interval column and preview."""
        t = self._monitor_time_from_event(event)
        if t is None:
            if self._monitor_hover_col is not None:
                self._monitor_hover_col = None
                self._update_monitor_selection_band(self._selected_col)
            return
        self._monitor_last_hover_t = t
        col = self._interval_index_at_monitor_time(t)
        if col == self._monitor_hover_col:
            return
        self._monitor_hover_col = col
        if col is not None:
            self._select_column(col)
        else:
            self._update_monitor_selection_band(self._selected_col)

    def _on_monitor_leave(self, _event=None) -> None:
        """Keep the selection band when the pointer leaves the plot."""
        self._monitor_hover_col = None
        self._update_monitor_selection_band(self._selected_col)

    def _build_light_monitor_series(self) -> tuple[np.ndarray, np.ndarray] | None:
        """Build step-series arrays from interval intensity/duration values."""
        segments = self._iter_monitor_interval_segments()
        if not segments:
            return None
        times: list[float] = []
        intensities: list[float] = []
        for t_start, t_end, _col_idx, intensity in segments:
            if not times:
                times.append(t_start)
                intensities.append(float(intensity))
            elif t_start == times[-1]:
                # steps-post: y[i] applies on [x[i], x[i+1]); update boundary value
                intensities[-1] = float(intensity)
            elif t_start > times[-1]:
                times.append(t_start)
                intensities.append(float(intensity))
            times.append(t_end)
            intensities.append(float(intensity))
        if len(times) < 2:
            return None
        return np.asarray(times, dtype=float), np.asarray(intensities, dtype=float)

    def _apply_monitor_axis_limits(
        self, ax, t_hi: float, log_light: np.ndarray
    ) -> None:
        """Apply uniform padding so the trace is not flush against the axis spines."""
        t_hi = max(float(t_hi), 1e-9)
        xpad = max(0.05 * t_hi, 1e-6)
        ax.set_xlim(-xpad, t_hi + xpad)

        finite_light = log_light[np.isfinite(log_light)]
        if finite_light.size == 0:
            ylo, yhi = 0.0, 1.0
        else:
            ylo = float(np.min(finite_light))
            yhi = float(np.max(finite_light))
        yspan = yhi - ylo
        if yspan < 1e-10:
            ypad = max(1.0, abs(ylo) * 0.1)
        else:
            ypad = 0.12 * yspan
        ax.set_ylim(ylo - ypad, yhi + ypad)
        self._monitor_fig.subplots_adjust(left=0.16, right=0.99, bottom=0.22, top=0.90)
        self._monitor_xlim_full = tuple(ax.get_xlim())

    def _render_light_monitor(self) -> None:
        ax = self._monitor_ax
        ax.clear()
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"$\mathrm{Log}_{10}$ photons/$\mathrm{\mu m}^{2}$/s")
        self._monitor_selection_patch = None

        self._monitor_interval_bounds = self._build_light_monitor_interval_bounds()
        valid_cols = {col for _, _, col in self._monitor_interval_bounds}

        series = self._build_light_monitor_series()
        if series is None:
            ax.set_title("Light monitor")
            ax.text(
                0.5,
                0.5,
                "Add interval Intensity and Duration to display light monitor.",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            self._monitor_xlim_full = None
            self._update_monitor_selection_band(None)
            return

        x, intensity = series
        log_light = np.log10(np.maximum(intensity, 0.0) + 1.0)
        ax.plot(x, log_light, color="k", drawstyle="steps-post")
        t_hi = float(np.max(x))
        self._apply_monitor_axis_limits(ax, t_hi, log_light)
        ax.set_title(
            "Light monitor (scroll to zoom, Shift+scroll to pan; "
            "←/→ move between intervals)"
        )
        if self._selected_col is not None and self._selected_col in valid_cols:
            self._update_monitor_selection_band(self._selected_col)
        else:
            self._monitor_canvas.draw_idle()

    def _compute_total_duration(self) -> float:
        total = 0.0
        for interval in self._intervals:
            raw_d = interval["duration_var"].get().strip()
            if raw_d == "":
                continue
            try:
                d = float(raw_d)
            except ValueError:
                continue
            if d > 0:
                total += d
        return total

    def _sync_duration_to_parent(self) -> None:
        total = self._compute_total_duration()
        if total > 0:
            self._duration_var.set(f"{total:g}")
            self._parent._stimulus_duration_override = float(total)
        else:
            self._duration_var.set("")
            self._parent._stimulus_duration_override = None

    def _add_dark_interval(self, duration: float) -> None:
        self._add_interval()
        interval = self._intervals[-1]
        interval["spectrum_var"].set("Dark")
        interval["intensity_var"].set("0")
        interval["duration_var"].set(f"{duration:g}")
        self._apply_dark_intensity_lock(interval, "Dark")
        self._render_light_monitor()
        self._sync_duration_to_parent()

    def _on_declared_duration_commit(self, _event=None) -> None:
        raw = self._duration_var.get().strip()
        if raw == "":
            self._sync_duration_to_parent()
            return
        declared = self._parse_scientific_text(raw)
        if declared is None or not np.isfinite(declared) or declared <= 0:
            messagebox.showerror(
                "Duration",
                "Duration must be a positive number.",
                parent=self,
            )
            self._sync_duration_to_parent()
            return
        current = self._compute_total_duration()
        tol = 1e-9 * max(1.0, abs(current), abs(declared))
        if declared + tol < current:
            messagebox.showerror(
                "Duration",
                (
                    "Declared duration is shorter than the total interval duration.\n"
                    f"Declared: {declared:g} s\n"
                    f"Current total: {current:g} s"
                ),
                parent=self,
            )
            self._sync_duration_to_parent()
            return
        if declared > current + tol:
            self._add_dark_interval(declared - current)
            self._parent._status.config(
                text=(
                    f"Added Dark interval ({declared - current:g} s) "
                    "to match declared duration."
                )
            )
            return
        self._sync_duration_to_parent()

    # ── spectrum picker ──────────────────────────────────────────────

    def _open_spectrum_picker(self, spectrum_var: tk.StringVar, col_idx: int) -> None:
        def _on_select(name: str) -> None:
            spectrum_var.set(name)
            try:
                idx = col_idx if 0 <= col_idx < len(self._intervals) else None
            except Exception:
                idx = None
            if idx is not None:
                self._apply_dark_intensity_lock(self._intervals[idx], name)
            self._select_column(col_idx)
            self._mark_stimulus_dirty()

        _SpectrumPickerPopup(
            parent=self,
            app=self._parent,
            current=spectrum_var.get(),
            on_select=_on_select,
        )

    def _apply_dark_intensity_lock(self, interval: dict, spectrum_name: str) -> None:
        """For a Dark spectrum, force intensity to 0 and disable the entry.

        For any other spectrum, re-enable the intensity entry. If the entry
        was previously locked at 0 by this method, clear it so the user can
        enter the real intensity for the new spectrum.
        """
        intensity_var = interval["intensity_var"]
        was_dark = bool(interval.get("dark_locked", False))

        if spectrum_name == "Dark":
            intensity_var.set("0")
            interval["dark_locked"] = True
        else:
            interval["dark_locked"] = False
            if was_dark and intensity_var.get().strip() == "0":
                intensity_var.set("")

        if self._interval_is_compact(interval):
            self._render_light_monitor()
            return

        intensity_entry = interval["widgets"][1]
        if spectrum_name == "Dark":
            try:
                intensity_entry.configure(
                    state="disabled",
                    disabledbackground=self._default_cell_bg,
                    disabledforeground="black",
                )
            except tk.TclError:
                intensity_entry.configure(state="disabled")
        else:
            intensity_entry.configure(state="normal")

        self._render_light_monitor()

    # ── validation & close ───────────────────────────────────────────

    @staticmethod
    def _validate_float(var: tk.StringVar) -> None:
        raw = var.get().strip()
        if raw == "":
            return
        try:
            float(raw)
        except ValueError:
            var.set("")

    @staticmethod
    def _parse_scientific_text(raw: str) -> float | None:
        """Parse regular/scientific number text and a few common x10^ forms."""
        s = raw.strip()
        if not s:
            return None
        try:
            return float(s)
        except ValueError:
            pass

        # Also support strings like "2.5x10^3", "2.5 × 10^3", and "10^3".
        m = re.fullmatch(
            r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))?\s*(?:x|\*|×)\s*10\s*\^\s*([+-]?\d+)",
            s,
            flags=re.IGNORECASE,
        )
        if m is None:
            return None
        coeff_txt = m.group(1)
        exp_txt = m.group(2)
        coeff = 1.0 if coeff_txt is None or coeff_txt == "" else float(coeff_txt)
        exponent = int(exp_txt)
        return coeff * (10.0 ** exponent)

    def _interval_for_var(self, var: tk.StringVar) -> dict | None:
        for interval in self._intervals:
            if interval.get("intensity_var") is var:
                return interval
        return None

    def _validate_integer_like_intensity(self, var: tk.StringVar) -> None:
        raw = var.get().strip()
        if raw == "":
            return

        val = self._parse_scientific_text(raw)
        if val is None or not np.isfinite(val):
            var.set("")
            return

        nearest = round(val)
        tol = 1e-9 * max(1.0, abs(val))
        if abs(val - nearest) > tol:
            var.set("")
            return
        var.set(str(int(nearest)))

    @staticmethod
    def _looks_like_linear_scientific(raw: str) -> bool:
        """True when text is scientific notation meant as a linear value.

        Plain log10 exponents are numbers like ``12`` or ``12.5``. Forms such as
        ``10e5``, ``1e12``, or ``2.5x10^3`` are almost always linear intensity.
        """
        s = raw.strip()
        if not s:
            return False
        if re.fullmatch(
            r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)[eE][+-]?\d+",
            s,
        ):
            return True
        return (
            re.search(r"(?:x|\*|×)\s*10\s*\^", s, flags=re.IGNORECASE) is not None
        )

    def _intensity_from_raw_with_mode(self, raw: str, log_mode: bool) -> float | None:
        """Parse raw intensity text under a specific input mode into linear intensity."""
        if raw.strip() == "":
            return None
        # Scientific notation (1e12, 10e5, 2.5x10^3) is linear even in log mode.
        use_log = bool(log_mode) and not self._looks_like_linear_scientific(raw)
        if use_log:
            exponent = self._parse_scientific_text(raw)
            if exponent is None or not np.isfinite(exponent):
                return None
            try:
                val = 10.0 ** float(exponent)
            except (OverflowError, ValueError):
                return None
            if not np.isfinite(val):
                return None
            return float(val)
        val = self._parse_scientific_text(raw)
        if val is None or not np.isfinite(val):
            return None
        return float(val)

    def _parse_intensity_text(self, raw: str, interval: dict | None) -> float | None:
        """Parse entry intensity to linear units, forcing Dark intervals to zero."""
        if interval is not None:
            spectrum = str(interval["spectrum_var"].get()).strip()
            if spectrum == "Dark":
                return 0.0
        return self._intensity_from_raw_with_mode(raw, self._intensity_log_mode.get())

    def _format_intensity_for_display(self, linear_val: float) -> str:
        """Format linear intensity according to current input mode."""
        linear_val = float(linear_val)
        if self._intensity_log_mode.get():
            if linear_val <= 0:
                return ""
            return f"{np.log10(linear_val):g}"
        return f"{linear_val:g}"

    def _on_intensity_mode_toggle(self) -> None:
        """Convert non-dark entries between linear and log10 display values."""
        new_mode = bool(self._intensity_log_mode.get())
        old_mode = not new_mode
        intensity_label = getattr(self, "_intensity_row_label", None)
        if intensity_label is not None:
            intensity_label.configure(text=self._intensity_row_label_text())
            self.update_idletasks()
            self._row_header_frame.configure(
                width=max(int(intensity_label.winfo_reqwidth()), 1)
            )
        was_tracking = self._tracking_edits
        self._tracking_edits = False
        try:
            for interval in self._intervals:
                if str(interval["spectrum_var"].get()).strip() == "Dark":
                    interval["intensity_var"].set("0")
                    continue
                raw = interval["intensity_var"].get().strip()
                linear_val = self._intensity_from_raw_with_mode(raw, old_mode)
                if linear_val is None:
                    continue
                interval["intensity_var"].set(self._format_intensity_for_display(linear_val))
        finally:
            self._tracking_edits = was_tracking
        self._render_light_monitor()

    def _on_interval_numeric_change(
        self,
        var: tk.StringVar,
        interval: dict | None = None,
    ) -> None:
        if interval is not None:
            raw = var.get().strip()
            if raw != "":
                linear_val = self._parse_intensity_text(raw, interval)
                if linear_val is None or linear_val < 0:
                    var.set("")
                elif str(interval["spectrum_var"].get()).strip() == "Dark":
                    var.set("0")
                else:
                    var.set(self._format_intensity_for_display(linear_val))
        else:
            self._validate_float(var)
        self._render_light_monitor()
        if self._selected_col is not None and self._selected_col < len(self._intervals):
            self._render_interval_preview(self._intervals[self._selected_col])
        self._sync_duration_to_parent()
        self._mark_stimulus_dirty()

    # ── stimulus save / load ─────────────────────────────────────────

    def _collect_blocks_for_save(self) -> list[dict] | None:
        """Read the builder state into a v2 block list, or None if invalid.

        Manual columns serialize to ``interval`` blocks, one per entry.
        """
        out: list[dict] = []
        for i, interval in enumerate(self._intervals):
            spectrum = interval["spectrum_var"].get().strip()
            if not spectrum or spectrum == self._NONE_LABEL:
                messagebox.showerror(
                    "Save stimulus",
                    f"Interval {i + 1}: choose a spectrum before saving.",
                    parent=self,
                )
                return None
            raw_i = interval["intensity_var"].get().strip()
            raw_d = interval["duration_var"].get().strip()
            intensity = self._parse_intensity_text(raw_i, interval)
            if intensity is None:
                intensity = 0.0 if raw_i == "" else None
            if intensity is None:
                messagebox.showerror(
                    "Save stimulus",
                    f"Interval {i + 1}: intensity is not a number.",
                    parent=self,
                )
                return None
            try:
                duration = float(raw_d) if raw_d else 0.0
            except ValueError:
                messagebox.showerror(
                    "Save stimulus",
                    f"Interval {i + 1}: duration is not a number.",
                    parent=self,
                )
                return None
            if duration <= 0:
                messagebox.showerror(
                    "Save stimulus",
                    f"Interval {i + 1}: duration must be positive.",
                    parent=self,
                )
                return None
            if intensity < 0:
                messagebox.showerror(
                    "Save stimulus",
                    f"Interval {i + 1}: intensity must be non-negative.",
                    parent=self,
                )
                return None
            out.append(
                library_storage.make_interval_block(spectrum, intensity, duration)
            )
        if not out:
            messagebox.showerror(
                "Save stimulus",
                "Add at least one valid interval before saving.",
                parent=self,
            )
            return None
        return out

    def _name_is_taken(self, name: str) -> bool:
        if name in LABELS:
            return True
        if name in self._parent._custom_stimuli:
            return True
        return False

    def _prompt_for_unique_name(self, default: str = "") -> str | None:
        while True:
            name = simpledialog.askstring(
                "Save stimulus",
                "Enter a name for this stimulus:",
                initialvalue=default,
                parent=self,
            )
            if name is None:
                return None
            name = name.strip()
            if not name:
                messagebox.showerror(
                    "Save stimulus", "Name cannot be empty.", parent=self
                )
                continue
            if self._name_is_taken(name):
                messagebox.showerror(
                    "Save stimulus",
                    f"A stimulus named {name!r} already exists. Choose a different name.",
                    parent=self,
                )
                continue
            return name

    def _on_save_stimulus(self) -> None:
        blocks = self._collect_blocks_for_save()
        if blocks is None:
            return
        default_name = ""
        if self._loaded_source_name:
            if self._loaded_source_name in LABELS:
                default_name = f"{self._loaded_source_name}_edited"
            else:
                default_name = self._loaded_source_name
        name = self._prompt_for_unique_name(default=default_name)
        if name is None:
            return
        total_duration = self._compute_total_duration()
        try:
            for block in blocks:
                if block.get("type") == "interval":
                    self._parent._ensure_custom_spectrum_persisted(
                        str(block.get("spectrum_ref", "")).strip()
                    )
            spec = _serialize_stimulus_spec(
                name, blocks, total_duration, self._parent._custom_spectra
            )
            self._parent._register_custom_stimulus(name, spec)
        except (ValueError, OSError) as exc:
            messagebox.showerror("Save stimulus", str(exc), parent=self)
            return
        self._parent._set_selected_stimulus(name)
        self._loaded_source_name = name
        self._mark_stimulus_clean()
        messagebox.showinfo(
            "Save stimulus",
            f"Saved stimulus {name!r}. It now appears in the run dropdown.",
            parent=self,
        )

    def _has_configured_intervals(self) -> bool:
        """True when the builder already has intervals beyond the empty default."""
        if len(self._intervals) > 1:
            return True
        if not self._intervals:
            return False
        interval = self._intervals[0]
        spectrum = str(interval["spectrum_var"].get()).strip()
        intensity = str(interval["intensity_var"].get()).strip()
        duration = str(interval["duration_var"].get()).strip()
        return (
            spectrum != self._NONE_LABEL
            or bool(intensity)
            or bool(duration)
        )

    def _prompt_replace_or_append_load(self) -> str | None:
        """Ask how to apply a load when intervals already exist."""
        result: dict[str, str | None] = {"choice": None}
        picker = tk.Toplevel(self)
        picker.title("Load stimulus")
        picker.transient(self)
        picker.resizable(False, False)

        outer = ttk.Frame(picker, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(
            outer,
            text=(
                "The builder already has stimulus intervals.\n"
                "How should the load be applied?"
            ),
            justify=tk.LEFT,
        ).pack(anchor=tk.W, pady=(0, 8))

        def _choose(choice: str) -> None:
            result["choice"] = choice
            picker.destroy()

        btn_row = ttk.Frame(outer)
        btn_row.pack(fill=tk.X)
        ttk.Button(
            btn_row, text="Replace existing", command=lambda: _choose("replace")
        ).pack(side=tk.LEFT)
        ttk.Button(
            btn_row, text="Append to end", command=lambda: _choose("append")
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(btn_row, text="Cancel", command=picker.destroy).pack(side=tk.RIGHT)

        picker.protocol("WM_DELETE_WINDOW", picker.destroy)
        picker.grab_set()
        picker.wait_window()
        return result["choice"]

    def _prompt_large_load_widget_display(self, n_intervals: int) -> bool | None:
        """Ask whether to mount a column widget per interval for a large load.

        Returns ``True`` to show columns, ``False`` for compact view, or ``None``
        if the user cancelled.
        """
        result: dict[str, bool | None] = {"choice": None}
        picker = tk.Toplevel(self)
        picker.title("Large stimulus load")
        picker.transient(self)
        picker.resizable(False, False)

        outer = ttk.Frame(picker, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(
            outer,
            text=(
                f"This load has {n_intervals:,} intervals "
                f"(more than {self._DEFER_REFRESH_INTERVAL_THRESHOLD:,}).\n\n"
                "Show a spreadsheet column for every interval?\n"
                "Compact view keeps all interval data for the monitor, save, "
                "and model run, but hides the per-interval columns for speed."
            ),
            justify=tk.LEFT,
        ).pack(anchor=tk.W, pady=(0, 8))

        def _choose(show_widgets: bool) -> None:
            result["choice"] = show_widgets
            picker.destroy()

        btn_row = ttk.Frame(outer)
        btn_row.pack(fill=tk.X)
        ttk.Button(
            btn_row,
            text="Show interval columns",
            command=lambda: _choose(True),
        ).pack(side=tk.LEFT)
        ttk.Button(
            btn_row,
            text="Compact view (faster)",
            command=lambda: _choose(False),
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(btn_row, text="Cancel", command=picker.destroy).pack(side=tk.RIGHT)

        picker.protocol("WM_DELETE_WINDOW", picker.destroy)
        picker.grab_set()
        picker.wait_window()
        return result["choice"]

    def _apply_library_spec(
        self,
        spec: dict,
        source_name: str | None = None,
        *,
        error_title: str = "Load from library",
    ) -> str | None:
        """Load a stimulus spec into the grid; return 'replace', 'append', or None."""
        append = False
        if self._has_configured_intervals():
            choice = self._prompt_replace_or_append_load()
            if choice is None:
                return None
            append = choice == "append"
        try:
            if not self._populate_grid_from_spec(
                spec, source_name=source_name, append=append
            ):
                return None
        except ValueError as exc:
            messagebox.showerror(error_title, str(exc), parent=self)
            return None
        return "append" if append else "replace"

    def _on_load_from_library(self) -> None:
        self._parent._reload_stimulus_library_from_disk()
        names = sorted(self._parent._custom_stimuli.keys())
        if not names:
            messagebox.showinfo(
                "Load from library",
                "No stimuli in the library yet. Use Save to add one.",
                parent=self,
            )
            return

        choice = self._show_load_library_picker(names)
        if choice is None:
            return
        self._load_library_stimulus(choice)

    def _show_load_library_picker(self, names: list[str]) -> str | None:
        """Picker for saved stimuli. Returns the chosen name, or None if cancelled."""
        result: dict[str, str | None] = {"choice": None}
        picker = tk.Toplevel(self)
        picker.title("Load from library")
        picker.transient(self)

        outer = ttk.Frame(picker, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(outer, text="Select a stimulus from the library:").pack(
            anchor=tk.W, pady=(0, 6)
        )
        list_frame = ttk.Frame(outer)
        list_frame.pack(fill=tk.BOTH, expand=True)
        listbox = tk.Listbox(
            list_frame, height=min(10, max(1, len(names))),
            exportselection=False, width=30,
        )
        listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=listbox.yview
        )
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        listbox.configure(yscrollcommand=sb.set)
        for n in names:
            listbox.insert(tk.END, n)
        if names:
            listbox.selection_set(0)

        def _confirm(_event=None) -> None:
            sel = listbox.curselection()
            if sel:
                result["choice"] = str(listbox.get(sel[0]))
                picker.destroy()

        listbox.bind("<Double-Button-1>", _confirm)

        btn_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        btn_row.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Button(btn_row, text="Load", command=_confirm).pack(side=tk.LEFT)
        ttk.Button(btn_row, text="Cancel", command=picker.destroy).pack(side=tk.RIGHT)
        picker.protocol("WM_DELETE_WINDOW", picker.destroy)
        picker.grab_set()
        picker.wait_window()
        return result["choice"]

    def _load_library_stimulus(self, chosen: str) -> None:
        spec = self._parent._custom_stimuli.get(chosen)
        if spec is None:
            messagebox.showerror(
                "Load from library",
                f"Stimulus {chosen!r} could not be loaded.",
                parent=self,
            )
            return
        action = self._apply_library_spec(spec, source_name=chosen)
        if action is None:
            return
        verb = "Appended" if action == "append" else "Loaded"
        self._parent._status.config(text=f"{verb} stimulus {chosen!r} from library.")

    def _on_load_manuscript_stimulus(self) -> None:
        if self._parent._env is None:
            messagebox.showwarning(
                "Load manuscript stimulus",
                "Spectra are still loading. Please try again in a moment.",
                parent=self,
            )
            return
        choice: dict[str, str | None] = {"value": None}
        picker = tk.Toplevel(self)
        picker.title("Load manuscript stimulus")
        picker.transient(self)
        picker.resizable(False, False)

        outer = ttk.Frame(picker, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(outer, text="Select a manuscript stimulus:").pack(
            anchor=tk.W, pady=(0, 6)
        )

        list_frame = ttk.Frame(outer)
        list_frame.pack(fill=tk.BOTH, expand=True)
        listbox = tk.Listbox(
            list_frame, height=min(10, len(LABELS)), exportselection=False, width=32
        )
        listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        for n in LABELS:
            listbox.insert(tk.END, n)
        listbox.selection_set(0)

        def _confirm(_event=None) -> None:
            sel = listbox.curselection()
            if sel:
                choice["value"] = str(listbox.get(sel[0]))
            picker.destroy()

        listbox.bind("<Double-Button-1>", _confirm)
        btn_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        btn_row.pack(fill=tk.X)
        ttk.Button(btn_row, text="Load", command=_confirm).pack(side=tk.LEFT)
        ttk.Button(btn_row, text="Cancel", command=picker.destroy).pack(side=tk.RIGHT)
        picker.protocol("WM_DELETE_WINDOW", picker.destroy)
        picker.grab_set()
        picker.wait_window()

        chosen = choice["value"]
        if chosen is None:
            return
        try:
            spec = self._parent._build_stimulus_spec_from_label(chosen)
        except ValueError as exc:
            messagebox.showerror("Load manuscript stimulus", str(exc), parent=self)
            return
        action = self._apply_library_spec(
            spec,
            source_name=chosen,
            error_title="Load manuscript stimulus",
        )
        if action is None:
            return
        verb = "Appended" if action == "append" else "Loaded"
        self._parent._status.config(
            text=f"{verb} manuscript stimulus {chosen!r}."
        )

    def _populate_grid_from_spec(
        self,
        spec: dict,
        source_name: str | None = None,
        *,
        append: bool = False,
    ) -> bool:
        """Load a v2 block-based stimulus spec into the grid (replace or append).

        Interval blocks render as columns or compact rows under the same
        threshold as before.

        Returns ``False`` if the user cancelled a large-load prompt.
        """
        was_tracking = self._tracking_edits
        self._tracking_edits = False
        try:
            return self._populate_grid_from_spec_impl(
                spec, source_name=source_name, append=append
            )
        finally:
            self._tracking_edits = was_tracking

    def _populate_grid_from_spec_impl(
        self,
        spec: dict,
        source_name: str | None = None,
        *,
        append: bool = False,
    ) -> bool:
        """Load a v2 block-based stimulus spec into the grid (replace or append).

        Interval blocks render as columns or compact rows under the same
        threshold as before.

        Returns ``False`` if the user cancelled a large-load prompt.
        """
        blocks = spec.get("blocks")
        if not isinstance(blocks, list) or not blocks:
            raise ValueError("Stimulus spec has no blocks.")

        custom_spectra = self._parent._custom_spectra

        manual_blocks = [b for b in blocks if b.get("type") == "interval"]

        for i, block in enumerate(manual_blocks):
            spectrum = str(block.get("spectrum_ref", "")).strip()
            if not spectrum:
                raise ValueError(f"Block {i + 1}: missing spectrum name.")
            if spectrum in _BUILTIN_STIMULUS_SPECTRA:
                continue
            if spectrum in custom_spectra:
                continue
            raise ValueError(
                f"Block {i + 1}: spectrum {spectrum!r} is not a built-in and "
                "is not present in the spectrum library."
            )

        n_manual = len(manual_blocks)
        show_interval_widgets = True
        if n_manual > self._DEFER_REFRESH_INTERVAL_THRESHOLD:
            widget_choice = self._prompt_large_load_widget_display(n_manual)
            if widget_choice is None:
                return False
            show_interval_widgets = widget_choice

        if not append:
            self._teardown_interval_widgets()
            self._intervals.clear()
            self._selected_col = None
            self._measured_col_width = None
            self._measured_col_pitch = None

        defer_ui_refresh = (
            show_interval_widgets
            and n_manual > self._DEFER_REFRESH_INTERVAL_THRESHOLD
        )

        for block in blocks:
            kind = block.get("type")
            if kind == "interval":
                spectrum = str(block["spectrum_ref"]).strip()
                intensity = float(block["intensity"])
                duration = float(block["duration"])
                if show_interval_widgets:
                    self._add_interval(defer_ui_refresh=defer_ui_refresh)
                    interval = self._intervals[-1]
                    interval["spectrum_var"].set(spectrum)
                    if spectrum == "Dark":
                        interval["intensity_var"].set("0")
                    else:
                        interval["intensity_var"].set(
                            self._format_intensity_for_display(float(intensity))
                        )
                    interval["duration_var"].set(f"{duration:g}")
                    self._apply_dark_intensity_lock(interval, spectrum)
                else:
                    self._add_interval_data(spectrum, intensity, duration)
            else:
                raise ValueError(f"Unknown block type: {kind!r}")

        self._loaded_source_name = source_name
        self._update_compact_grid_display()
        self._refresh_interval_ui()
        self._render_preview(None)
        # Replace-load matches library/manuscript content; append is an unsaved edit.
        if append:
            self._stimulus_dirty = True
        else:
            self._mark_stimulus_clean()
        return True


    def _on_close(self) -> None:
        if self._has_unsaved_stimulus_changes():
            discard = messagebox.askyesno(
                "Unsaved changes",
                "The stimulus has unsaved changes.\n\n"
                "Close the Stimulus Builder and discard them?",
                parent=self,
            )
            if not discard:
                return
        try:
            plt.close(self._monitor_fig)
        except Exception:
            pass
        try:
            plt.close(self._preview_fig)
        except Exception:
            pass
        self._parent._stimulus_creator_dialog = None
        self.destroy()


# kept for backward compatibility — alias points at new class
StimulusCreatorDialog = StimulusBuilderDialog


class AutoSaveOptionsDialog(tk.Toplevel):
    """Checkboxes for optional auto-save of figure and CSV data after each run."""

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Auto-save options")
        self._parent = parent
        self.transient(parent)

        outer = ttk.Frame(self, padding=12)
        outer.pack(fill=tk.BOTH, expand=True)
        outer.columnconfigure(0, weight=1)

        r = 0
        self._cb_img = ttk.Checkbutton(
            outer,
            text="Save figure after each run",
            variable=parent._autosave_images,
        )
        self._cb_img.grid(row=r, column=0, sticky=tk.W, pady=(0, 4))
        r += 1

        self._img_frame = ttk.Frame(outer)
        ttk.Label(self._img_frame, text="Save folder (figures):").grid(
            row=0, column=0, columnspan=2, sticky=tk.W, pady=(0, 2)
        )
        self._img_path_lbl = ttk.Label(
            self._img_frame,
            text=_format_autosave_path_for_display(parent._autosave_images_dir),
            wraplength=480,
        )
        self._img_path_lbl.grid(row=1, column=0, columnspan=2, sticky=tk.W)
        img_btn = ttk.Frame(self._img_frame)
        img_btn.grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(4, 0))
        ttk.Button(
            img_btn, text="Choose folder…", command=self._on_browse_images
        ).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(
            img_btn, text="Reset to default", command=self._on_reset_images_dir
        ).pack(side=tk.LEFT)
        self._img_frame.grid(row=r, column=0, sticky=tk.EW, pady=(0, 8))
        r += 1

        self._cb_data = ttk.Checkbutton(
            outer,
            text="Save simulation data (CSV) after each run",
            variable=parent._autosave_data,
        )
        self._cb_data.grid(row=r, column=0, sticky=tk.W, pady=(0, 4))
        r += 1

        self._data_frame = ttk.Frame(outer)
        ttk.Label(self._data_frame, text="Save folder (CSV data):").grid(
            row=0, column=0, columnspan=2, sticky=tk.W, pady=(0, 2)
        )
        self._data_path_lbl = ttk.Label(
            self._data_frame,
            text=_format_autosave_path_for_display(parent._autosave_data_dir),
            wraplength=480,
        )
        self._data_path_lbl.grid(row=1, column=0, columnspan=2, sticky=tk.W)
        data_btn = ttk.Frame(self._data_frame)
        data_btn.grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(4, 0))
        ttk.Button(
            data_btn, text="Choose folder…", command=self._on_browse_data
        ).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(
            data_btn, text="Reset to default", command=self._on_reset_data_dir
        ).pack(side=tk.LEFT)
        self._data_frame.grid(row=r, column=0, sticky=tk.EW, pady=(0, 8))
        r += 1

        btn_row = ttk.Frame(outer, padding=(0, 4, 0, 0))
        btn_row.grid(row=r, column=0, sticky=tk.E)
        ttk.Button(btn_row, text="Close", command=self._on_close).pack(side=tk.RIGHT)

        self._trace_id_autosave_img: str = parent._autosave_images.trace_add(
            "write", self._refresh_path_sections
        )
        self._trace_id_autosave_data: str = parent._autosave_data.trace_add(
            "write", self._refresh_path_sections
        )
        self._refresh_path_sections()

        self.minsize(520, 200)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _initial_dir_for_dialog(self, p: Path) -> str:
        if p.exists() and p.is_dir():
            return str(p.resolve())
        if p.parent.exists() and p.parent.is_dir():
            return str(p.parent.resolve())
        return str(Path.cwd())

    def _refresh_path_sections(self, *_args: object) -> None:
        p = self._parent
        self._img_path_lbl.config(
            text=_format_autosave_path_for_display(p._autosave_images_dir)
        )
        self._data_path_lbl.config(
            text=_format_autosave_path_for_display(p._autosave_data_dir)
        )
        if p._autosave_images.get():
            self._img_frame.grid()
        else:
            self._img_frame.grid_remove()
        if p._autosave_data.get():
            self._data_frame.grid()
        else:
            self._data_frame.grid_remove()

    def _on_browse_images(self) -> None:
        start = self._initial_dir_for_dialog(self._parent._autosave_images_dir)
        d = filedialog.askdirectory(
            parent=self,
            title="Select folder for auto-saved figures",
            initialdir=start,
        )
        if d:
            self._parent._autosave_images_dir = Path(d)
            self._refresh_path_sections()

    def _on_browse_data(self) -> None:
        start = self._initial_dir_for_dialog(self._parent._autosave_data_dir)
        d = filedialog.askdirectory(
            parent=self,
            title="Select folder for auto-saved CSV data",
            initialdir=start,
        )
        if d:
            self._parent._autosave_data_dir = Path(d)
            self._refresh_path_sections()

    def _on_reset_images_dir(self) -> None:
        self._parent._autosave_images_dir = Path(DEFAULT_AUTOSAVE_IMAGES_DIR)
        self._refresh_path_sections()

    def _on_reset_data_dir(self) -> None:
        self._parent._autosave_data_dir = Path(DEFAULT_AUTOSAVE_DATA_DIR)
        self._refresh_path_sections()

    def _on_close(self) -> None:
        p = self._parent
        try:
            p._autosave_images.trace_remove("write", self._trace_id_autosave_img)
        except (tk.TclError, ValueError):
            pass
        try:
            p._autosave_data.trace_remove("write", self._trace_id_autosave_data)
        except (tk.TclError, ValueError):
            pass
        p._autosave_options_dialog = None
        self.destroy()


_GITHUB_REPO_URL = "https://github.com/Do-Laboratory/melanopsin-model"


class StimulusDeleteProgressDialog(tk.Toplevel):
    """Progress window shown while a stimulus is deleted from disk."""

    def __init__(self, parent: tk.Misc, label: str) -> None:
        super().__init__(parent)
        self.title("Deleting stimulus")
        self.transient(parent)
        self.resizable(False, False)
        self.protocol("WM_DELETE_WINDOW", lambda: None)

        outer = ttk.Frame(self, padding=12)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(
            outer,
            text=f"Deleting {label!r}…",
            wraplength=320,
        ).pack(anchor=tk.W)
        self._status = ttk.Label(outer, text="Starting…")
        self._status.pack(anchor=tk.W, pady=(8, 4))
        self._progress = ttk.Progressbar(
            outer, orient=tk.HORIZONTAL, mode="determinate", maximum=100, value=0
        )
        self._progress.pack(fill=tk.X, pady=(0, 4))
        self.update_idletasks()
        self.grab_set()

    def set_progress(self, message: str, percent: int) -> None:
        self._status.config(text=message)
        self._progress.configure(value=max(0, min(100, int(percent))))
        self.update_idletasks()


def _names_matching_search(names: list[str], query: str) -> list[str]:
    """Return names whose text contains ``query`` (case-insensitive)."""
    needle = query.strip().casefold()
    if not needle:
        return list(names)
    return [name for name in names if needle in name.casefold()]


def _listbox_fill_filtered(
    listbox: tk.Listbox,
    names: list[str],
    *,
    select_name: str | None = None,
) -> str | None:
    """Repopulate a listbox and optionally select one visible name."""
    listbox.delete(0, tk.END)
    for name in names:
        listbox.insert(tk.END, name)
    if select_name is not None and select_name in names:
        idx = names.index(select_name)
        listbox.selection_set(idx)
        listbox.see(idx)
        return select_name
    return None


def _install_listbox_search_row(parent: ttk.Frame) -> tk.StringVar:
    """Pack a labeled search entry above a list and return its StringVar."""
    search_row = ttk.Frame(parent)
    search_row.pack(side=tk.TOP, fill=tk.X, pady=(0, 4))
    ttk.Label(search_row, text="Search:").pack(side=tk.LEFT)
    search_var = tk.StringVar(master=parent)
    ttk.Entry(search_row, textvariable=search_var).pack(
        side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0)
    )
    return search_var


class StimulusLibraryPopup(tk.Toplevel):
    """Popup showing the library of known stimuli with a per-interval preview."""

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Stimulus Library")
        self._parent = parent
        self.transient(parent)

        self._preview_frames: list[dict] = []
        self._preview_ymax: float | None = None
        self._current_label: str | None = None
        self._preview_deferred = False

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        btn_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        btn_row.pack(side=tk.BOTTOM, fill=tk.X)
        self._delete_btn = ttk.Button(
            btn_row, text="Delete stimulus", command=self._on_delete_stimulus,
            state=tk.DISABLED,
        )
        self._delete_btn.pack(side=tk.LEFT)
        ttk.Button(
            btn_row,
            text="Import Stimulus",
            command=self._parent._on_file_import_stimulus,
        ).pack(side=tk.LEFT, padx=(8, 0))
        ttk.Button(
            btn_row,
            text="Make Custom Stimulus",
            command=self._parent._open_stimulus_creator,
        ).pack(side=tk.LEFT, padx=(8, 0))

        panes = ttk.Frame(outer)
        panes.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        left_col = ttk.Frame(panes)
        left_col.pack(side=tk.LEFT, fill=tk.BOTH)
        self._all_labels: list[str] = []
        self._search_var = _install_listbox_search_row(left_col)
        self._search_var.trace_add("write", lambda *_: self._apply_search_filter())

        list_frame = ttk.Frame(left_col)
        list_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self._listbox = tk.Listbox(list_frame, height=12, exportselection=False, width=28)
        self._listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self._listbox.yview
        )
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self._listbox.configure(yscrollcommand=scrollbar.set)

        plot_frame = ttk.Frame(panes, padding=(8, 0, 0, 0))
        plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        preview_controls = ttk.Frame(plot_frame)
        preview_controls.pack(side=tk.BOTTOM, fill=tk.X, pady=(6, 0))

        self._preview_btn_frame = ttk.Frame(preview_controls)
        self._preview_btn = ttk.Button(
            self._preview_btn_frame,
            text="Preview stimulus",
            command=self._on_preview_stimulus,
        )

        self._scrubber_frame = ttk.Frame(preview_controls)
        self._scrubber_var = tk.IntVar(master=self, value=0)
        self._scrubber_scale = ttk.Scale(
            self._scrubber_frame,
            from_=0,
            to=0,
            orient=tk.HORIZONTAL,
            variable=self._scrubber_var,
            command=self._on_scrubber_moved,
        )
        self._scrubber_scale.pack(side=tk.TOP, fill=tk.X)
        self._scrubber_label = ttk.Label(self._scrubber_frame, text="", anchor=tk.W)
        self._scrubber_label.pack(side=tk.TOP, fill=tk.X, pady=(4, 0))

        # Do not enable Figure(tight_layout=True): each redraw ratchets subplot
        # margins inward and progressively squishes axes/labels when browsing.
        self._fig = plt.Figure(figsize=(5, 3.2))
        self._ax = self._fig.add_subplot(111)
        self._preview_title_text = ""
        self._preview_title_fontsize = 11
        self._preview_title_wrap_width = 0
        self._preview_canvas_size = (0, 0)
        self._apply_preview_layout()
        self._canvas = FigureCanvasTkAgg(self._fig, master=plot_frame)
        canvas_widget = self._canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True)
        # Use add="+": FigureCanvasTkAgg already binds Configure->resize to grow
        # its PhotoImage. A plain bind would replace that and leave white space.
        canvas_widget.bind(
            "<Configure>", self._on_preview_canvas_configure, add="+"
        )

        self._listbox.bind("<<ListboxSelect>>", self._on_selection_changed)

        self._refresh_list()

        self.minsize(720, 360)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _preview_title_char_width(self) -> int:
        """Estimate wrap width from the current preview pane size."""
        fontsize = max(float(self._preview_title_fontsize), 1.0)
        width_px = 0
        try:
            width_px = int(self._canvas.get_tk_widget().winfo_width())
        except Exception:
            width_px = 0
        if width_px <= 1:
            width_px = int(float(self._fig.get_size_inches()[0]) * float(self._fig.dpi))
        # ~0.6 * fontsize pixels per character; leave side padding for axis labels.
        usable_px = max(80, width_px - 40)
        return max(18, int(usable_px / max(6.0, 0.6 * fontsize)))

    def _set_preview_title(self, title: str, *, fontsize: float = 11) -> None:
        """Set a wrapped title and reserve top margin so it stays visible."""
        self._preview_title_text = str(title or "")
        self._preview_title_fontsize = float(fontsize)
        self._apply_preview_title()

    def _apply_preview_title(self) -> None:
        """Wrap the current title to the figure width and adjust subplot top."""
        width_chars = self._preview_title_char_width()
        self._preview_title_wrap_width = width_chars
        wrapped = self._wrap_preview_text(
            self._preview_title_text, width=width_chars
        )
        n_lines = wrapped.count("\n") + 1 if wrapped else 1
        # Leave room for each title line without collapsing the plot area.
        top = min(0.90, max(0.68, 0.96 - 0.07 * n_lines))
        self._fig.set_tight_layout(False)
        self._fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=top)
        self._ax.set_title(wrapped, fontsize=self._preview_title_fontsize)

    def _apply_preview_layout(self) -> None:
        """Keep stable axis margins across stimulus selection changes."""
        if self._preview_title_text:
            self._apply_preview_title()
            return
        self._fig.set_tight_layout(False)
        self._fig.subplots_adjust(left=0.16, right=0.98, bottom=0.16, top=0.86)

    def _on_preview_canvas_configure(self, event: tk.Event) -> None:
        """After matplotlib resizes the PhotoImage, re-wrap the title to fit."""
        w, h = int(event.width), int(event.height)
        if w < 50 or h < 50:
            return
        if (w, h) == self._preview_canvas_size:
            return
        self._preview_canvas_size = (w, h)
        self._apply_preview_layout()
        self._canvas.draw_idle()

    def _apply_search_filter(self, *, preferred: str | None = None) -> None:
        previous = preferred if preferred is not None else self._selected_stimulus_label()
        filtered = _names_matching_search(self._all_labels, self._search_var.get())
        if previous is not None and previous not in filtered:
            if previous not in self._all_labels:
                previous = None
        selected = _listbox_fill_filtered(
            self._listbox, filtered, select_name=previous
        )
        self._prepare_preview_for_label(selected)
        self._update_delete_button_state()

    def _refresh_list(self) -> None:
        """Reload library from disk and refresh the stimulus list."""
        previous = self._selected_stimulus_label()
        if previous is None:
            current = self._parent._selected_stimulus_var.get()
            if current in self._parent._known_stimuli:
                previous = current
        self._parent._reload_stimulus_library_from_disk()
        self._all_labels = list(self._parent._known_stimuli)
        self._apply_search_filter(preferred=previous)

    def _selected_stimulus_label(self) -> str | None:
        sel = self._listbox.curselection()
        if not sel:
            return None
        return str(self._listbox.get(sel[0]))

    @staticmethod
    def _is_manuscript_stimulus(label: str) -> bool:
        return label in LABELS

    def _update_delete_button_state(self) -> None:
        label = self._selected_stimulus_label()
        can_delete = label is not None and not self._is_manuscript_stimulus(label)
        self._delete_btn.config(state=tk.NORMAL if can_delete else tk.DISABLED)

    def _hide_preview_controls(self) -> None:
        self._preview_btn_frame.pack_forget()
        self._preview_btn.pack_forget()
        self._scrubber_frame.pack_forget()

    def _show_preview_button(self) -> None:
        self._scrubber_frame.pack_forget()
        self._preview_btn.pack(side=tk.LEFT)
        self._preview_btn_frame.pack(fill=tk.X)

    def _show_scrubber(self, n_frames: int) -> None:
        self._preview_btn_frame.pack_forget()
        self._preview_btn.pack_forget()
        max_idx = max(0, n_frames - 1)
        self._scrubber_scale.configure(from_=0, to=max_idx)
        self._scrubber_var.set(0)
        self._scrubber_frame.pack(fill=tk.X)

    def _stimulus_spec_for_label(self, label: str) -> dict:
        if label in self._parent._custom_stimuli:
            return self._parent._custom_stimuli[label]
        if label in LABELS:
            return self._parent._build_stimulus_spec_from_label(label)
        raise ValueError(f"Unknown stimulus {label!r}.")

    def _interval_count_for_label(self, label: str) -> int:
        if label in self._parent._custom_stimuli:
            return library_storage.total_intervals(self._parent._custom_stimuli[label])
        if label in LABELS:
            stim = _build_stimulus(label, self._parent._env)
            timings = np.asarray(stim["timings"], dtype=float).reshape(-1)
            return int(timings.size // 2)
        return 0

    def _build_preview_frames(self, label: str) -> list[dict]:
        protocol = self._parent._build_run_protocol(label)
        intensities = np.asarray(protocol["intensities"], dtype=float)
        wlen = np.asarray(protocol["wlen"], dtype=float)
        timings = np.asarray(protocol["timings"], dtype=float).reshape(-1)
        n = int(intensities.shape[0])
        if n == 0:
            return []

        try:
            names = _interval_names_from_spec(self._stimulus_spec_for_label(label))
        except Exception:
            names = []
        if len(names) != n:
            names = [f"Interval {i + 1}" for i in range(n)]

        frames: list[dict] = []
        for i in range(n):
            t0 = float(timings[2 * i])
            t1 = float(timings[2 * i + 1])
            duration = max(0.0, t1 - t0)
            frames.append(
                {
                    "wlen": np.asarray(wlen[i], dtype=float),
                    "intensity": np.asarray(intensities[i], dtype=float),
                    "caption": f"Interval {i + 1}/{n} · {names[i]} · {duration:g} s",
                    "title": f"Stimulus: {label}",
                }
            )
        return frames

    @staticmethod
    def _preview_peak_intensity(frames: list[dict]) -> float:
        """Return the maximum spectral intensity across all preview intervals."""
        peak = 0.0
        for frame in frames:
            intensity = np.asarray(frame.get("intensity", []), dtype=float)
            if intensity.size:
                peak = max(peak, float(np.max(intensity)))
        return peak if peak > 0 else 1.0

    @staticmethod
    def _wrap_preview_text(text: str, width: int = 40) -> str:
        """Wrap placeholder copy so it stays readable inside the preview axes."""
        lines: list[str] = []
        for paragraph in text.splitlines():
            if not paragraph.strip():
                lines.append("")
            else:
                lines.append(textwrap.fill(paragraph, width=width))
        return "\n".join(lines)

    def _render_preview_placeholder(
        self, title: str, message: str | None = None
    ) -> None:
        ax = self._ax
        ax.clear()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Intensity (photons/\u00b5m\u00b2/nm/s)")
        ax.set_xlim(100, 900)
        ax.ticklabel_format(
            axis="y", style="scientific", scilimits=(0, 0), useMathText=True
        )
        if message:
            ax.text(
                0.5,
                0.5,
                self._wrap_preview_text(message),
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=10,
                linespacing=1.4,
            )
        # Interval details live on the scrubber label; keep the axes title short
        # and wrap it to the current pane width so it stays visible on resize.
        self._set_preview_title(title, fontsize=11)
        self._canvas.draw_idle()

    def _render_preview_frame(self, index: int) -> None:
        if not self._preview_frames:
            return
        index = max(0, min(index, len(self._preview_frames) - 1))
        frame = self._preview_frames[index]
        ax = self._ax
        ax.clear()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Intensity (photons/\u00b5m\u00b2/nm/s)")
        ax.set_xlim(100, 900)
        ax.ticklabel_format(
            axis="y", style="scientific", scilimits=(0, 0), useMathText=True
        )
        wlen = frame["wlen"]
        intensity = frame["intensity"]
        if wlen.size and intensity.size:
            ax.plot(wlen, intensity, color="k")
        ymax = self._preview_ymax
        if ymax is not None and ymax > 0:
            ax.set_ylim(0.0, ymax)
        self._scrubber_label.config(text=frame["caption"])
        self._set_preview_title(frame["title"], fontsize=11)
        self._canvas.draw_idle()

    def _prepare_preview_for_label(self, label: str | None) -> None:
        self._preview_frames = []
        self._preview_ymax = None
        self._preview_deferred = False
        self._current_label = label
        self._hide_preview_controls()

        if label is None:
            self._render_preview_placeholder("Select a stimulus")
            return

        if self._parent._env is None:
            self._render_preview_placeholder(
                f"Stimulus: {label}", "Spectra not loaded yet."
            )
            return

        try:
            n_intervals = self._interval_count_for_label(label)
        except Exception as exc:
            self._render_preview_placeholder(f"Stimulus: {label}", str(exc))
            return

        if n_intervals <= 0:
            self._render_preview_placeholder(f"Stimulus: {label}", "No intervals.")
            return

        if n_intervals >= _STIMULUS_PREVIEW_AUTO_MAX_INTERVALS:
            self._preview_deferred = True
            self._show_preview_button()
            self._render_preview_placeholder(
                f"Stimulus: {label}",
                (
                    f"This stimulus has {n_intervals:,} intervals,\n"
                    "which is too many to load automatically.\n\n"
                    "Click Preview stimulus below to view spectra."
                ),
            )
            return

        self._load_and_show_preview(label)

    def _load_and_show_preview(self, label: str) -> None:
        try:
            frames = self._build_preview_frames(label)
        except Exception as exc:
            self._hide_preview_controls()
            self._render_preview_placeholder(f"Stimulus: {label}", str(exc))
            return
        if not frames:
            self._hide_preview_controls()
            self._render_preview_placeholder(f"Stimulus: {label}", "No intervals to preview.")
            return

        self._preview_frames = frames
        self._preview_ymax = self._preview_peak_intensity(frames)
        self._preview_deferred = False
        self._show_scrubber(len(frames))
        self._scrubber_var.set(0)
        self._render_preview_frame(0)

    def _on_preview_stimulus(self) -> None:
        label = self._current_label
        if label is None or not self._preview_deferred:
            return
        self._load_and_show_preview(label)

    def _on_scrubber_moved(self, _value: str | None = None) -> None:
        if not self._preview_frames:
            return
        idx = int(round(float(self._scrubber_var.get())))
        idx = max(0, min(idx, len(self._preview_frames) - 1))
        if int(self._scrubber_var.get()) != idx:
            self._scrubber_var.set(idx)
        self._render_preview_frame(idx)

    def _on_selection_changed(self, _event=None) -> None:
        sel = self._listbox.curselection()
        if not sel:
            self._prepare_preview_for_label(None)
            self._update_delete_button_state()
            return
        label = str(self._listbox.get(sel[0]))
        self._prepare_preview_for_label(label)
        self._update_delete_button_state()

    def _on_delete_stimulus(self) -> None:
        sel = self._listbox.curselection()
        if not sel:
            messagebox.showinfo(
                "Delete stimulus", "Select a stimulus to delete.", parent=self
            )
            return
        idx = sel[0]
        label = str(self._listbox.get(idx))
        if self._is_manuscript_stimulus(label):
            messagebox.showinfo(
                "Delete stimulus",
                "Manuscript stimuli cannot be deleted.",
                parent=self,
            )
            return
        ok = messagebox.askyesno(
            "Delete stimulus",
            f"Delete {label!r} from the stimulus library?\n\n"
            "This permanently removes it from the user stimulus library.",
            parent=self,
        )
        if not ok:
            return
        self._parent.start_remove_known_stimulus(
            label, parent=self, on_success=self._refresh_list
        )

    def _on_close(self) -> None:
        try:
            plt.close(self._fig)
        except Exception:
            pass
        self._parent._stimulus_library_popup = None
        self.destroy()


class PredictionComparisonDialog(tk.Toplevel):
    """Compare saved prediction runs by overlaying selected fields."""

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Prediction Comparison")
        self._parent = parent
        self.transient(parent)
        self._runs: list[dict] = []

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        panes = ttk.Frame(outer)
        panes.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(panes)
        left.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(left, text="Loaded runs").pack(anchor=tk.W)
        run_frame = ttk.Frame(left)
        run_frame.pack(fill=tk.BOTH, expand=False, pady=(2, 8))
        self._run_list = tk.Listbox(run_frame, exportselection=False, width=34, height=10)
        self._run_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        run_sb = ttk.Scrollbar(run_frame, orient=tk.VERTICAL, command=self._run_list.yview)
        run_sb.pack(side=tk.RIGHT, fill=tk.Y)
        self._run_list.configure(yscrollcommand=run_sb.set)

        run_btns = ttk.Frame(left)
        run_btns.pack(fill=tk.X, pady=(0, 8))
        ttk.Button(run_btns, text="Add files…", command=self._on_add_files).pack(
            side=tk.LEFT
        )
        ttk.Button(run_btns, text="Load current run", command=self._on_load_current_run).pack(
            side=tk.LEFT, padx=(6, 0)
        )
        ttk.Button(
            run_btns, text="Remove selected", command=self._on_remove_selected
        ).pack(side=tk.LEFT, padx=(6, 0))
        ttk.Button(run_btns, text="Clear", command=self._on_clear).pack(
            side=tk.LEFT, padx=(6, 0)
        )

        ttk.Label(left, text="Fields to compare").pack(anchor=tk.W)
        field_frame = ttk.Frame(left)
        field_frame.pack(fill=tk.BOTH, expand=True, pady=(2, 0))
        self._field_list = tk.Listbox(
            field_frame, selectmode=tk.EXTENDED, exportselection=False, width=34, height=12
        )
        self._field_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        field_sb = ttk.Scrollbar(field_frame, orient=tk.VERTICAL, command=self._field_list.yview)
        field_sb.pack(side=tk.RIGHT, fill=tk.Y)
        self._field_list.configure(yscrollcommand=field_sb.set)

        right = ttk.Frame(panes, padding=(8, 0, 0, 0))
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._fig = plt.Figure(figsize=(8, 5), tight_layout=True)
        self._canvas = FigureCanvasTkAgg(self._fig, master=right)
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._render_placeholder("Add exported prediction files to compare.")

        bottom = ttk.Frame(outer, padding=(0, 8, 0, 0))
        bottom.pack(fill=tk.X)
        ttk.Button(bottom, text="Plot", command=self._on_plot).pack(side=tk.LEFT)
        ttk.Button(bottom, text="Close", command=self._on_close).pack(side=tk.RIGHT)

        self.minsize(980, 540)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _render_placeholder(self, text: str) -> None:
        self._fig.clear()
        ax = self._fig.add_subplot(111)
        ax.set_axis_off()
        ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes)
        self._canvas.draw_idle()

    def _common_fields(self) -> list[str]:
        if not self._runs:
            return []
        sets = [set(run["series"].keys()) for run in self._runs]
        common = set.intersection(*sets) if sets else set()
        return sorted(common)

    def _refresh_run_list(self) -> None:
        self._run_list.delete(0, tk.END)
        for run in self._runs:
            self._run_list.insert(
                tk.END, f"{run['label']} ({Path(run['source_path']).name})"
            )

    def _refresh_field_list(self) -> None:
        prior = {str(self._field_list.get(i)) for i in self._field_list.curselection()}
        fields = self._common_fields()
        self._field_list.delete(0, tk.END)
        for f in fields:
            self._field_list.insert(tk.END, f)
        for i, f in enumerate(fields):
            if f in prior:
                self._field_list.selection_set(i)

    def _on_add_files(self) -> None:
        initial_dir = DEFAULT_AUTOSAVE_DATA_DIR
        paths = filedialog.askopenfilenames(
            parent=self,
            title="Select saved prediction files",
            initialdir=str(initial_dir) if initial_dir.is_dir() else None,
            filetypes=[
                ("Prediction exports", "*.npz *.csv"),
                ("NPZ files", "*.npz"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )
        if not paths:
            return
        added = 0
        for path in paths:
            try:
                run = _load_saved_prediction(str(path))
            except Exception as exc:
                messagebox.showerror(
                    "Prediction compare",
                    f"Could not load {Path(path).name}:\n{exc}",
                    parent=self,
                )
                continue
            self._runs.append(run)
            added += 1
        if added == 0:
            return
        self._refresh_run_list()
        self._refresh_field_list()
        fields = self._common_fields()
        if not fields:
            self._render_placeholder(
                "No common numeric time-aligned fields across loaded runs."
            )
        else:
            self._render_placeholder("Select one or more fields, then click Plot.")

    def _on_remove_selected(self) -> None:
        sel = list(self._run_list.curselection())
        if not sel:
            return
        for idx in sorted(sel, reverse=True):
            self._runs.pop(idx)
        self._refresh_run_list()
        self._refresh_field_list()
        if not self._runs:
            self._render_placeholder("Add exported prediction files to compare.")
        elif not self._common_fields():
            self._render_placeholder(
                "No common numeric time-aligned fields across loaded runs."
            )

    def _on_load_current_run(self) -> None:
        if self._parent._last_result is None:
            messagebox.showinfo(
                "Prediction compare",
                "No completed run is available yet. Run the model first.",
                parent=self,
            )
            return
        try:
            run = _normalize_in_memory_result_for_compare(
                self._parent._last_result,
                self._parent._last_label,
            )
        except Exception as exc:
            messagebox.showerror(
                "Prediction compare",
                f"Could not load current run:\n{exc}",
                parent=self,
            )
            return
        self._runs.append(run)
        self._refresh_run_list()
        self._refresh_field_list()
        fields = self._common_fields()
        if not fields:
            self._render_placeholder(
                "No common numeric time-aligned fields across loaded runs."
            )
        else:
            self._render_placeholder("Select one or more fields, then click Plot.")

    def _on_clear(self) -> None:
        self._runs.clear()
        self._refresh_run_list()
        self._refresh_field_list()
        self._render_placeholder("Add exported prediction files to compare.")

    def _selected_fields(self) -> list[str]:
        return [str(self._field_list.get(i)) for i in self._field_list.curselection()]

    def _on_plot(self) -> None:
        if len(self._runs) < 2:
            messagebox.showinfo(
                "Prediction compare",
                "Load at least two runs to compare.",
                parent=self,
            )
            return
        fields = self._selected_fields()
        if not fields:
            messagebox.showinfo(
                "Prediction compare",
                "Select at least one field to plot.",
                parent=self,
            )
            return

        self._fig.clear()
        axes = self._fig.subplots(
            len(fields), 1, sharex=False, squeeze=False
        ).reshape(-1)
        for ax, field in zip(axes, fields):
            plotted = 0
            for run in self._runs:
                series = run["series"]
                if field not in series:
                    continue
                ax.plot(run["time_s"], series[field], label=run["label"])
                plotted += 1
            ax.set_ylabel(field)
            if plotted == 0:
                ax.text(
                    0.5,
                    0.5,
                    "No data for this field.",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
            else:
                ax.legend(loc="best", fontsize=8)
        axes[-1].set_xlabel("time (s)")
        self._fig.tight_layout()
        self._canvas.draw_idle()

    def _on_close(self) -> None:
        try:
            plt.close(self._fig)
        except Exception:
            pass
        self._parent._prediction_compare_dialog = None
        self.destroy()


class SamplingRateDialog(tk.Toplevel):
    """Modal prompt for the sampling rate of a dataset with no time column.

    The chosen rate in Hz is available in ``self.result`` after the dialog
    closes (``None`` when cancelled). The user may enter either Hz or a
    sample interval in seconds.
    """

    _MODE_HZ = "hz"
    _MODE_DT = "dt"

    def __init__(self, parent, default_hz: float = DATA_COMPARE_SAMPLE_RATE_HZ):
        super().__init__(parent)
        self.result: float | None = None
        self.title("Sampling rate")
        self.transient(parent)
        self.resizable(False, False)

        self._mode_var = tk.StringVar(master=self, value=self._MODE_HZ)
        self._value_var = tk.StringVar(master=self, value=f"{float(default_hz):g}")

        frame = ttk.Frame(self, padding=16)
        frame.grid(row=0, column=0, sticky=tk.NSEW)

        ttk.Label(
            frame,
            text=(
                "No time column was found in this dataset.\n"
                "Enter the sampling rate so a time base can be built."
            ),
            justify=tk.LEFT,
        ).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=(0, 12))

        ttk.Radiobutton(
            frame,
            text="Sampling rate (Hz)",
            value=self._MODE_HZ,
            variable=self._mode_var,
        ).grid(row=1, column=0, sticky=tk.W)
        ttk.Radiobutton(
            frame,
            text="Sample interval (s)",
            value=self._MODE_DT,
            variable=self._mode_var,
        ).grid(row=2, column=0, sticky=tk.W)

        ttk.Entry(frame, textvariable=self._value_var, width=18).grid(
            row=1, column=1, rowspan=2, sticky=tk.W, padx=(12, 0)
        )

        button_row = ttk.Frame(frame)
        button_row.grid(row=3, column=0, columnspan=2, sticky=tk.E, pady=(16, 0))
        ttk.Button(button_row, text="OK", command=self._on_ok).grid(
            row=0, column=0, padx=(0, 6)
        )
        ttk.Button(button_row, text="Cancel", command=self._on_cancel).grid(
            row=0, column=1
        )

        self.protocol("WM_DELETE_WINDOW", self._on_cancel)
        self.bind("<Return>", lambda _e: self._on_ok())
        self.bind("<Escape>", lambda _e: self._on_cancel())
        self.grab_set()
        self.wait_window()

    def _on_ok(self) -> None:
        raw = self._value_var.get().strip()
        try:
            value = float(raw)
        except ValueError:
            messagebox.showwarning(
                "Sampling rate",
                "Enter a positive number.",
                parent=self,
            )
            return
        if not np.isfinite(value) or value <= 0:
            messagebox.showwarning(
                "Sampling rate",
                "Enter a positive number.",
                parent=self,
            )
            return
        if self._mode_var.get() == self._MODE_DT:
            rate_hz = 1.0 / value
        else:
            rate_hz = value
        self.result = float(rate_hz)
        self.destroy()

    def _on_cancel(self) -> None:
        self.result = None
        self.destroy()


def _prompt_sampling_rate(
    parent, default_hz: float = DATA_COMPARE_SAMPLE_RATE_HZ
) -> float | None:
    """Show the sampling-rate dialog; return Hz or None if cancelled."""
    return SamplingRateDialog(parent, default_hz=default_hz).result


class DataComparatorDialog(tk.Toplevel):
    """Compare an experimental data trace to a model run via ``plotModelVsData``.

    Datasets may include a ``time`` column (units inferred from the header) or
    a user-supplied sampling rate. Multiple data columns are supported via a
    column picker.
    """

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Data Comparator")
        self._parent = parent
        self.transient(parent)

        self._runs: list[dict] = []
        self._dataset_columns: dict[str, np.ndarray] | None = None
        self._dataset_time_s: np.ndarray | None = None
        self._dataset_rate_hz: float | None = None
        self._dataset_path: str | None = None
        self._dataset_time_source: str | None = None  # "column" | "rate" | None
        self._dataset_time_column: str | None = None
        self._dataset_time_unit_label: str | None = None
        self._dataset_time_offset_s: float = 0.0
        self._compare_fig = None
        self._fig_placeholder: plt.Figure | None = None

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        panes = ttk.Frame(outer)
        panes.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(panes)
        left.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(left, text="Model runs").pack(anchor=tk.W)
        run_frame = ttk.Frame(left)
        run_frame.pack(fill=tk.BOTH, expand=False, pady=(2, 8))
        self._run_list = tk.Listbox(
            run_frame, exportselection=False, width=36, height=8
        )
        self._run_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        run_sb = ttk.Scrollbar(
            run_frame, orient=tk.VERTICAL, command=self._run_list.yview
        )
        run_sb.pack(side=tk.RIGHT, fill=tk.Y)
        self._run_list.configure(yscrollcommand=run_sb.set)
        self._run_list.bind("<<ListboxSelect>>", self._on_run_selection_changed)

        run_btns = ttk.Frame(left)
        run_btns.pack(fill=tk.X, pady=(0, 8))
        ttk.Button(run_btns, text="Add saved run…", command=self._on_add_saved_run).pack(
            side=tk.LEFT
        )
        ttk.Button(
            run_btns, text="Load current run", command=self._on_load_current_run
        ).pack(side=tk.LEFT, padx=(6, 0))
        ttk.Button(
            run_btns, text="Remove selected", command=self._on_remove_selected_run
        ).pack(side=tk.LEFT, padx=(6, 0))
        ttk.Button(run_btns, text="Clear runs", command=self._on_clear_runs).pack(
            side=tk.LEFT, padx=(6, 0)
        )

        ttk.Label(left, text="Model series (prediction)").pack(anchor=tk.W)
        self._series_var = tk.StringVar(master=self, value="")
        self._series_combo = ttk.Combobox(
            left,
            textvariable=self._series_var,
            state="readonly",
            width=34,
        )
        self._series_combo.pack(fill=tk.X, pady=(2, 8))

        ds_frame = ttk.LabelFrame(left, text="Dataset", padding=6)
        ds_frame.pack(fill=tk.X, pady=(0, 8))
        ttk.Button(
            ds_frame, text="Load dataset…", command=self._on_load_dataset
        ).pack(anchor=tk.W)

        col_row = ttk.Frame(ds_frame)
        col_row.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(col_row, text="Data column").pack(side=tk.LEFT)
        self._data_col_var = tk.StringVar(master=self, value="")
        self._data_col_combo = ttk.Combobox(
            col_row,
            textvariable=self._data_col_var,
            state="disabled",
            width=22,
        )
        self._data_col_combo.pack(side=tk.LEFT, padx=(6, 0), fill=tk.X, expand=True)

        rate_row = ttk.Frame(ds_frame)
        rate_row.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(rate_row, text="Sampling rate (Hz)").pack(side=tk.LEFT)
        self._rate_var = tk.StringVar(master=self, value="")
        self._rate_entry = ttk.Entry(
            rate_row, textvariable=self._rate_var, width=12, state="disabled"
        )
        self._rate_entry.pack(side=tk.LEFT, padx=(6, 0))

        self._dataset_status = tk.StringVar(
            master=self, value="No dataset loaded."
        )
        ttk.Label(ds_frame, textvariable=self._dataset_status, wraplength=280).pack(
            anchor=tk.W, pady=(4, 0)
        )

        ttk.Button(left, text="Plot model vs data", command=self._on_plot).pack(
            anchor=tk.W
        )

        right = ttk.Frame(panes, padding=(10, 0, 0, 0))
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._plot_host = ttk.Frame(right)
        self._plot_host.pack(fill=tk.BOTH, expand=True)
        self._fig_placeholder = None
        self._render_placeholder("Add a model run and dataset, then click Plot.")

        bottom = ttk.Frame(outer, padding=(0, 8, 0, 0))
        bottom.pack(fill=tk.X)
        ttk.Button(bottom, text="Close", command=self._on_close).pack(side=tk.RIGHT)

        self.minsize(1020, 560)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _render_placeholder(self, text: str) -> None:
        self._close_compare_figure()
        for w in self._plot_host.winfo_children():
            w.destroy()
        if self._fig_placeholder is None:
            self._fig_placeholder = plt.Figure(figsize=(8, 5), tight_layout=True)
        self._fig_placeholder.clear()
        ax = self._fig_placeholder.add_subplot(111)
        ax.set_axis_off()
        ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes)
        self._canvas = FigureCanvasTkAgg(self._fig_placeholder, master=self._plot_host)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _refresh_run_list(self) -> None:
        self._run_list.delete(0, tk.END)
        for run in self._runs:
            src = run.get("source_path", "")
            disp = Path(str(src)).name if src else ""
            self._run_list.insert(tk.END, f"{run['label']} ({disp})")

    def _selected_run_index(self) -> int | None:
        sel = self._run_list.curselection()
        if not sel:
            return None
        return int(sel[0])

    def _model_for_selected_run(self) -> dict | None:
        idx = self._selected_run_index()
        if idx is None or idx < 0 or idx >= len(self._runs):
            return None
        return self._runs[idx]["model"]

    def _refresh_series_combo(self, preserve_selection: bool = False) -> None:
        prior = self._series_var.get().strip() if preserve_selection else ""
        model = self._model_for_selected_run()
        if model is None:
            self._series_combo.configure(values=())
            self._series_var.set("")
            return
        names = _aligned_series_field_names(model)
        self._series_combo.configure(values=names)
        if prior in names:
            self._series_var.set(prior)
        else:
            self._series_var.set(_default_data_compare_series_name(model))

    def _on_run_selection_changed(self, _event=None) -> None:
        self._refresh_series_combo(preserve_selection=False)

    def _append_run(self, label: str, source_path: str, raw: dict) -> None:
        model = _coerce_model_dict_for_plot_model_vs_data(raw)
        self._runs.append(
            {"label": label, "source_path": source_path, "model": model}
        )
        self._refresh_run_list()
        self._run_list.selection_clear(0, tk.END)
        self._run_list.selection_set(tk.END)
        self._run_list.see(tk.END)
        self._refresh_series_combo(preserve_selection=False)

    def _on_add_saved_run(self) -> None:
        initial_dir = DEFAULT_AUTOSAVE_DATA_DIR
        path = filedialog.askopenfilename(
            parent=self,
            title="Select saved model run",
            initialdir=str(initial_dir) if initial_dir.is_dir() else None,
            filetypes=[
                ("Prediction exports", "*.npz *.csv"),
                ("NPZ files", "*.npz"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        try:
            loaded = _load_saved_prediction(str(path))
            label = str(loaded.get("label", Path(path).stem))
            self._append_run(label, str(path), loaded)
        except Exception as exc:
            messagebox.showerror(
                "Data Comparator",
                f"Could not load run:\n{exc}",
                parent=self,
            )

    def _on_load_current_run(self) -> None:
        if self._parent._last_result is None:
            messagebox.showinfo(
                "Data Comparator",
                "No completed run yet. Click Run model in the main window first.",
                parent=self,
            )
            return
        try:
            label = (self._parent._last_label or "Current run").strip() or "Current run"
            self._append_run(label, "<current run>", self._parent._last_result)
        except Exception as exc:
            messagebox.showerror(
                "Data Comparator",
                f"Could not use current run:\n{exc}",
                parent=self,
            )

    def _on_remove_selected_run(self) -> None:
        idx = self._selected_run_index()
        if idx is None:
            return
        self._runs.pop(idx)
        self._refresh_run_list()
        self._series_combo.configure(values=())
        self._series_var.set("")
        if not self._runs:
            self._render_placeholder("Add a model run and dataset, then click Plot.")
        else:
            self._run_list.selection_set(min(idx, len(self._runs) - 1))
            self._refresh_series_combo(preserve_selection=False)

    def _on_clear_runs(self) -> None:
        self._runs.clear()
        self._refresh_run_list()
        self._series_combo.configure(values=())
        self._series_var.set("")
        self._render_placeholder("Add a model run and dataset, then click Plot.")

    def _update_dataset_status(self) -> None:
        if not self._dataset_columns or not self._dataset_path:
            self._dataset_status.set("No dataset loaded.")
            return
        name = Path(self._dataset_path).name
        n_cols = len(self._dataset_columns)
        col = self._data_col_var.get().strip()
        n_samples = (
            int(np.sum(np.isfinite(self._dataset_columns[col])))
            if col in self._dataset_columns
            else 0
        )
        if self._dataset_time_source == "column" and self._dataset_time_s is not None:
            finite_t = self._dataset_time_s[np.isfinite(self._dataset_time_s)]
            t0 = float(finite_t[0]) if finite_t.size else 0.0
            t1 = float(finite_t[-1]) if finite_t.size else 0.0
            unit = self._dataset_time_unit_label or "s"
            tcol = self._dataset_time_column or "time"
            offset_note = ""
            if abs(self._dataset_time_offset_s) > 1e-12:
                offset_note = f", file start {self._dataset_time_offset_s:g} s"
            self._dataset_status.set(
                f"{name}: {n_cols} column(s), {n_samples} samples; "
                f"time from {tcol!r} ({unit}) "
                f"{t0:.2f}–{t1:.2f} s{offset_note}"
            )
        elif self._dataset_rate_hz is not None:
            self._dataset_status.set(
                f"{name}: {n_cols} column(s), {n_samples} samples @ "
                f"{self._dataset_rate_hz:g} Hz (entered)"
            )
        else:
            self._dataset_status.set(
                f"{name}: {n_cols} column(s), {n_samples} samples"
            )

    def _on_load_dataset(self) -> None:
        path = filedialog.askopenfilename(
            parent=self,
            title="Select dataset file",
            filetypes=[
                ("Spreadsheets", "*.csv *.tsv *.txt *.xlsx *.xls"),
                ("CSV files", "*.csv"),
                ("Excel files", "*.xlsx *.xls"),
                ("TSV / text files", "*.tsv *.txt"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        try:
            loaded = _read_data_trace(str(path))
        except Exception as exc:
            messagebox.showerror(
                "Data Comparator",
                f"Could not load dataset:\n{exc}",
                parent=self,
            )
            return

        time_s = loaded.get("time_s")
        rate_hz: float | None = None
        time_source: str

        if time_s is None:
            rate = _prompt_sampling_rate(self, default_hz=DATA_COMPARE_SAMPLE_RATE_HZ)
            if rate is None:
                # Cancel: leave any previously loaded dataset intact.
                return
            rate_hz = float(rate)
            time_source = "rate"
        else:
            time_source = "column"
            finite_t = np.asarray(time_s, dtype=float)
            finite_t = finite_t[np.isfinite(finite_t)]
            if finite_t.size >= 2:
                dt = np.diff(finite_t)
                dt = dt[np.isfinite(dt) & (dt > 0)]
                if dt.size:
                    rate_hz = float(1.0 / np.mean(dt))

        # Commit loaded state only after a successful rate prompt (if needed).
        self._dataset_columns = loaded["columns"]
        self._dataset_time_s = (
            np.asarray(time_s, dtype=float) if time_s is not None else None
        )
        self._dataset_rate_hz = rate_hz
        self._dataset_path = str(path)
        self._dataset_time_source = time_source
        self._dataset_time_column = loaded.get("time_column")
        self._dataset_time_unit_label = loaded.get("time_unit_label")
        self._dataset_time_offset_s = float(loaded.get("time_offset_s") or 0.0)

        col_names = list(self._dataset_columns.keys())
        self._data_col_combo.configure(values=col_names, state="readonly")
        self._data_col_var.set(col_names[0] if col_names else "")

        if time_source == "column":
            self._rate_entry.configure(state="normal")
            self._rate_var.set(f"{rate_hz:g}" if rate_hz is not None else "")
            self._rate_entry.configure(state="disabled")
        else:
            self._rate_entry.configure(state="normal")
            self._rate_var.set(f"{rate_hz:g}" if rate_hz is not None else "")

        self._update_dataset_status()

    def _selected_dataset_trace(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (values, time_s) for the currently selected data column.

        Applies a joint finite mask over time and values so rows stay aligned.
        """
        if not self._dataset_columns:
            raise ValueError("No dataset loaded.")
        col = self._data_col_var.get().strip()
        if not col or col not in self._dataset_columns:
            raise ValueError("Choose a data column from the dropdown.")
        values = np.asarray(self._dataset_columns[col], dtype=float).ravel()

        if self._dataset_time_source == "column" and self._dataset_time_s is not None:
            time_s = np.asarray(self._dataset_time_s, dtype=float).ravel()
            if time_s.shape != values.shape:
                raise ValueError("Time and data columns have different lengths.")
        else:
            rate_txt = self._rate_var.get().strip()
            try:
                rate_hz = float(rate_txt) if rate_txt else float(self._dataset_rate_hz or 0)
            except ValueError as exc:
                raise ValueError("Sampling rate must be a positive number.") from exc
            if not np.isfinite(rate_hz) or rate_hz <= 0:
                raise ValueError("Sampling rate must be a positive number.")
            self._dataset_rate_hz = rate_hz
            time_s = np.arange(values.size, dtype=float) / rate_hz

        mask = np.isfinite(values) & np.isfinite(time_s)
        values = values[mask]
        time_s = time_s[mask]
        if values.size == 0:
            raise ValueError("No finite samples in the selected column.")
        return values, time_s

    def _close_compare_figure(self) -> None:
        for fig in (self._compare_fig,):
            if fig is not None:
                try:
                    plt.close(fig)
                except Exception:
                    pass
        self._compare_fig = None

    def _on_plot(self) -> None:
        model = self._model_for_selected_run()
        if model is None:
            messagebox.showinfo(
                "Data Comparator",
                "Select a model run in the list (add a run first).",
                parent=self,
            )
            return
        series = self._series_var.get().strip()
        if not series or series not in model:
            messagebox.showinfo(
                "Data Comparator",
                "Choose a model series from the dropdown.",
                parent=self,
            )
            return
        if not self._dataset_columns:
            messagebox.showinfo(
                "Data Comparator",
                "Load a dataset file first.",
                parent=self,
            )
            return

        try:
            raw, data_timebase = self._selected_dataset_trace()
        except ValueError as exc:
            messagebox.showinfo("Data Comparator", str(exc), parent=self)
            return

        xaxis = np.asarray(model["xaxis"], dtype=float)
        pred = np.asarray(model[series], dtype=float)
        if pred.shape != xaxis.shape:
            messagebox.showerror(
                "Data Comparator",
                f"Series {series!r} length does not match xaxis.",
                parent=self,
            )
            return

        try:
            synced = syncTimeseries(raw, data_timebase, pred, xaxis)
        except Exception as exc:
            messagebox.showerror(
                "Data Comparator",
                f"Synchronization failed:\n{exc}",
                parent=self,
            )
            return

        synced = np.asarray(synced, dtype=float)
        datatb = np.asarray(xaxis, dtype=float)

        idx = self._selected_run_index()
        run_label = self._runs[idx]["label"] if idx is not None else "Model"
        title = f"{run_label} vs data ({series})"

        self._close_compare_figure()
        if self._fig_placeholder is not None:
            try:
                plt.close(self._fig_placeholder)
            except Exception:
                pass
            self._fig_placeholder = None
        for w in self._plot_host.winfo_children():
            w.destroy()

        try:
            fig = plotModelVsData(
                pred,
                synced,
                xaxis,
                model,
                title=title,
                datatimebase=datatb,
                figsize=(10, 5),
            )
        except Exception as exc:
            messagebox.showerror(
                "Data Comparator",
                f"Plot failed:\n{exc}",
                parent=self,
            )
            self._render_placeholder("Plot failed. Adjust inputs and try again.")
            return

        self._compare_fig = fig
        self._canvas = FigureCanvasTkAgg(fig, master=self._plot_host)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _on_close(self) -> None:
        self._close_compare_figure()
        if self._fig_placeholder is not None:
            try:
                plt.close(self._fig_placeholder)
            except Exception:
                pass
            self._fig_placeholder = None
        self._parent._data_comparator_dialog = None
        self.destroy()


_SPECTRUM_LIBRARY_ENTRIES = (
    "440 nm",
    "560 nm",
    "Xenon",
    "Xenon (eye)",
    "Dark",
)

_DARK_WAVELENGTHS = np.linspace(200.0, 800.0, int(round((800.0 - 200.0) / 0.1)) + 1)
_DARK_INTENSITY = np.zeros_like(_DARK_WAVELENGTHS, dtype=float)


class SpectrumBuilderDialog(tk.Toplevel):
    """Build monochromatic impulse spectra or linear mixtures; save to the custom library."""

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Spectrum Builder")
        self._parent = parent
        self.transient(parent)
        self._mono_peak = tk.StringVar(master=self, value="550.0")
        self._mono_resolved = tk.StringVar(master=self, value="—")
        self._mix_norm_integral = tk.BooleanVar(master=self, value=False)
        self._mix_norm_peak = tk.BooleanVar(master=self, value=True)
        self._mix_rows: list[dict] = []
        self._notebook: ttk.Notebook | None = None
        self._mono_fig = None
        self._mono_ax = None
        self._mono_canvas = None
        self._mix_fig = None
        self._mix_ax = None
        self._mix_canvas = None

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        if parent._env is None:
            ttk.Label(
                outer,
                text=(
                    "Spectra are not loaded. Run the app from the repository root with "
                    "data/ present, wait until the status shows Ready, then reopen Spectrum Builder."
                ),
                wraplength=520,
            ).pack(anchor=tk.W, pady=(0, 8))
            self.minsize(440, 140)
            self.protocol("WM_DELETE_WINDOW", self._on_close)
            return

        nb = ttk.Notebook(outer)
        nb.pack(fill=tk.BOTH, expand=True)
        self._notebook = nb

        mono = ttk.Frame(nb, padding=6)
        nb.add(mono, text="Monochromatic Spectrum")
        row_peak = ttk.Frame(mono)
        row_peak.pack(fill=tk.X)
        ttk.Label(row_peak, text="Target wavelength (nm):").pack(side=tk.LEFT)
        ent_peak = ttk.Entry(row_peak, textvariable=self._mono_peak, width=14)
        ent_peak.pack(side=tk.LEFT, padx=(8, 0))
        ttk.Label(mono, textvariable=self._mono_resolved, wraplength=520).pack(
            anchor=tk.W, pady=(6, 0)
        )
        self._mono_peak.trace_add("write", lambda *_: self._try_mono_preview())
        mono_plot = ttk.Frame(mono)
        mono_plot.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        self._mono_fig = plt.Figure(figsize=(5.2, 3.2), tight_layout=True)
        self._mono_ax = self._mono_fig.add_subplot(111)
        self._mono_canvas = FigureCanvasTkAgg(self._mono_fig, master=mono_plot)
        self._mono_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        ttk.Button(mono, text="Save to library…", command=self._on_save_mono).pack(
            anchor=tk.E, pady=(8, 0)
        )

        mix = ttk.Frame(nb, padding=6)
        nb.add(mix, text="Mixed Spectrum")
        mix_top = ttk.Frame(mix)
        mix_top.pack(fill=tk.X)
        ttk.Checkbutton(
            mix_top,
            text="Normalize rows to integrated photon count",
            variable=self._mix_norm_integral,
            command=lambda: self._on_mix_norm_toggle(_MIX_NORM_INTEGRAL),
        ).pack(side=tk.LEFT)
        ttk.Checkbutton(
            mix_top,
            text="Normalize rows to peak",
            variable=self._mix_norm_peak,
            command=lambda: self._on_mix_norm_toggle(_MIX_NORM_PEAK),
        ).pack(side=tk.LEFT, padx=(12, 0))
        mix_btns = ttk.Frame(mix)
        mix_btns.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(mix_btns, text="Add row", command=self._add_mix_row).pack(
            side=tk.LEFT
        )
        ttk.Button(mix_btns, text="Remove last row", command=self._remove_mix_row).pack(
            side=tk.LEFT, padx=(4, 0)
        )
        self._mix_rows_host = ttk.Frame(mix)
        self._mix_rows_host.pack(fill=tk.X, pady=(8, 0))
        mix_plot = ttk.Frame(mix)
        mix_plot.pack(fill=tk.BOTH, expand=True, pady=(8, 0))
        self._mix_fig = plt.Figure(figsize=(5.2, 3.2), tight_layout=True)
        self._mix_ax = self._mix_fig.add_subplot(111)
        self._mix_canvas = FigureCanvasTkAgg(self._mix_fig, master=mix_plot)
        self._mix_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        ttk.Button(mix, text="Save to library…", command=self._on_save_mix).pack(
            anchor=tk.E, pady=(8, 0)
        )

        self._add_mix_row()
        self._try_mono_preview()

        nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)

        self.minsize(640, 480)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _on_tab_changed(self, _event=None) -> None:
        if self._notebook is None:
            return
        try:
            tab = self._notebook.index(self._notebook.select())
        except tk.TclError:
            return
        if tab == 0:
            self._try_mono_preview()
        else:
            self._update_mix_preview()

    def _try_mono_preview(self) -> None:
        if self._mono_ax is None or self._parent._env is None:
            return
        try:
            peak = float(self._mono_peak.get().strip())
        except ValueError:
            self._mono_resolved.set("—")
            self._mono_ax.clear()
            self._mono_ax.set_title("Monochromatic impulse")
            self._mono_canvas.draw_idle()
            return
        wlen = np.asarray(self._parent._env["wlenshared"], dtype=float)
        try:
            wg, intens = _impulse_monochromatic(wlen, peak)
        except ValueError:
            self._mono_resolved.set("Target is outside the model wavelength grid.")
            self._mono_ax.clear()
            self._mono_ax.set_title("Monochromatic impulse")
            self._mono_canvas.draw_idle()
            return
        idx = int(np.argmax(intens))
        self._mono_resolved.set(
            f"Resolved on model grid: {float(wg[idx]):.4f} nm (sample index {idx})"
        )
        ax = self._mono_ax
        ax.clear()
        ax.plot(wg, intens, color="k")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Intensity (relative)")
        ax.set_title("Monochromatic impulse")
        self._mono_canvas.draw_idle()

    def _on_save_mono(self) -> None:
        env = self._parent._env
        if env is None or self._mono_ax is None:
            messagebox.showwarning("Spectrum Builder", "Spectra are not loaded.", parent=self)
            return
        try:
            peak = float(self._mono_peak.get().strip())
        except ValueError:
            messagebox.showerror("Monochromatic", "Enter a numeric wavelength (nm).", parent=self)
            return
        wlen = np.asarray(env["wlenshared"], dtype=float)
        try:
            wg, intens = _impulse_monochromatic(wlen, peak)
        except ValueError as exc:
            messagebox.showerror("Monochromatic", str(exc), parent=self)
            return
        name = simpledialog.askstring(
            "Save spectrum",
            "Name for this spectrum in the library:",
            parent=self,
        )
        if not name or not name.strip():
            return
        self._parent._register_custom_spectrum(name.strip(), wg, intens, parent=self)

    def _refresh_mix_row_combos(self) -> None:
        vals = _spectrum_picker_entries(self._parent)
        for r in self._mix_rows:
            cb: ttk.Combobox = r["combo"]
            cb["values"] = vals
            cur = cb.get().strip()
            if cur not in vals and vals:
                cb.set(vals[0])

    def _add_mix_row(self) -> None:
        self._refresh_mix_row_combos()
        vals = _spectrum_picker_entries(self._parent)
        rowf = ttk.Frame(self._mix_rows_host)
        rowf.pack(fill=tk.X, pady=2)
        cb = ttk.Combobox(rowf, width=24, state="readonly")
        cb["values"] = vals
        if vals:
            cb.set(vals[0])
        cb.pack(side=tk.LEFT)
        cb.bind("<<ComboboxSelected>>", lambda _e: self._update_mix_preview())
        ttk.Label(rowf, text="Weight:").pack(side=tk.LEFT, padx=(8, 2))
        wvar = tk.StringVar(master=self, value="1.0")
        ttk.Entry(rowf, textvariable=wvar, width=10).pack(side=tk.LEFT)
        wvar.trace_add("write", lambda *_: self._update_mix_preview())
        self._mix_rows.append({"frame": rowf, "combo": cb, "weight": wvar})
        self._refresh_mix_row_combos()
        self._update_mix_preview()

    def _remove_mix_row(self) -> None:
        if len(self._mix_rows) <= 1:
            messagebox.showinfo(
                "Mixture", "At least one spectrum row is required.", parent=self
            )
            return
        last = self._mix_rows.pop()
        last["frame"].destroy()
        self._update_mix_preview()

    def _collect_mix_rows(self) -> list[tuple[str, float]]:
        out: list[tuple[str, float]] = []
        for r in self._mix_rows:
            name = str(r["combo"].get()).strip()
            if not name:
                continue
            try:
                w = float(str(r["weight"].get()).strip())
            except ValueError as exc:
                raise ValueError(
                    f"Invalid weight for spectrum {name!r}."
                ) from exc
            out.append((name, w))
        if not out:
            raise ValueError("Add at least one spectrum row with a name and weight.")
        return out

    def _mix_norm_mode(self) -> str:
        if self._mix_norm_peak.get():
            return _MIX_NORM_PEAK
        if self._mix_norm_integral.get():
            return _MIX_NORM_INTEGRAL
        return _MIX_NORM_NONE

    def _on_mix_norm_toggle(self, mode: str) -> None:
        if mode == _MIX_NORM_PEAK and self._mix_norm_peak.get():
            self._mix_norm_integral.set(False)
        elif mode == _MIX_NORM_INTEGRAL and self._mix_norm_integral.get():
            self._mix_norm_peak.set(False)
        self._update_mix_preview()

    def _update_mix_preview(self) -> None:
        """Redraw the mixture preview; ignore transient invalid edits silently."""
        if self._mix_ax is None:
            return
        try:
            rows = self._collect_mix_rows()
            mode = self._mix_norm_mode()
            wlen, comb = _mix_weighted_spectra(self._parent, rows, mode)
        except ValueError:
            return
        ax = self._mix_ax
        ax.clear()
        ax.plot(wlen, comb, color="k")
        ax.set_xlabel("Wavelength (nm)")
        if mode == _MIX_NORM_INTEGRAL:
            ax.set_ylabel("Normalized intensity (1/nm)")
        elif mode == _MIX_NORM_PEAK:
            ax.set_ylabel("Normalized intensity (peak = 1)")
        else:
            ax.set_ylabel("Intensity (photons/\u00b5m\u00b2/nm/s)")
            ax.ticklabel_format(
                axis="y", style="scientific", scilimits=(0, 0), useMathText=True
            )
        ax.set_title("Weighted mixture")
        self._mix_canvas.draw_idle()

    def _on_save_mix(self) -> None:
        if self._mix_ax is None:
            return
        try:
            rows = self._collect_mix_rows()
            wlen, comb = _mix_weighted_spectra(
                self._parent, rows, self._mix_norm_mode()
            )
        except ValueError as exc:
            messagebox.showerror("Mixture", str(exc), parent=self)
            return
        if not np.any(comb > 0):
            messagebox.showwarning(
                "Mixture",
                "The combined spectrum is all zeros. Adjust weights or spectra.",
                parent=self,
            )
            return
        name = simpledialog.askstring(
            "Save spectrum",
            "Name for this mixture in the library:",
            parent=self,
        )
        if not name or not name.strip():
            return
        self._parent._register_custom_spectrum(name.strip(), wlen, comb, parent=self)
        self._refresh_mix_row_combos()
        self._update_mix_preview()

    def _on_close(self) -> None:
        for fig in (self._mono_fig, self._mix_fig):
            if fig is not None:
                try:
                    plt.close(fig)
                except Exception:
                    pass
        self._parent._spectrum_builder_dialog = None
        self.destroy()


class SpectrumDeleteConflictDialog(tk.Toplevel):
    """Warn that deleting a spectrum will break stimulus protocols; offer resolution."""

    def __init__(
        self,
        parent: tk.Misc,
        app: "ManuscriptSimApp",
        spectrum_name: str,
        affected: list[str],
    ) -> None:
        super().__init__(parent)
        self.title("Delete Spectrum")
        self.transient(parent)
        self.resizable(False, False)
        self.result: tuple[str, str | None] | None = None
        self._spectrum_name = spectrum_name
        self._app = app

        replacements: list[str] = []
        for name in _SPECTRUM_LIBRARY_ENTRIES:
            if name in app._hidden_builtin_spectra:
                continue
            if name == spectrum_name:
                continue
            replacements.append(name)
        for name in app._custom_spectra:
            if name == spectrum_name:
                continue
            if name not in replacements:
                replacements.append(name)

        outer = ttk.Frame(self, padding=10)
        outer.pack(fill=tk.BOTH, expand=True)

        ttk.Label(
            outer,
            text=(
                f"Spectrum {spectrum_name!r} is used by the stimulus protocol(s) "
                "listed below. Deleting it without resolving these references will "
                "break those protocols."
            ),
            wraplength=420,
            justify=tk.LEFT,
        ).pack(anchor=tk.W)

        list_frame = ttk.Frame(outer, padding=(0, 8, 0, 0))
        list_frame.pack(fill=tk.BOTH, expand=True)
        self._listbox = tk.Listbox(
            list_frame, height=min(8, max(3, len(affected))), exportselection=False, width=42
        )
        self._listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(list_frame, orient=tk.VERTICAL, command=self._listbox.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self._listbox.configure(yscrollcommand=sb.set)
        for name in affected:
            self._listbox.insert(tk.END, name)

        self._action_var = tk.StringVar(master=self, value="delete_stimuli")
        opts = ttk.Frame(outer, padding=(0, 8, 0, 0))
        opts.pack(fill=tk.X)

        ttk.Radiobutton(
            opts,
            text="Delete these stimulus protocols",
            variable=self._action_var,
            value="delete_stimuli",
            command=self._on_action_changed,
        ).pack(anchor=tk.W)

        replace_row = ttk.Frame(opts)
        replace_row.pack(anchor=tk.W, fill=tk.X, pady=(4, 0))
        ttk.Radiobutton(
            replace_row,
            text="Replace it in these protocols with:",
            variable=self._action_var,
            value="replace",
            command=self._on_action_changed,
        ).pack(side=tk.LEFT)
        self._replace_var = tk.StringVar(master=self, value=replacements[0] if replacements else "")
        self._replace_combo = ttk.Combobox(
            replace_row,
            textvariable=self._replace_var,
            values=replacements,
            state="disabled",
            width=22,
        )
        self._replace_combo.pack(side=tk.LEFT, padx=(8, 0))
        self._replacements = replacements

        btn_row = ttk.Frame(outer, padding=(0, 12, 0, 0))
        btn_row.pack(fill=tk.X)
        self._confirm_btn = ttk.Button(
            btn_row, text="Delete Spectrum", command=self._on_confirm
        )
        self._confirm_btn.pack(side=tk.LEFT)
        ttk.Button(btn_row, text="Cancel", command=self._on_cancel).pack(side=tk.RIGHT)

        self.protocol("WM_DELETE_WINDOW", self._on_cancel)
        self.grab_set()
        self.wait_window(self)

    def _on_action_changed(self) -> None:
        if self._action_var.get() == "replace":
            self._confirm_btn.configure(text="Replace and Delete")
            if self._replacements:
                self._replace_combo.configure(state="readonly")
            else:
                self._replace_combo.configure(state="disabled")
        else:
            self._confirm_btn.configure(text="Delete Spectrum")
            self._replace_combo.configure(state="disabled")

    def _on_confirm(self) -> None:
        action = self._action_var.get()
        if action == "delete_stimuli":
            self.result = ("delete_stimuli", None)
            self.destroy()
            return
        if action == "replace":
            replacement = str(self._replace_var.get()).strip()
            if not replacement or replacement not in self._replacements:
                messagebox.showinfo(
                    "Delete Spectrum",
                    "Select a replacement spectrum.",
                    parent=self,
                )
                return
            self.result = ("replace", replacement)
            self.destroy()
            return
        messagebox.showinfo(
            "Delete Spectrum",
            "Choose how to resolve the affected stimulus protocols.",
            parent=self,
        )

    def _on_cancel(self) -> None:
        self.result = None
        self.destroy()


class SpectrumLibraryPopup(tk.Toplevel):
    """Popup showing individual loaded spectra with a single-trace preview."""

    def __init__(self, parent: "ManuscriptSimApp") -> None:
        super().__init__(parent)
        self.title("Spectrum Library")
        self._parent = parent
        self.transient(parent)

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        bottom_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        bottom_row.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Button(bottom_row, text="Close", command=self._on_close).pack(side=tk.RIGHT)

        library_body = ttk.Frame(outer, padding=8)
        library_body.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self._build_library_tab(library_body)

        self.minsize(760, 420)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _build_library_tab(self, container) -> None:
        btn_row = ttk.Frame(container, padding=(0, 8, 0, 0))
        btn_row.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Button(btn_row, text="Import Spectrum", command=self._on_add_spectrum).pack(side=tk.LEFT)
        self._delete_btn = ttk.Button(
            btn_row, text="Delete Spectrum", command=self._on_delete_spectrum,
            state=tk.DISABLED,
        )
        self._delete_btn.pack(side=tk.LEFT, padx=(8, 0))
        self._crop_btn = ttk.Button(
            btn_row, text="Crop Spectrum", command=self._on_crop_spectrum,
            state=tk.DISABLED,
        )
        self._crop_btn.pack(side=tk.LEFT, padx=(8, 0))

        panes = ttk.Frame(container)
        panes.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        left_col = ttk.Frame(panes)
        left_col.pack(side=tk.LEFT, fill=tk.BOTH)
        self._all_spectrum_names: list[str] = []

        list_frame = ttk.Frame(left_col)
        list_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        scrollbar = ttk.Scrollbar(list_frame, orient=tk.VERTICAL)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        list_inner = ttk.Frame(list_frame)
        list_inner.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._spectrum_search_var = _install_listbox_search_row(list_inner)
        self._spectrum_search_var.trace_add(
            "write", lambda *_: self._apply_spectrum_search_filter()
        )
        self._listbox = tk.Listbox(list_inner, height=12, exportselection=False, width=28)
        self._listbox.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self._listbox.configure(yscrollcommand=scrollbar.set)
        scrollbar.configure(command=self._listbox.yview)

        plot_frame = ttk.Frame(panes, padding=(8, 0, 0, 0))
        plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._fig = plt.Figure(figsize=(5, 3.2), tight_layout=True)
        self._ax = self._fig.add_subplot(111)
        self._canvas = FigureCanvasTkAgg(self._fig, master=plot_frame)
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self._listbox.bind("<<ListboxSelect>>", self._on_selection_changed)
        self._rebuild_spectrum_name_list()
        self._apply_spectrum_search_filter()


    def _rebuild_spectrum_name_list(self) -> None:
        parent = self._parent
        names: list[str] = []
        for name in _SPECTRUM_LIBRARY_ENTRIES:
            if name in parent._hidden_builtin_spectra:
                continue
            names.append(name)
        for name in parent._custom_spectra:
            if name not in names:
                names.append(name)
        self._all_spectrum_names = names

    def _apply_spectrum_search_filter(self, *, select_name: str | None = None) -> None:
        previous = (
            select_name
            if select_name is not None
            else self._selected_library_spectrum()
        )
        filtered = _names_matching_search(
            self._all_spectrum_names, self._spectrum_search_var.get()
        )
        if previous is not None and previous not in filtered:
            if previous not in self._all_spectrum_names:
                previous = None
        selected = _listbox_fill_filtered(
            self._listbox, filtered, select_name=previous
        )
        self._render_spectrum(selected)
        self._update_crop_button_state()

    def _spectrum_for_name(
        self, name: str
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Return (wlen, intensity) for a single named spectrum, or None if env not loaded."""
        return _library_spectrum_arrays(self._parent, name)

    def append_entry_if_missing(self, name: str) -> None:
        """Insert ``name`` into the library list if it is not already listed."""
        if name not in self._all_spectrum_names:
            self._all_spectrum_names.append(name)
        self._apply_spectrum_search_filter(select_name=name)

    def _render_spectrum(self, name: str | None) -> None:
        ax = self._ax
        ax.clear()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Intensity (photons/\u00b5m\u00b2/nm/s)")
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0), useMathText=True)
        if name is None:
            ax.set_title("Select a spectrum")
            ax.set_xlim(100, 1000)
        else:
            result = self._spectrum_for_name(name)
            if result is None:
                ax.set_title(f"Spectrum: {name}")
                ax.set_xlim(100, 1000)
                ax.text(
                    0.5, 0.5,
                    "Spectra not loaded yet.",
                    ha="center", va="center",
                    transform=ax.transAxes,
                )
            else:
                wlen, intensity = result
                ax.plot(wlen, intensity, color="k")
                ax.set_title(f"Spectrum: {name}")
                wmin = float(np.min(wlen))
                wmax = float(np.max(wlen))
                if wmin < wmax:
                    ax.set_xlim(wmin, wmax)
                else:
                    ax.set_xlim(100, 1000)
        self._canvas.draw_idle()

    def _selected_library_spectrum(self) -> str | None:
        sel = self._listbox.curselection()
        if not sel:
            return None
        return str(self._listbox.get(sel[0]))

    @staticmethod
    def _is_manuscript_spectrum(name: str) -> bool:
        return name in _SPECTRUM_LIBRARY_ENTRIES

    def _update_crop_button_state(self) -> None:
        name = self._selected_library_spectrum()
        can_edit = name is not None and not self._is_manuscript_spectrum(name)
        state = tk.NORMAL if can_edit else tk.DISABLED
        self._crop_btn.config(state=state)
        self._delete_btn.config(state=state)

    def _on_selection_changed(self, _event=None) -> None:
        sel = self._listbox.curselection()
        if not sel:
            self._render_spectrum(None)
            self._update_crop_button_state()
            return
        self._render_spectrum(str(self._listbox.get(sel[0])))
        self._update_crop_button_state()

    def _on_add_spectrum(self) -> None:
        name = self._parent._add_spectrum_via_dialog(parent=self)
        if name is None:
            return
        if name not in self._all_spectrum_names:
            self._all_spectrum_names.append(name)
        self._apply_spectrum_search_filter(select_name=name)

    def _on_delete_spectrum(self) -> None:
        sel = self._listbox.curselection()
        if not sel:
            messagebox.showinfo("Delete Spectrum", "Select a spectrum to delete.", parent=self)
            return
        idx = sel[0]
        name = str(self._listbox.get(idx))
        if self._is_manuscript_spectrum(name):
            messagebox.showinfo(
                "Delete Spectrum",
                "Manuscript spectra cannot be deleted or modified.",
                parent=self,
            )
            return
        affected = self._parent._stimuli_using_spectrum(name)
        if affected:
            dlg = SpectrumDeleteConflictDialog(self, self._parent, name, affected)
            if dlg.result is None:
                return
            if not self._parent._apply_spectrum_delete_resolution(name, affected, dlg.result):
                return
        else:
            ok = messagebox.askyesno(
                "Delete Spectrum",
                f"Are you sure you want to delete spectrum {name!r}?",
                parent=self,
            )
            if not ok:
                return
        if not self._parent._delete_spectrum(name):
            return
        if name in self._all_spectrum_names:
            self._all_spectrum_names.remove(name)
        filtered = _names_matching_search(
            self._all_spectrum_names, self._spectrum_search_var.get()
        )
        next_name = filtered[min(idx, len(filtered) - 1)] if filtered else None
        self._apply_spectrum_search_filter(select_name=next_name)

    def _on_crop_spectrum(self) -> None:
        sel = self._listbox.curselection()
        if not sel:
            return
        name = str(self._listbox.get(sel[0]))
        if self._is_manuscript_spectrum(name):
            messagebox.showinfo(
                "Crop Spectrum",
                "Manuscript spectra cannot be cropped or modified.",
                parent=self,
            )
            return
        result = self._spectrum_for_name(name)
        if result is None:
            messagebox.showinfo(
                "Crop Spectrum", "Spectra not loaded yet.", parent=self
            )
            return
        existing = getattr(self, "_crop_dialog", None)
        if existing is not None:
            try:
                if existing.winfo_exists():
                    existing.lift()
                    return
            except tk.TclError:
                pass
        wlen, intensity = result
        self._crop_dialog = CropSpectrumDialog(
            self, self._parent, name,
            np.asarray(wlen, dtype=float), np.asarray(intensity, dtype=float),
        )


    def _on_close(self) -> None:
        try:
            plt.close(self._fig)
        except Exception:
            pass
        self._parent._spectrum_library_popup = None
        self.destroy()


class CropSpectrumDialog(tk.Toplevel):
    """Interactive wavelength-axis crop for a single spectrum.

    The user drags the edges of a wavelength window (synced to numeric min/max
    entries). Intensities inside the window are kept; everything outside is set
    to zero. The result is saved as a new custom spectrum; the original arrays
    are never modified.
    """

    def __init__(
        self,
        popup: "SpectrumLibraryPopup",
        app: "ManuscriptSimApp",
        name: str,
        wlen: np.ndarray,
        intensity: np.ndarray,
        *,
        title: str | None = None,
        instruction: str | None = None,
        save_text: str = "Save as new spectrum…",
        on_apply=None,
    ) -> None:
        super().__init__(popup)
        self._popup = popup
        self._app = app
        self._name = name
        # Optional hook: callable(lo, hi) -> bool | None. When provided, it is
        # invoked instead of the default single-spectrum save (used by the
        # long-form movie crop). Returning False keeps the dialog open.
        self._on_apply = on_apply
        # Defensive copies so the source spectrum is never mutated.
        self._wlen = np.asarray(wlen, dtype=float).copy()
        self._intensity = np.asarray(intensity, dtype=float).copy()
        self._wmin = float(self._wlen.min())
        self._wmax = float(self._wlen.max())
        self._lo = self._wmin
        self._hi = self._wmax

        self.title(title or f"Crop Spectrum — {name}")
        self.transient(popup)
        self.minsize(560, 420)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        outer = ttk.Frame(self, padding=8)
        outer.pack(fill=tk.BOTH, expand=True)

        ttk.Label(
            outer,
            text=(
                instruction
                or (
                    "Drag the edges of the highlighted window (or edit min/max). "
                    "Intensity outside the window is set to 0."
                )
            ),
            wraplength=520,
            justify=tk.LEFT,
        ).pack(side=tk.TOP, anchor=tk.W, pady=(0, 6))

        plot_frame = ttk.Frame(outer)
        plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self._fig = plt.Figure(figsize=(5.5, 3.4), tight_layout=True)
        self._ax = self._fig.add_subplot(111)
        self._ax.set_xlabel("Wavelength (nm)")
        self._ax.set_ylabel("Intensity (photons/\u00b5m\u00b2/nm/s)")
        self._ax.ticklabel_format(
            axis="y", style="scientific", scilimits=(0, 0), useMathText=True
        )
        self._ax.plot(
            self._wlen, self._intensity, color="0.75", linewidth=1.0,
            label="Original",
        )
        (self._preview_line,) = self._ax.plot(
            self._wlen, self._intensity, color="C0", linewidth=1.4,
            label="Cropped",
        )
        self._ax.legend(fontsize=8, loc="upper right")
        self._canvas = FigureCanvasTkAgg(self._fig, master=plot_frame)
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self._span = SpanSelector(
            self._ax,
            self._on_span,
            "horizontal",
            interactive=True,
            drag_from_anywhere=True,
            useblit=False,
            props=dict(alpha=0.15, facecolor="C0"),
        )
        self._span.extents = (self._lo, self._hi)

        entry_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        entry_row.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(entry_row, text="Min (nm):").pack(side=tk.LEFT)
        self._min_var = tk.StringVar(master=self, value=f"{self._lo:g}")
        min_entry = ttk.Entry(entry_row, textvariable=self._min_var, width=10)
        min_entry.pack(side=tk.LEFT, padx=(4, 12))
        min_entry.bind("<Return>", self._on_entry_commit)
        min_entry.bind("<FocusOut>", self._on_entry_commit)
        ttk.Label(entry_row, text="Max (nm):").pack(side=tk.LEFT)
        self._max_var = tk.StringVar(master=self, value=f"{self._hi:g}")
        max_entry = ttk.Entry(entry_row, textvariable=self._max_var, width=10)
        max_entry.pack(side=tk.LEFT, padx=(4, 0))
        max_entry.bind("<Return>", self._on_entry_commit)
        max_entry.bind("<FocusOut>", self._on_entry_commit)
        ttk.Button(
            entry_row, text="Reset", command=self._on_reset
        ).pack(side=tk.RIGHT)

        btn_row = ttk.Frame(outer, padding=(0, 8, 0, 0))
        btn_row.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Button(
            btn_row, text=save_text, command=self._on_save
        ).pack(side=tk.LEFT)
        ttk.Button(
            btn_row, text="Cancel", command=self._on_close
        ).pack(side=tk.RIGHT)

        self._refresh_preview()

    def _set_window(self, lo: float, hi: float, from_span: bool = False) -> None:
        """Clamp/normalize the window and sync span, entries, and preview."""
        lo = float(lo)
        hi = float(hi)
        if hi < lo:
            lo, hi = hi, lo
        lo = min(max(lo, self._wmin), self._wmax)
        hi = min(max(hi, self._wmin), self._wmax)
        self._lo = lo
        self._hi = hi
        self._min_var.set(f"{lo:g}")
        self._max_var.set(f"{hi:g}")
        if not from_span:
            try:
                self._span.extents = (lo, hi)
            except Exception:
                pass
        self._refresh_preview()

    def _refresh_preview(self) -> None:
        mask = (self._wlen >= self._lo) & (self._wlen <= self._hi)
        cropped = np.where(mask, self._intensity, 0.0)
        self._preview_line.set_ydata(cropped)
        self._canvas.draw_idle()

    def _on_span(self, lo: float, hi: float) -> None:
        self._set_window(lo, hi, from_span=True)

    def _on_entry_commit(self, _event=None) -> None:
        try:
            lo = float(self._min_var.get())
            hi = float(self._max_var.get())
        except ValueError:
            # Revert to current values on invalid input.
            self._min_var.set(f"{self._lo:g}")
            self._max_var.set(f"{self._hi:g}")
            return
        self._set_window(lo, hi)

    def _on_reset(self) -> None:
        self._set_window(self._wmin, self._wmax)

    def _on_save(self) -> None:
        if self._on_apply is not None:
            if self._on_apply(self._lo, self._hi) is False:
                return
            self._on_close()
            return
        mask = (self._wlen >= self._lo) & (self._wlen <= self._hi)
        cropped = np.where(mask, self._intensity, 0.0)
        default = f"{self._name} cropped"
        new_name = simpledialog.askstring(
            "Crop Spectrum",
            "Name for the cropped spectrum:",
            initialvalue=default,
            parent=self,
        )
        if not new_name or not new_name.strip():
            return
        registered = self._app._register_custom_spectrum(
            new_name.strip(), self._wlen, cropped, parent=self
        )
        if registered is None:
            return
        lb = self._popup._listbox
        self._popup.append_entry_if_missing(registered)
        for i in range(lb.size()):
            if str(lb.get(i)) == registered:
                lb.selection_clear(0, tk.END)
                lb.selection_set(i)
                lb.see(i)
                break
        self._popup._render_spectrum(registered)
        self._popup._update_crop_button_state()
        self._on_close()

    def _on_close(self) -> None:
        try:
            self._span.disconnect_events()
        except Exception:
            pass
        try:
            plt.close(self._fig)
        except Exception:
            pass
        if getattr(self._popup, "_crop_dialog", None) is self:
            self._popup._crop_dialog = None
        self.destroy()





class ManuscriptSimApp(tk.Tk):
    def __init__(self, splash: StartupSplash | None = None) -> None:
        super().__init__()
        # Keep the main window hidden until startup work finishes.
        self.withdraw()

        if splash is not None:
            self._startup_splash = splash
        else:
            self._startup_splash = StartupSplash(version=__version__)

        self.title(f"Melanopsin Model v{__version__}")
        self._env: dict | None = None
        self._plot_frame: ttk.Frame | None = None
        self._canvas: FigureCanvasTkAgg | None = None
        self._last_fig = None
        self._last_result: dict | None = None
        self._last_label: str | None = None
        self._model_config: dict = dict(get_predict_melanopsin_defaults())
        self._stimulus_duration_override: float | None = None
        self._config_dialog: ModelConfigDialog | None = None
        self._custom_stim_dialog: CustomStimulusDialog | None = None
        self._stimulus_creator_dialog: StimulusCreatorDialog | None = None
        self._stimulus_library_popup: StimulusLibraryPopup | None = None
        self._spectrum_library_popup: SpectrumLibraryPopup | None = None
        self._spectrum_builder_dialog: SpectrumBuilderDialog | None = None
        self._prediction_compare_dialog: PredictionComparisonDialog | None = None
        self._data_comparator_dialog: DataComparatorDialog | None = None

        self._startup_splash.update("Preparing user library…", 15)
        try:
            library_storage.reset_legacy_paths_if_needed()
        except Exception:
            # The wipe is best-effort; failures should not prevent startup.
            library_storage.ensure_library_tree()

        self._startup_splash.update("Loading custom spectra…", 25)
        self._custom_spectra: dict[str, tuple[np.ndarray, np.ndarray]] = (
            library_storage.load_all_spectra()
        )
        self._hidden_builtin_spectra: set[str] = set()

        self._startup_splash.update("Loading stimulus library…", 40)
        self._custom_stimuli: dict[str, dict] = library_storage.load_all_stimuli()
        self._known_stimuli: list[str] = list(LABELS) + sorted(
            self._custom_stimuli.keys()
        )
        self._selected_stimulus_var = tk.StringVar(
            master=self, value=NO_STIMULUS_SELECTED_LABEL
        )
        self._run_lock = threading.Lock()
        self._cancel_event = threading.Event()
        self._progress_queue: Queue[tuple[int, int]] = Queue()
        self._progress_poll_id: str | None = None
        self._autosave_images = tk.BooleanVar(master=self, value=False)
        self._autosave_data = tk.BooleanVar(master=self, value=False)
        self._autosave_images_dir: Path = Path(DEFAULT_AUTOSAVE_IMAGES_DIR)
        self._autosave_data_dir: Path = Path(DEFAULT_AUTOSAVE_DATA_DIR)
        self._autosave_options_dialog: AutoSaveOptionsDialog | None = None

        self._startup_splash.update("Building interface…", 60)

        menubar = tk.Menu(self, tearoff=0)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(
            label="Export Model Run",
            command=self._on_file_export_data,
            accelerator="Ctrl+S",
        )
        file_menu.add_command(
            label="Import Spectrum",
            command=self._on_file_import_spectrum,
        )
        file_menu.add_command(
            label="Import Stimulus",
            command=self._on_file_import_stimulus,
        )
        file_menu.add_separator()
        file_menu.add_command(
            label="Exit",
            command=self._quit_application,
            accelerator="Ctrl+Q",
        )
        menubar.add_cascade(label="File", menu=file_menu)

        configure_menu = tk.Menu(menubar, tearoff=0)
        configure_menu.add_command(
            label="Model parameters…",
            command=self._on_configure_model,
        )
        menubar.add_cascade(label="Configure", menu=configure_menu)

        stimulus_menu = tk.Menu(menubar, tearoff=0)
        stimulus_menu.add_command(
            label="Stimulus Library",
            command=self._on_stimulus_spectrum_library,
        )
        stimulus_menu.add_command(
            label="Stimulus Builder",
            command=self._on_stimulus_creator,
        )
        stimulus_menu.add_command(
            label="Spectrum Library",
            command=self._on_spectrum_library,
        )
        stimulus_menu.add_command(
            label="Spectrum Builder",
            command=self._on_spectrum_builder,
        )
        menubar.add_cascade(label="Data", menu=stimulus_menu)

        plot_menu = tk.Menu(menubar, tearoff=0)
        plot_menu.add_command(
            label="Compare saved predictions…",
            command=self._on_prediction_compare,
        )
        plot_menu.add_command(
            label="Data Comparator…",
            command=self._on_data_comparator,
        )
        menubar.add_cascade(label="Plot", menu=plot_menu)

        settings_menu = tk.Menu(menubar, tearoff=0)
        settings_menu.add_command(
            label="Auto-save options",
            command=self._open_autosave_options,
        )
        menubar.add_cascade(label="Settings", menu=settings_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(
            label="Open README…",
            command=self._open_readme_locally,
        )
        help_menu.add_command(
            label="Open GitHub…",
            command=self._on_open_github,
        )
        menubar.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menubar)

        self.bind_all("<Control-s>", self._accel_export_data)
        self.bind_all("<Control-q>", self._accel_exit)
        self.protocol("WM_DELETE_WINDOW", self._quit_application)

        top = ttk.Frame(self, padding=8)
        top.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(top, text="Stimulus:").pack(side=tk.LEFT, padx=(0, 8))
        self._selected_display = ttk.Combobox(
            top,
            textvariable=self._selected_stimulus_var,
            values=tuple(self._known_stimuli),
            state="readonly",
            width=32,
        )
        self._selected_display.pack(side=tk.LEFT)
        self._selected_display.bind(
            "<<ComboboxSelected>>", self._on_stimulus_combobox_selected
        )
        self._selected_display.bind(
            "<Button-1>", self._on_stimulus_combobox_click, add="+"
        )

        self._run_btn = tk.Button(
            top,
            text="\u25b6 Run model",
            command=self._on_run,
            compound=tk.LEFT,
            fg="#1a9f1a",
        )
        self._run_btn.pack(side=tk.LEFT, padx=(12, 0))
        self._stop_btn = tk.Button(
            top,
            text="\u25a0 Stop",
            command=self._on_stop,
            compound=tk.LEFT,
            fg="#cc0000",
            state=tk.DISABLED,
        )
        self._stop_btn.pack(side=tk.LEFT, padx=(6, 0))

        self._progress = ttk.Progressbar(
            top, mode="determinate", length=150, orient=tk.HORIZONTAL, maximum=100
        )
        self._progress.pack(side=tk.LEFT, padx=(12, 0))

        self._status = ttk.Label(top, text="Loading…")
        self._status.pack(side=tk.LEFT, padx=(16, 0))

        self._plot_frame = ttk.Frame(self, padding=(12, 4, 12, 16))
        self._plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        if self._startup_splash is not None:
            self._startup_splash.update("Loading manuscript spectra…", 80)
        self.after(50, self._load_env)

    def _finish_startup(self) -> None:
        """Close the startup splash (if any) and reveal the main window."""
        if self._startup_splash is not None:
            self._startup_splash.close()
            self._startup_splash = None
        # run_melanopsin_gui.py creates the splash's Tk root before this root exists,
        # so destroying the splash leaves tkinter without a default root. Point it at
        # the app so both launch paths behave the same for any implicit-master call.
        if getattr(tk, "_support_default_root", False) and tk._default_root is None:
            tk._default_root = self
        self.deiconify()

    def _load_env(self) -> None:
        def _worker() -> None:
            try:
                env = load()
                err: Exception | None = None
            except Exception as exc:  # noqa: BLE001 - surfaced to the user
                env = None
                err = exc
            self.after(0, lambda: self._on_env_loaded(env, err))

        threading.Thread(target=_worker, daemon=True).start()

    def _on_env_loaded(self, env: dict | None, err: Exception | None) -> None:
        if err is None:
            self._env = env
            self._status.config(text="Ready.")
            if self._startup_splash is not None:
                self._startup_splash.update("Ready.", 100)
        else:
            self._status.config(text="Failed to load spectra.")
            if self._startup_splash is not None:
                self._startup_splash.update("Failed to load spectra.", 100)
            messagebox.showerror(
                "Load error",
                "Could not load spectra (run from repo root with data/ present):\n"
                f"{err}",
            )
        self.after(200, self._finish_startup)

    def _on_run(self) -> None:
        if self._env is None:
            messagebox.showwarning("Not ready", "Spectra are still loading.")
            return
        label = self._selected_stimulus_var.get()
        if label not in self._known_stimuli:
            messagebox.showwarning("No stimulus selected", "Please select a stimulus before running.")
            return
        if not self._run_lock.acquire(blocking=False):
            return
        self._cancel_event.clear()
        self._run_btn.config(state=tk.DISABLED)
        self._stop_btn.config(state=tk.NORMAL)
        self._progress.configure(maximum=100, value=0)
        self._status.config(text="Running model… 0% (0/0)")
        self._start_progress_poll()

        def work() -> None:
            def on_progress(current_step: int, total_steps: int) -> None:
                self._progress_queue.put((current_step, total_steps))

            try:
                stim = self._build_run_protocol(label)
                if (
                    label in LABELS
                    and self._stimulus_duration_override is not None
                ):
                    stim["ti"] = float(self._stimulus_duration_override)
                result = predictMelanopsin(
                    stim,
                    config=self._model_config,
                    progress_cb=on_progress,
                    progress_every=200,
                    cancel_event=self._cancel_event,
                )
                self.after(0, lambda: self._on_done_ok(result, label))
            except SimulationCancelled:
                self.after(0, self._on_done_cancelled)
            except Exception as exc:
                self.after(0, lambda err=exc: self._on_done_err(err))
            finally:
                self.after(0, self._finish_run)

        threading.Thread(target=work, daemon=True).start()

    def _on_stop(self) -> None:
        self._cancel_event.set()

    def _finish_run(self) -> None:
        self._stop_progress_poll()
        self._run_btn.config(state=tk.NORMAL)
        self._stop_btn.config(state=tk.DISABLED)
        self._run_lock.release()

    def _on_done_cancelled(self) -> None:
        self._progress.configure(value=0)
        self._status.config(text="Stopped.")

    def _on_done_err(self, exc: Exception) -> None:
        self._progress.configure(value=0)
        self._status.config(text="Error.")
        messagebox.showerror("Simulation error", str(exc))

    def _on_done_ok(self, result: dict, label: str) -> None:
        self._progress.configure(value=100)
        self._status.config(text="Done.")
        self._last_result = result
        self._last_label = label
        if self._last_fig is not None:
            plt.close(self._last_fig)
            self._last_fig = None

        for w in self._plot_frame.winfo_children():
            w.destroy()

        fig, _axs = plotModelRun(
            result,
            figsize=(13, 9),
            title=str(label),
            show=False,
        )
        self._last_fig = fig
        _ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        if self._autosave_images.get():
            _out_dir = self._autosave_images_dir
            _out_dir.mkdir(parents=True, exist_ok=True)
            try:
                fig.savefig(_out_dir / f"{label}_{_ts}.png")
            except OSError as exc:
                self._status.config(text="Error saving figure to disk.")
                messagebox.showerror("Auto-save", f"Could not save figure:\n{exc}")
        if self._autosave_data.get():
            out_data = self._autosave_data_dir
            out_data.mkdir(parents=True, exist_ok=True)
            data_path = out_data / f"{label}_{_ts}.csv"
            try:
                _export_result_csv(str(data_path), result)
            except (OSError, ValueError) as exc:
                self._status.config(text="Error saving data to disk.")
                messagebox.showerror("Auto-save", f"Could not save data:\n{exc}")
        canvas = FigureCanvasTkAgg(fig, master=self._plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=8, pady=8)
        self._canvas = canvas

    def _accel_export_data(self, event=None) -> str | None:
        self._on_file_export_data()
        return "break"

    def _accel_exit(self, event=None) -> str | None:
        self._quit_application()
        return "break"

    def _quit_application(self) -> None:
        """Exit ``mainloop`` and release resources so the process returns to the shell."""
        self._cancel_event.set()
        try:
            self._stop_progress_poll()
        except tk.TclError:
            pass
        try:
            if self._last_fig is not None:
                plt.close(self._last_fig)
        except Exception:
            pass
        self._last_fig = None
        try:
            plt.close("all")
        except Exception:
            pass
        try:
            self.quit()
        except tk.TclError:
            pass
        try:
            self.destroy()
        except tk.TclError:
            pass

    def _on_configure_model(self) -> None:
        if self._config_dialog is not None:
            try:
                if self._config_dialog.winfo_exists():
                    self._config_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._config_dialog = ModelConfigDialog(self)

    def _open_autosave_options(self) -> None:
        if self._autosave_options_dialog is not None:
            try:
                if self._autosave_options_dialog.winfo_exists():
                    self._autosave_options_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._autosave_options_dialog = AutoSaveOptionsDialog(self)

    def _open_readme_locally(self) -> None:
        """Open the repository README.md in the OS default app."""
        readme = project_root() / "README.md"
        if not readme.is_file():
            messagebox.showinfo(
                "Help",
                f"README not found at:\n{readme}\n\n"
                "Place the app inside the cloned repository (next to README.md).",
                parent=self,
            )
            return
        try:
            if sys.platform == "win32":
                os.startfile(readme)  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.run(["open", str(readme)], check=False)
            else:
                subprocess.run(["xdg-open", str(readme)], check=False)
        except OSError as exc:
            messagebox.showerror(
                "Help",
                f"Could not open README:\n{exc}",
                parent=self,
            )

    def _on_open_github(self) -> None:
        """Open the project GitHub repository in the default browser."""
        webbrowser.open(_GITHUB_REPO_URL)

    def _on_stimulus_spectrum_library(self) -> None:
        self._open_stimulus_library_popup()

    def _on_spectrum_library(self) -> None:
        self._open_spectrum_library_popup()

    def _on_spectrum_builder(self) -> None:
        self._open_spectrum_builder()

    def _on_prediction_compare(self) -> None:
        self._open_prediction_compare_dialog()

    def _on_data_comparator(self) -> None:
        self._open_data_comparator_dialog()

    def _on_stimulus_creator(self) -> None:
        if self._stimulus_creator_dialog is not None:
            try:
                if self._stimulus_creator_dialog.winfo_exists():
                    self._stimulus_creator_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._stimulus_creator_dialog = StimulusBuilderDialog(self)

    def _open_custom_stimulus_builder(self) -> None:
        if self._custom_stim_dialog is not None:
            try:
                if self._custom_stim_dialog.winfo_exists():
                    self._custom_stim_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._custom_stim_dialog = CustomStimulusDialog(self)

    def _open_stimulus_creator(self) -> None:
        if self._stimulus_creator_dialog is not None:
            try:
                if self._stimulus_creator_dialog.winfo_exists():
                    self._stimulus_creator_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._stimulus_creator_dialog = StimulusCreatorDialog(self)

    def _open_stimulus_library_popup(self) -> None:
        if self._stimulus_library_popup is not None:
            try:
                if self._stimulus_library_popup.winfo_exists():
                    self._stimulus_library_popup.lift()
                    return
            except tk.TclError:
                pass
        self._stimulus_library_popup = StimulusLibraryPopup(self)

    def _open_spectrum_library_popup(self) -> None:
        if self._spectrum_library_popup is not None:
            try:
                if self._spectrum_library_popup.winfo_exists():
                    self._spectrum_library_popup.lift()
                    return
            except tk.TclError:
                pass
        self._spectrum_library_popup = SpectrumLibraryPopup(self)

    def _open_spectrum_builder(self) -> None:
        if self._spectrum_builder_dialog is not None:
            try:
                if self._spectrum_builder_dialog.winfo_exists():
                    self._spectrum_builder_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._spectrum_builder_dialog = SpectrumBuilderDialog(self)

    def _open_prediction_compare_dialog(self) -> None:
        if self._prediction_compare_dialog is not None:
            try:
                if self._prediction_compare_dialog.winfo_exists():
                    self._prediction_compare_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._prediction_compare_dialog = PredictionComparisonDialog(self)

    def _open_data_comparator_dialog(self) -> None:
        if self._data_comparator_dialog is not None:
            try:
                if self._data_comparator_dialog.winfo_exists():
                    self._data_comparator_dialog.lift()
                    return
            except tk.TclError:
                pass
        self._data_comparator_dialog = DataComparatorDialog(self)

    def _set_selected_stimulus(self, label: str) -> None:
        self._selected_stimulus_var.set(label)

    def _refresh_stimulus_dropdown_values(self) -> None:
        """Keep the main-window stimulus combobox values in sync with the library."""
        try:
            self._selected_display["values"] = tuple(self._known_stimuli)
        except (tk.TclError, AttributeError):
            pass

    def _on_stimulus_combobox_selected(self, _event=None) -> None:
        label = self._selected_stimulus_var.get()
        if label in self._known_stimuli:
            self._set_selected_stimulus(label)
        try:
            self._selected_display.selection_clear()
        except tk.TclError:
            pass
        self.focus_set()

    def _on_stimulus_combobox_click(self, _event=None) -> None:
        self._refresh_stimulus_dropdown_values()

    def _add_spectrum_via_dialog(self, parent=None) -> str | None:
        path = filedialog.askopenfilename(
            parent=parent if parent is not None else self,
            title="Select spectrum file",
            filetypes=[
                ("Spreadsheets", "*.csv *.xlsx *.xls *.tsv *.txt"),
                ("CSV files", "*.csv"),
                ("Excel files", "*.xlsx *.xls"),
                ("TSV / text files", "*.tsv *.txt"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return None
        try:
            wlen, intensity, note = _read_spectrum_file(path)
        except Exception as exc:
            messagebox.showerror("Add Spectrum", f"Could not read file:\n{exc}", parent=parent)
            return None
        if note:
            messagebox.showwarning("Add Spectrum", note, parent=parent)
        unit_key = _prompt_spectrum_units(parent if parent is not None else self)
        if unit_key is None:
            return None
        try:
            intensity = _convert_spectrum_to_photon_flux(wlen, intensity, unit_key)
        except Exception as exc:
            messagebox.showerror(
                "Add Spectrum", f"Could not convert units:\n{exc}", parent=parent
            )
            return None
        name = simpledialog.askstring(
            "Add Spectrum",
            "Enter a name for this spectrum:",
            parent=parent if parent is not None else self,
        )
        if not name or not name.strip():
            return None
        return self._register_custom_spectrum(
            name.strip(), wlen, intensity, parent=parent
        )

    def _ensure_custom_spectrum_persisted(self, name: str) -> None:
        """Persist an in-memory custom spectrum if it has not yet been written.

        Used when saving a stimulus that references custom spectra by name
        (for example, spectra brought in by loading another stimulus). Built-
        ins and spectra already in the registry are no-ops.
        """
        name = (name or "").strip()
        if not name or name in _BUILTIN_STIMULUS_SPECTRA:
            return
        data = self._custom_spectra.get(name)
        if data is None:
            return
        if library_storage.has_spectrum(name):
            return
        wlen, intensity = data
        library_storage.save_spectrum(name, wlen, intensity)

    def _register_custom_spectrum(
        self,
        name: str,
        wlen: np.ndarray,
        intensity: np.ndarray,
        parent=None,
    ) -> str | None:
        """Persist a spectrum to the user library and register it in memory."""
        name = (name or "").strip()
        if not name:
            messagebox.showerror(
                "Spectrum", "Enter a non-empty name.", parent=parent or self
            )
            return None
        existing = set(_SPECTRUM_LIBRARY_ENTRIES) | set(self._custom_spectra)
        if name in existing:
            messagebox.showerror(
                "Spectrum",
                f"A spectrum named {name!r} already exists. Choose a different name.",
                parent=parent or self,
            )
            return None
        try:
            library_storage.save_spectrum(name, wlen, intensity)
        except Exception as exc:
            messagebox.showerror(
                "Spectrum", f"Could not save spectrum:\n{exc}", parent=parent or self
            )
            return None
        self._custom_spectra[name] = (
            np.asarray(wlen, dtype=float),
            np.asarray(intensity, dtype=float),
        )
        self._status.config(text=f"Added spectrum {name!r}.")
        pop = self._spectrum_library_popup
        if pop is not None:
            try:
                if pop.winfo_exists():
                    pop.append_entry_if_missing(name)
            except tk.TclError:
                pass
        return name

    def _finish_remove_known_stimulus(self, label: str) -> None:
        """Apply in-memory/UI updates after a stimulus has been removed from disk."""
        self._custom_stimuli.pop(label, None)
        self._known_stimuli = list(LABELS) + sorted(self._custom_stimuli.keys())
        if self._selected_stimulus_var.get() == label:
            self._selected_stimulus_var.set(NO_STIMULUS_SELECTED_LABEL)
        self._status.config(text=f"Removed stimulus {label!r}.")
        self._refresh_stimulus_dropdown_values()
        if self._stimulus_library_popup is not None:
            try:
                if self._stimulus_library_popup.winfo_exists():
                    self._stimulus_library_popup._refresh_list()
            except (tk.TclError, AttributeError):
                pass

    def start_remove_known_stimulus(
        self,
        label: str,
        parent=None,
        on_success=None,
    ) -> None:
        """Delete a custom stimulus with a progress dialog (disk I/O in background)."""
        if label in LABELS:
            messagebox.showinfo(
                "Delete stimulus",
                "Manuscript stimuli cannot be deleted.",
                parent=parent or self,
            )
            return
        if label not in self._custom_stimuli:
            messagebox.showerror(
                "Delete stimulus",
                f"Could not find stimulus {label!r}.",
                parent=parent or self,
            )
            return

        dlg = StimulusDeleteProgressDialog(parent or self, label)
        progress_q: Queue[tuple[str, int]] = Queue()
        result_q: Queue[tuple[str, str | None]] = Queue()

        def worker() -> None:
            try:
                progress_q.put(("Removing stimulus…", 30))
                library_storage.delete_stimulus(label)
                progress_q.put(("Complete", 100))
                result_q.put(("ok", None))
            except Exception as exc:
                result_q.put(("error", str(exc)))

        def poll() -> None:
            try:
                while True:
                    message, percent = progress_q.get_nowait()
                    dlg.set_progress(message, percent)
            except Empty:
                pass
            try:
                kind, payload = result_q.get_nowait()
            except Empty:
                dlg.after(50, poll)
                return
            try:
                dlg.grab_release()
            except tk.TclError:
                pass
            dlg.destroy()
            if kind == "error":
                messagebox.showerror(
                    "Delete stimulus",
                    f"Could not delete stimulus:\n{payload}",
                    parent=parent or self,
                )
                return
            self._finish_remove_known_stimulus(label)
            if on_success is not None:
                on_success()

        threading.Thread(target=worker, daemon=True).start()
        poll()

    def _remove_known_stimulus(self, label: str) -> bool:
        if label in LABELS:
            messagebox.showinfo(
                "Delete stimulus",
                "Manuscript stimuli cannot be deleted.",
                parent=self,
            )
            return False
        if label not in self._custom_stimuli:
            return False
        try:
            library_storage.delete_stimulus(label)
        except OSError:
            messagebox.showerror(
                "Delete stimulus",
                f"Could not delete stimulus {label!r} from disk.",
                parent=self,
            )
            return False
        self._finish_remove_known_stimulus(label)
        return True

    def _reload_stimulus_library_from_disk(self) -> None:
        """Reload custom stimuli from the user library and refresh known list."""
        self._custom_stimuli = library_storage.load_all_stimuli()
        self._known_stimuli = list(LABELS) + sorted(self._custom_stimuli.keys())
        self._refresh_stimulus_dropdown_values()

    def _register_custom_stimulus(self, name: str, spec: dict) -> None:
        """Persist a stimulus spec to disk and add it to the run dropdown."""
        library_storage.save_stimulus(spec)
        self._custom_stimuli[name] = spec
        if name not in self._known_stimuli:
            self._known_stimuli = list(LABELS) + sorted(self._custom_stimuli.keys())
        self._refresh_stimulus_dropdown_values()
        if self._stimulus_library_popup is not None:
            try:
                if self._stimulus_library_popup.winfo_exists():
                    self._stimulus_library_popup._refresh_list()
            except (tk.TclError, AttributeError):
                pass
        self._status.config(text=f"Saved stimulus {name!r} to library.")

    def _infer_builtin_spectrum_name(self, profile: np.ndarray) -> str:
        """Infer which built-in spectrum a profile most closely matches."""
        env = self._env
        if env is None:
            raise ValueError("Environment not loaded; cannot infer manuscript spectrum.")
        p = np.asarray(profile, dtype=float)
        p = np.clip(p, 0.0, None)
        if np.allclose(p, 0.0):
            return "Dark"

        candidates = {
            "440 nm": np.asarray(env["intensity440"], dtype=float),
            "560 nm": np.asarray(env["intensity560"], dtype=float),
            "Xenon": np.asarray(env["xenon"]["photons"], dtype=float),
            "Xenon (eye)": np.asarray(env["xenon_eye"]["photons"], dtype=float),
        }
        p_norm = p / max(float(np.max(p)), 1e-12)
        best_name = "Xenon"
        best_score = float("inf")
        for name, c in candidates.items():
            c = np.clip(c, 0.0, None)
            c_norm = c / max(float(np.max(c)), 1e-12)
            score = float(np.mean(np.abs(p_norm - c_norm)))
            if score < best_score:
                best_name = name
                best_score = score
        return best_name

    def _build_stimulus_spec_from_label(self, label: str) -> dict:
        """Convert a hardcoded manuscript label into an editable stimulus spec."""
        if self._env is None:
            raise ValueError("Spectra are still loading.")
        if label not in LABELS:
            raise ValueError(f"{label!r} is not a manuscript stimulus label.")

        stim = _build_stimulus(label, self._env)
        timings = np.asarray(stim["timings"], dtype=float).reshape(-1)
        if timings.size % 2 != 0:
            raise ValueError(
                f"Stimulus {label!r} has malformed timing array with odd length."
            )
        intensities = np.asarray(stim["intensities"], dtype=float)
        wlen = np.asarray(stim["wlen"], dtype=float)
        n_intervals = timings.size // 2
        if intensities.shape[0] != n_intervals or wlen.shape[0] != n_intervals:
            raise ValueError(
                f"Stimulus {label!r} interval mismatch: timings={n_intervals}, "
                f"intensities={intensities.shape[0]}, wlen={wlen.shape[0]}."
            )

        pulse_intervals: list[dict] = []
        for idx in range(n_intervals):
            t_start = float(timings[2 * idx])
            t_end = float(timings[2 * idx + 1])
            duration = t_end - t_start
            if duration <= 0:
                continue
            profile = np.asarray(intensities[idx], dtype=float)
            wrow = np.asarray(wlen[idx], dtype=float)
            intensity = float(np.trapezoid(profile, wrow))
            spectrum_name = self._infer_builtin_spectrum_name(profile)
            pulse_intervals.append(
                {
                    "start": t_start,
                    "end": t_end,
                    "spectrum": spectrum_name,
                    "intensity": max(0.0, intensity),
                    "duration": duration,
                }
            )

        if not pulse_intervals:
            raise ValueError(f"Stimulus {label!r} produced no positive-duration intervals.")

        total_duration = float(stim.get("ti", timings[-1]))
        blocks: list[dict] = []
        cursor = 0.0
        for pulse in pulse_intervals:
            gap = float(pulse["start"]) - cursor
            if gap > 0:
                blocks.append(
                    library_storage.make_interval_block("Dark", 0.0, gap)
                )
            blocks.append(
                library_storage.make_interval_block(
                    str(pulse["spectrum"]),
                    float(pulse["intensity"]),
                    float(pulse["duration"]),
                )
            )
            cursor = float(pulse["end"])

        trailing_dark = total_duration - cursor
        if trailing_dark > 0:
            blocks.append(
                library_storage.make_interval_block("Dark", 0.0, trailing_dark)
            )

        return {
            "version": _STIMULUS_SCHEMA_VERSION,
            "name": label,
            "total_duration": total_duration,
            "blocks": blocks,
        }

    def _build_run_protocol(self, label: str) -> dict:
        """Return a ``predictMelanopsin``-ready protocol dict for ``label``."""
        if label in LABELS:
            return _build_stimulus(label, self._env)
        if label in self._custom_stimuli:
            return _build_custom_protocol(
                self._custom_stimuli[label], self._env, self._custom_spectra
            )
        raise ValueError(f"Unknown stimulus label: {label!r}")

    def _stimuli_using_spectrum(self, name: str) -> list[str]:
        """Return custom stimulus names whose blocks reference ``name``."""
        self._reload_stimulus_library_from_disk()
        target = str(name).strip()
        return sorted(
            stim_name
            for stim_name, spec in self._custom_stimuli.items()
            if library_storage.spec_references_spectrum(spec, target)
        )

    def _apply_spectrum_delete_resolution(
        self,
        name: str,
        affected: list[str],
        result: tuple[str, str | None],
    ) -> bool:
        """Resolve stimulus references before deleting spectrum ``name``.

        ``result`` is ``(\"delete_stimuli\", None)`` or ``(\"replace\", replacement)``.
        Returns False if resolution fails (spectrum should not be deleted).
        """
        action, replacement = result
        try:
            if action == "delete_stimuli":
                for label in list(affected):
                    if not self._remove_known_stimulus(label):
                        if label in self._custom_stimuli:
                            messagebox.showerror(
                                "Delete Spectrum",
                                f"Could not delete stimulus protocol {label!r}. "
                                "The spectrum was not deleted.",
                                parent=self,
                            )
                            return False
                self._status.config(
                    text=(
                        f"Removed {len(affected)} stimulus protocol(s) that used "
                        f"spectrum {name!r}."
                    )
                )
            elif action == "replace":
                if not replacement:
                    messagebox.showerror(
                        "Delete Spectrum",
                        "No replacement spectrum was selected.",
                        parent=self,
                    )
                    return False
                for label in affected:
                    spec = self._custom_stimuli.get(label)
                    if spec is None:
                        continue
                    new_spec = library_storage.replace_spectrum_ref(
                        spec, name, replacement
                    )
                    self._register_custom_stimulus(label, new_spec)
                self._status.config(
                    text=(
                        f"Replaced spectrum {name!r} with {replacement!r} in "
                        f"{len(affected)} stimulus protocol(s)."
                    )
                )
            else:
                messagebox.showerror(
                    "Delete Spectrum",
                    f"Unknown resolution action {action!r}.",
                    parent=self,
                )
                return False
        except Exception as exc:
            messagebox.showerror(
                "Delete Spectrum",
                f"Could not resolve stimulus references:\n{exc}\n"
                "The spectrum was not deleted.",
                parent=self,
            )
            return False

        if self._stimulus_creator_dialog is not None:
            try:
                if self._stimulus_creator_dialog.winfo_exists():
                    messagebox.showinfo(
                        "Delete Spectrum",
                        "The open Stimulus Builder may still reference the deleted "
                        "spectrum in its unsaved grid. Review its intervals before saving.",
                        parent=self,
                    )
            except tk.TclError:
                pass
        return True

    def _delete_spectrum(self, name: str) -> bool:
        if name in _SPECTRUM_LIBRARY_ENTRIES:
            messagebox.showinfo(
                "Delete Spectrum",
                "Manuscript spectra cannot be deleted or modified.",
                parent=self,
            )
            return False
        if name in self._custom_spectra:
            self._custom_spectra.pop(name, None)
            library_storage.delete_spectrum(name)
            self._status.config(text=f"Deleted custom spectrum {name!r}.")
            return True
        messagebox.showerror(
            "Delete Spectrum",
            f"Could not find spectrum {name!r}.",
            parent=self,
        )
        return False

    def _on_file_import_spectrum(self) -> None:
        """Import a spectrum file into the user spectrum library."""
        self._add_spectrum_via_dialog(parent=self)

    def _stimulus_name_is_taken(self, name: str) -> bool:
        return name in LABELS or name in self._custom_stimuli

    def _prompt_unique_stimulus_import_name(self, default: str = "") -> str | None:
        while True:
            name = simpledialog.askstring(
                "Import Stimulus",
                "Enter a name for this stimulus:",
                initialvalue=default,
                parent=self,
            )
            if name is None:
                return None
            name = name.strip()
            if not name:
                messagebox.showerror(
                    "Import Stimulus", "Name cannot be empty.", parent=self
                )
                continue
            if self._stimulus_name_is_taken(name):
                messagebox.showerror(
                    "Import Stimulus",
                    f"A stimulus named {name!r} already exists. "
                    "Choose a different name.",
                    parent=self,
                )
                continue
            return name

    def _on_file_import_stimulus(self) -> None:
        """Import a v2 stimulus JSON into the user stimulus library."""
        path = filedialog.askopenfilename(
            parent=self,
            title="Import Stimulus",
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            spec = _load_stimulus_spec_from_path(path)
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            messagebox.showerror(
                "Import Stimulus",
                f"Could not read stimulus:\n{exc}",
                parent=self,
            )
            return

        name = str(spec.get("name", "")).strip()
        if not name:
            messagebox.showerror(
                "Import Stimulus",
                "Stimulus file is missing a non-empty name.",
                parent=self,
            )
            return

        if self._stimulus_name_is_taken(name):
            rename = messagebox.askyesno(
                "Import Stimulus",
                f"A stimulus named {name!r} already exists.\n\n"
                "Choose a different name to import a copy?",
                parent=self,
            )
            if not rename:
                return
            new_name = self._prompt_unique_stimulus_import_name(default=name)
            if new_name is None:
                return
            spec = dict(spec)
            spec["name"] = new_name
            name = new_name

        try:
            for block in spec.get("blocks", []):
                if isinstance(block, dict) and block.get("type") == "interval":
                    self._ensure_custom_spectrum_persisted(
                        str(block.get("spectrum_ref", "")).strip()
                    )
            self._register_custom_stimulus(name, spec)
        except (ValueError, OSError) as exc:
            messagebox.showerror("Import Stimulus", str(exc), parent=self)
            return

        self._set_selected_stimulus(name)
        messagebox.showinfo(
            "Import Stimulus",
            f"Imported stimulus {name!r} into the stimulus library.",
            parent=self,
        )

    def _on_file_load_stimulus(self) -> None:
        path = filedialog.askopenfilename(
            title="Load stimulus preset",
            filetypes=[
                ("JSON", "*.json"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            messagebox.showerror("Load stimulus", f"Could not read JSON:\n{exc}")
            return
        label = _parse_stimulus_label_from_json(data)
        if label is None:
            messagebox.showerror(
                "Load stimulus",
                'JSON must contain a string "label" or "stimulus" key.',
            )
            return
        if label not in self._known_stimuli:
            messagebox.showerror(
                "Load stimulus",
                f"Unknown label {label!r}. Use one of:\n{', '.join(self._known_stimuli)}",
            )
            return
        self._set_selected_stimulus(label)
        self._status.config(text=f"Loaded stimulus from {path}.")

    def _on_file_export_data(self) -> None:
        if self._last_result is None or self._last_label is None:
            messagebox.showwarning(
                "Export Model Run", "Run the model first, then export."
            )
            return
        export_dir = Path(DEFAULT_AUTOSAVE_DATA_DIR)
        try:
            export_dir.mkdir(parents=True, exist_ok=True)
            initial_dir = str(export_dir.resolve())
        except OSError:
            parent = export_dir.parent
            initial_dir = (
                str(parent.resolve()) if parent.exists() else str(Path.cwd())
            )
        path = filedialog.asksaveasfilename(
            title="Export Model Run",
            defaultextension=".csv",
            initialdir=initial_dir,
            filetypes=[
                ("CSV table", "*.csv"),
                ("NumPy archive", "*.npz"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        lower = path.lower()
        try:
            if lower.endswith(".npz"):
                _export_result_npz(path, self._last_result, self._last_label)
            else:
                if not lower.endswith(".csv"):
                    path = path + ".csv"
                _export_result_csv(path, self._last_result)
        except (OSError, ValueError) as exc:
            messagebox.showerror("Export Model Run", str(exc))
            return
        self._status.config(text=f"Exported data to {path}.")

    def _start_progress_poll(self) -> None:
        self._drain_progress_queue()
        self._progress_poll_id = self.after(75, self._start_progress_poll)

    def _stop_progress_poll(self) -> None:
        if self._progress_poll_id is not None:
            self.after_cancel(self._progress_poll_id)
            self._progress_poll_id = None
        self._drain_progress_queue()

    def _drain_progress_queue(self) -> None:
        latest: tuple[int, int] | None = None
        while True:
            try:
                latest = self._progress_queue.get_nowait()
            except Empty:
                break
        if latest is None:
            return
        current_step, total_steps = latest
        if total_steps <= 0:
            self._progress.configure(value=0)
            self._status.config(text="Running model… 0% (0/0)")
            return
        percent = int(round((100.0 * current_step) / total_steps))
        self._progress.configure(value=max(0, min(100, percent)))
        self._status.config(
            text=f"Running model… {percent}% ({current_step}/{total_steps})"
        )


def main(splash: StartupSplash | None = None) -> None:
    app = ManuscriptSimApp(splash=splash)
    app.minsize(640, 480)
    try:
        app.mainloop()
    finally:
        try:
            plt.close("all")
        except Exception:
            pass


if __name__ == "__main__":
    main()
