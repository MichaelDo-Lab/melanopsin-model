"""Filesystem layout helpers for the Melanopsin Model application.

This module is the single source of truth for where the app reads its bundled
assets (``data/``) and writes user output (``outputs/``). Every other module
should resolve those locations through here instead of using relative paths or
``__file__`` directly, so behavior stays identical whether the code runs from a
source checkout (notebooks, ``python -m myutils.melanopsin_gui``) or from a
frozen PyInstaller executable.

Distribution contract for the ``.exe`` / ``.app``
-------------------------------------------------
The packaged app is intentionally *not* fully self-contained: it ships alongside
a cloned copy of this repository and reads ``data/`` and writes ``outputs/``
relative to that repository -- not relative to the temporary folder PyInstaller
unpacks itself into. The intended workflow for academics is:

    git clone <repo>
    # double-click the .exe or .app placed in the clone (or in dist/)

``project_root()`` therefore resolves the repository, not the bundle:

* Frozen build: start at the executable's directory and walk upward looking for
  a folder that contains the expected asset directories (so the same binary
  works whether it lives in the repo root or in ``dist/``). On macOS the
  executable sits at ``Something.app/Contents/MacOS/...``; walking up from there
  still reaches the clone. If none is found we fall back to the executable's
  own directory.
* Source checkout: the parent of the ``myutils/`` package.
"""

from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

# Directories that, together, identify a melanopsin-model checkout. Keep this in
# sync with the folders committed to the repository root. Use folders that are
# always present after a fresh clone (empty dirs are not tracked by git).
_ROOT_MARKERS = ("data", "myutils")


def is_frozen() -> bool:
    """Return ``True`` when running inside a PyInstaller-built executable."""
    return bool(getattr(sys, "frozen", False))


def _looks_like_root(path: Path) -> bool:
    """True if ``path`` contains all of the marker asset directories."""
    return all((path / marker).is_dir() for marker in _ROOT_MARKERS)


@lru_cache(maxsize=1)
def project_root() -> Path:
    """Absolute path to the repository the app should read from / write to."""
    if is_frozen():
        start = Path(sys.executable).resolve().parent
        for candidate in (start, *start.parents):
            if _looks_like_root(candidate):
                return candidate
        # Asset folders not found next to the exe; best effort is its own dir.
        return start
    # Source checkout: <root>/myutils/app_paths.py -> <root>
    return Path(__file__).resolve().parent.parent


def data_dir() -> Path:
    """Path to the repository ``data/`` directory (bundled spectral assets)."""
    return project_root() / "data"


def manuscript_dir() -> Path:
    """Path to ``data/manuscript/`` (container for recordings, spectra, irradiance data)."""
    return data_dir() / "manuscript"


def manuscript_spectra_dir() -> Path:
    """Path to ``data/manuscript/spectra/`` (bundled measured spectra)."""
    return manuscript_dir() / "spectra"


def manuscript_recordings_dir() -> Path:
    """Path to ``data/manuscript/recordings/`` (example melanopsin recording CSVs)."""
    return manuscript_dir() / "recordings"


def configs_dir() -> Path:
    """Path to the legacy ``configs/`` directory (pre-overhaul stimulus JSON).

    Runtime storage lives under ``data/user_library/``. This path is retained only
    so startup can delete an old ``configs/config.json`` if one is still present.
    """
    return project_root() / "configs"


def outputs_dir() -> Path:
    """Path to the repository ``outputs/`` directory (user-generated output)."""
    return project_root() / "outputs"


def ensure_dir(path: Path) -> Path:
    """Create ``path`` (and parents) if missing and return it for chaining."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
