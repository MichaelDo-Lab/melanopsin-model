# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the Melanopsin Model desktop GUI.

Build (from the repository root):

    pyinstaller packaging/melanopsin_gui.spec

Produces a single windowed executable in ``dist/`` named
``MelanopsinModel-v<version>.exe`` (the extension is added automatically on
Windows). The executable is meant to live inside a cloned copy of this
repository: it reads ``data/`` and ``configs/`` and writes ``outputs/``
relative to its own location, so the large data assets are NOT bundled.
"""

import os
import sys

from PyInstaller.utils.hooks import collect_submodules

# ``SPECPATH`` is injected by PyInstaller and points at the directory holding
# this spec file (``<repo>/packaging``); the repo root is its parent.
PROJECT_ROOT = os.path.abspath(os.path.join(SPECPATH, os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from myutils._version import __version__  # noqa: E402

APP_NAME = f"MelanopsinModel-v{__version__}"
ENTRY_SCRIPT = os.path.join(PROJECT_ROOT, "run_melanopsin_gui.py")
ICON_PATH = os.path.join(PROJECT_ROOT, "packaging", "app_icon.ico")

def _collect_conda_runtime_dlls():
    """Bundle conda ``Library/bin`` DLLs that stdlib extensions need.

    A venv created from Anaconda/Miniconda still loads ``pyexpat.pyd`` from the
    base prefix, which depends on ``libexpat.dll`` in ``Library/bin``. PyInstaller
    does not treat that venv as a conda env (no ``conda-meta``), so those DLLs
    are skipped unless collected here. Harmless on a stock CPython install:
    the folder simply does not exist.
    """
    collected = []
    if os.name != "nt":
        return collected
    prefixes = []
    for prefix in (sys.prefix, getattr(sys, "base_prefix", sys.prefix)):
        if prefix and prefix not in prefixes:
            prefixes.append(prefix)
    dll_names = (
        "libexpat.dll",
        "expat.dll",
        "libffi.dll",
        "ffi.dll",
        "liblzma.dll",
        "zlib.dll",
        "libz.dll",
        "sqlite3.dll",
        "libsqlite3.dll",
    )
    seen = set()
    for prefix in prefixes:
        bin_dir = os.path.join(prefix, "Library", "bin")
        if not os.path.isdir(bin_dir):
            continue
        for name in dll_names:
            path = os.path.join(bin_dir, name)
            key = name.lower()
            if key not in seen and os.path.isfile(path):
                collected.append((path, "."))
                seen.add(key)
    return collected


datas = []
binaries = _collect_conda_runtime_dlls()

# numpy / scipy / pandas / matplotlib are handled by PyInstaller's built-in
# hooks (collecting their data + DLLs automatically); collect_all on them just
# bloats the build with test suites. We only spell out the dynamically imported
# bits that hooks miss: the Tk matplotlib backend and skopt's submodules.
hiddenimports = [
    "matplotlib.backends.backend_tkagg",
    "PyPDF2",
    "openpyxl",
    "joblib",
    "tqdm",
    "seaborn",
    # setuptools 71+ / pkg_resources on Python < 3.12:
    # pyi_rth_pkgres -> pkg_resources -> jaraco.context -> backports.tarfile
    "backports",
    "backports.tarfile",
    "jaraco",
    "jaraco.context",
    "jaraco.text",
    "jaraco.functools",
    "setuptools._vendor.backports",
    "setuptools._vendor.backports.tarfile",
]
hiddenimports += collect_submodules("skopt")
hiddenimports += collect_submodules("backports")

# Trim large, unused modules that get pulled in transitively. ipywidgets and the
# wider Jupyter stack are only used by the notebook UI; myutils.init imports them
# optionally, so excluding them here keeps the executable small without breaking
# the GUI.
excludes = [
    "IPython",
    "ipywidgets",
    "ipykernel",
    "notebook",
    "jupyter",
    "jupyter_client",
    "jupyter_core",
    "zmq",
    "tornado",
    "pytest",
    "numpy.tests",
    "scipy.tests",
    "pandas.tests",
    "matplotlib.tests",
]

block_cipher = None

a = Analysis(
    [ENTRY_SCRIPT],
    pathex=[PROJECT_ROOT],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name=APP_NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=ICON_PATH if os.path.exists(ICON_PATH) else None,
)
