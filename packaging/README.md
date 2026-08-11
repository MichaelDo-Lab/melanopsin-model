# Packaging the Melanopsin Model GUI

This folder contains everything needed to turn the desktop GUI
(`myutils/melanopsin_gui.py`) into a single, one-click Windows executable.

| File | Purpose |
| --- | --- |
| `melanopsin_gui.spec` | PyInstaller build recipe (entry point, hidden imports, excludes, icon, output name). |
| `build_windows.ps1` | Convenience script: makes a clean venv, installs deps, runs PyInstaller. |
| `requirements-build.txt` | Pinned build dependencies (app runtime deps + PyInstaller). |
| `app_icon.ico` | Optional. If present, used as the executable icon. |

## How distribution works

The executable is **not** fully self-contained. It is designed to live inside a
cloned copy of this repository and reads `data/` and writes `outputs/` relative
to the repository, not relative to itself. This keeps the large spectral assets
in version control (so they can be updated independently) and keeps the
executable small.

Path resolution is centralized in `myutils/app_paths.py`. When frozen, it locates
the repository by walking up from the executable's location until it finds a
folder containing `data/` and `myutils/`. That means the same `.exe` works
whether it sits in the repository root or in `dist/`.

## Build locally

From the repository root, in PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File packaging\build_windows.ps1 -Clean
```

This produces `dist/MelanopsinModel-v<version>.exe`. The version string comes
from `myutils/_version.py`.

To build with your current Python environment (no isolated venv):

```powershell
powershell -ExecutionPolicy Bypass -File packaging\build_windows.ps1 -SkipVenv
```

Or invoke PyInstaller directly:

```powershell
pip install -r packaging\requirements-build.txt
pyinstaller --noconfirm --clean packaging\melanopsin_gui.spec
```

## Cutting a new release

1. Bump `__version__` in `myutils/_version.py`.
2. Commit and tag, e.g. `git tag v0.2.0 && git push --tags`.
3. The GitHub Actions workflow (`.github/workflows/build-release.yml`) builds the
   Windows executable and attaches it to a GitHub Release for that tag.

## Adding the GUI's dependencies

If the GUI grows to import a new third-party package, add it to
`requirements.txt`. If PyInstaller fails to detect it automatically (common for
packages imported dynamically or via plugins), also add it to `hiddenimports` in
`melanopsin_gui.spec`. Conversely, large notebook-only dependencies should be
added to `excludes` and imported lazily so the executable stays lean.
