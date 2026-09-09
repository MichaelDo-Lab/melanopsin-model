# Packaging the Melanopsin Model GUI

This folder contains everything needed to turn the desktop GUI
(`myutils/melanopsin_gui.py`) into a one-click Windows executable or a macOS
`.app` bundle.

| File | Purpose |
| --- | --- |
| `melanopsin_gui.spec` | PyInstaller build recipe (entry point, hidden imports, excludes, icon, output name). On macOS it emits an onedir `.app`; on Windows a one-file `.exe`. |
| `build_windows.ps1` | Convenience script: makes a clean venv, installs deps, runs PyInstaller. |
| `build_macos.sh` | Same for macOS; also zips the `.app` for GitHub Releases. Builds the host Mac's architecture (`arm64` or `x86_64`). |
| `requirements-build.txt` | Pinned build dependencies (app runtime deps + PyInstaller). |
| `app_icon.ico` | Optional. If present, used as the Windows executable icon. |
| `app_icon.icns` | Optional. If present, used as the macOS app icon. |

## How distribution works

The packaged app is **not** fully self-contained. It is designed to live inside a
cloned copy of this repository and reads `data/` and writes `outputs/` relative
to the repository, not relative to itself. This keeps the large spectral assets
in version control (so they can be updated independently) and keeps the
binary small.

Path resolution is centralized in `myutils/app_paths.py`. When frozen, it locates
the repository by walking up from the executable's location until it finds a
folder containing `data/` and `myutils/`. That means the same `.exe` or `.app`
works whether it sits in the repository root or in `dist/`. On macOS the
binary lives at `Something.app/Contents/MacOS/...`; walking up from there still
reaches the clone.

## Build locally (Windows)

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

## Build locally (macOS)

From the repository root:

```bash
bash packaging/build_macos.sh --clean
```

This produces `dist/MelanopsinModel-v<version>-macos-<arch>.app` and
`dist/MelanopsinModel-v<version>-macos-<arch>.zip`, where `<arch>` is `arm64`
(Apple silicon) or `x86_64` (Intel) matching the Mac you built on. A local
build only produces one architecture; GitHub Actions builds both.

To build with your current Python environment (no isolated venv):

```bash
bash packaging/build_macos.sh --skip-venv
```

The first time users open an unsigned `.app` downloaded from the internet,
macOS Gatekeeper may block it. Right-click the app, choose **Open**, and
confirm, or run `xattr -dr com.apple.quarantine MelanopsinModel-v*.app`.

## Cutting a new release

1. Bump `__version__` in `myutils/_version.py`.
2. Commit and tag, e.g. `git tag v0.2.0 && git push --tags`.
3. The GitHub Actions workflow (`.github/workflows/build-release.yml`) builds
   the Windows executable, an Apple silicon `.app`, and an Intel `.app`, then
   attaches all three to a GitHub Release for that tag.

Release artifacts:

| File | For |
| --- | --- |
| `MelanopsinModel-v<version>.exe` | Windows |
| `MelanopsinModel-v<version>-macos-arm64.zip` | Apple silicon Macs |
| `MelanopsinModel-v<version>-macos-x86_64.zip` | Intel Macs (also runs on Apple silicon via Rosetta 2) |

The Intel job uses GitHub's `macos-15-intel` runner. GitHub plans to drop
hosted Intel macOS runners around Fall 2027; the Apple silicon job is the
long-term path.

## Adding the GUI's dependencies

If the GUI grows to import a new third-party package, add it to
`requirements.txt`. If PyInstaller fails to detect it automatically (common for
packages imported dynamically or via plugins), also add it to `hiddenimports` in
`melanopsin_gui.spec`. Conversely, large notebook-only dependencies should be
added to `excludes` and imported lazily so the executable stays lean.
