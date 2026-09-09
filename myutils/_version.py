"""Single source of truth for the application version.

The version is consumed by the PyInstaller spec (to name the Windows
``MelanopsinModel-v<version>.exe`` and the macOS
``MelanopsinModel-v<version>-macos-<arch>.app``) and is shown in the GUI's
window title. Bump this string and tag the commit (e.g. ``v0.2.0``) to cut a
new release.
"""

__version__ = "0.1.1"
