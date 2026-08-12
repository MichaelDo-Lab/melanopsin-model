"""Live regression: launcher-order startup must leave a usable default root.

``run_melanopsin_gui.py`` (and the released ``.exe``) create ``StartupSplash``
before ``ManuscriptSimApp``. Destroying that splash used to clear
``tkinter._default_root``, so post-startup dialogs crashed. This test reproduces
that launch order and asserts the dialogs open.
"""

from __future__ import annotations

import time
import tkinter as tk

import pytest

from myutils._version import __version__
from myutils.melanopsin_gui import (
    AutoSaveOptionsDialog,
    DataComparatorDialog,
    ManuscriptSimApp,
    ModelConfigDialog,
    PredictionComparisonDialog,
    SpectrumBuilderDialog,
    SpectrumLibraryPopup,
    StimulusBuilderDialog,
    StimulusLibraryPopup,
)
from myutils.startup_splash import StartupSplash


def _display_available() -> bool:
    try:
        root = tk.Tk()
    except tk.TclError:
        return False
    try:
        root.destroy()
    except tk.TclError:
        pass
    return True


pytestmark = pytest.mark.skipif(
    not _display_available(),
    reason="No Tk display available",
)


def _run_until_startup_finished(app: ManuscriptSimApp, *, timeout_s: float = 120.0) -> None:
    """Enter ``mainloop`` until the splash closes (required for threaded ``after``)."""
    deadline = time.monotonic() + timeout_s
    timed_out = False

    def _poll() -> None:
        nonlocal timed_out
        if app._startup_splash is None:
            app.quit()
            return
        if time.monotonic() > deadline:
            timed_out = True
            app.quit()
            return
        app.after(50, _poll)

    app.after(50, _poll)
    app.mainloop()
    if timed_out:
        raise TimeoutError(
            "Startup splash did not close within "
            f"{timeout_s:.0f}s (spectra load may have hung)"
        )


def test_launcher_order_keeps_default_root_and_opens_dialogs() -> None:
    # Match run_melanopsin_gui.py: splash Tk is created before the app Tk.
    splash = StartupSplash(version=__version__)
    splash.update("Loading application modules…", None)
    app = ManuscriptSimApp(splash=splash)
    try:
        _run_until_startup_finished(app)
        assert tk._default_root is not None
        assert tk._default_root is app

        dialogs: list[tk.Toplevel] = []
        try:
            for factory in (
                ModelConfigDialog,
                StimulusBuilderDialog,
                StimulusLibraryPopup,
                SpectrumLibraryPopup,
                SpectrumBuilderDialog,
                DataComparatorDialog,
                PredictionComparisonDialog,
                AutoSaveOptionsDialog,
            ):
                dialogs.append(factory(app))
                app.update_idletasks()
        finally:
            for dialog in dialogs:
                try:
                    dialog.destroy()
                except tk.TclError:
                    pass
    finally:
        try:
            app.destroy()
        except tk.TclError:
            pass
