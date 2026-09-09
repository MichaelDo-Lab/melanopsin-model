#!/usr/bin/env python3
"""One-click launcher for the Melanopsin Model desktop GUI.

This is the script PyInstaller bundles into the distributed Windows
``.exe`` and macOS ``.app`` (see ``packaging/melanopsin_gui.spec``). It is also
the recommended way to start the app from a source checkout::

    python run_melanopsin_gui.py

Heavy imports are deferred until after a friendly error handler is in place so
that a missing asset or dependency surfaces as a dialog box rather than a silent
crash in the windowed (no-console) build.
"""

from __future__ import annotations

import sys
import traceback


def _ensure_std_streams() -> None:
    """Guarantee ``sys.stdout``/``sys.stderr`` are writable.

    A windowed (``console=False``) PyInstaller build runs with no console, so the
    interpreter sets ``sys.stdout`` and ``sys.stderr`` to ``None``. Any library
    that writes to them (for example tqdm's progress bar, or a stray ``print``)
    then raises ``'NoneType' object has no attribute 'write'``. Redirect them to
    the null device so those writes are harmless no-ops.
    """
    import io
    import os

    for name in ("stdout", "stderr"):
        if getattr(sys, name, None) is None:
            try:
                setattr(sys, name, open(os.devnull, "w", encoding="utf-8"))
            except Exception:
                setattr(sys, name, io.StringIO())


def _show_fatal_error(title: str, message: str) -> None:
    """Best-effort error dialog; falls back to stderr when no display exists."""
    try:
        import tkinter as tk
        from tkinter import messagebox

        root = tk.Tk()
        root.withdraw()
        messagebox.showerror(title, message)
        root.destroy()
    except Exception:
        print(f"{title}\n\n{message}", file=sys.stderr)


def main() -> int:
    _ensure_std_streams()
    splash = None
    try:
        from myutils._version import __version__
        from myutils.startup_splash import StartupSplash

        splash = StartupSplash(version=__version__)
        splash.update("Loading application modules…", None)

        from myutils.melanopsin_gui import main as gui_main
    except Exception:
        if splash is not None:
            splash.close()
        _show_fatal_error(
            "Melanopsin Model - startup error",
            "The application failed to load its modules.\n\n"
            "If you are running the packaged app (.exe or .app), make sure it "
            "sits inside the cloned repository (next to the data/ and myutils/ "
            "folders).\n\n"
            + traceback.format_exc(),
        )
        return 1

    try:
        gui_main(splash=splash)
    except Exception:
        if splash is not None:
            try:
                splash.close()
            except Exception:
                pass
        _show_fatal_error(
            "Melanopsin Model - unexpected error",
            "The application stopped unexpectedly.\n\n" + traceback.format_exc(),
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
