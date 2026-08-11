"""Minimal startup splash window (tkinter only).

Used by ``run_melanopsin_gui.py`` during heavy imports and by
``ManuscriptSimApp`` while the interface and spectral assets load. Intentionally
free of heavy dependencies so it can appear before ``matplotlib``/``numpy`` are
imported.

The splash is its own ``tk.Tk()`` root (not a ``Toplevel``) so it stays visible
even while the main application window is withdrawn.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk


class StartupSplash:
    """Standalone loading window with status text and a progress bar."""

    def __init__(
        self,
        *,
        title: str = "Melanopsin Model",
        version: str = "",
    ) -> None:
        self._root = tk.Tk()
        self._root.title(title)
        self._root.resizable(False, False)
        self._root.protocol("WM_DELETE_WINDOW", lambda: None)

        outer = ttk.Frame(self._root, padding=24)
        outer.pack()

        ttk.Label(outer, text=title, font=("", 14, "bold")).pack()
        if version:
            ttk.Label(outer, text=f"v{version}", foreground="#555555").pack(
                pady=(0, 12)
            )

        self._message = tk.StringVar(self._root, value="Starting…")
        ttk.Label(outer, textvariable=self._message, width=44).pack(pady=(0, 8))

        self._progress = ttk.Progressbar(
            outer, mode="determinate", length=340, maximum=100
        )
        self._progress.pack()
        self._indeterminate = False

        self._center_on_screen()
        self._root.update_idletasks()

    def _center_on_screen(self) -> None:
        self._root.update_idletasks()
        w = self._root.winfo_width()
        h = self._root.winfo_height()
        sw = self._root.winfo_screenwidth()
        sh = self._root.winfo_screenheight()
        x = max(0, (sw - w) // 2)
        y = max(0, (sh - h) // 2)
        self._root.geometry(f"+{x}+{y}")

    def set_message(self, message: str) -> None:
        self._message.set(message)
        self._root.update_idletasks()

    def set_progress(self, percent: float | None) -> None:
        """Set determinate progress (0-100), or ``None`` for an animated bar."""
        if percent is None:
            if not self._indeterminate:
                self._progress.configure(mode="indeterminate")
                self._progress.start(12)
                self._indeterminate = True
        else:
            if self._indeterminate:
                self._progress.stop()
                self._progress.configure(mode="determinate")
                self._indeterminate = False
            self._progress["value"] = max(0.0, min(100.0, float(percent)))
        self._root.update_idletasks()

    def update(self, message: str, percent: float | None = None) -> None:
        self.set_message(message)
        self.set_progress(percent)
        self._root.update()

    def close(self) -> None:
        if self._indeterminate:
            try:
                self._progress.stop()
            except tk.TclError:
                pass
        try:
            self._root.destroy()
        except tk.TclError:
            pass
