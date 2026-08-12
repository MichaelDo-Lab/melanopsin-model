"""Static guard: tkinter Vars and nametofont must pass an explicit master/root.

``run_melanopsin_gui.py`` creates a splash ``Tk`` before the app, so destroying
the splash clears ``tkinter._default_root``. Dialogs that rely on the implicit
default then raise ``RuntimeError``. This AST check catches regressions without
needing a display or heavy imports.
"""

from __future__ import annotations

import ast
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
GUI_PATH = REPO_ROOT / "myutils" / "melanopsin_gui.py"

_VAR_NAMES = frozenset({"StringVar", "BooleanVar", "IntVar", "DoubleVar"})


def _attr_chain(node: ast.AST) -> list[str]:
    parts: list[str] = []
    cur: ast.AST | None = node
    while isinstance(cur, ast.Attribute):
        parts.append(cur.attr)
        cur = cur.value
    if isinstance(cur, ast.Name):
        parts.append(cur.id)
    return list(reversed(parts))


def _call_has_master(call: ast.Call) -> bool:
    """True if master is given as a keyword or as the first positional arg."""
    for kw in call.keywords:
        if kw.arg == "master":
            return True
    return len(call.args) >= 1


def _call_has_root(call: ast.Call) -> bool:
    return any(kw.arg == "root" for kw in call.keywords)


def test_melanopsin_gui_vars_and_fonts_use_explicit_master() -> None:
    source = GUI_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(GUI_PATH))
    failures: list[str] = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        chain = _attr_chain(node.func)
        if not chain:
            continue
        name = chain[-1]
        lineno = getattr(node, "lineno", "?")

        if name in _VAR_NAMES:
            if not _call_has_master(node):
                failures.append(
                    f"{GUI_PATH.name}:{lineno}: {'.'.join(chain)}(...) "
                    "missing explicit master"
                )
        elif name == "nametofont":
            if not _call_has_root(node):
                failures.append(
                    f"{GUI_PATH.name}:{lineno}: {'.'.join(chain)}(...) "
                    "missing explicit root="
                )

    assert not failures, "Implicit default-root usage:\n" + "\n".join(failures)
