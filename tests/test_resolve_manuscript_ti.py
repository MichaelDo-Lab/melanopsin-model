"""Unit tests for manuscript simulation-duration override resolution.

Loads ``resolve_manuscript_ti`` from the GUI source by AST so ``myutils.__init__``
(heavy scipy/pandas) is not imported — same idea as other lightweight tests.
"""

from __future__ import annotations

import ast
from pathlib import Path

GUI_PATH = Path(__file__).resolve().parents[1] / "myutils" / "melanopsin_gui.py"


def _load_resolve_manuscript_ti():
    source = GUI_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(GUI_PATH))
    keep: list[ast.stmt] = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "LABELS":
                    keep.append(node)
        elif (
            isinstance(node, ast.FunctionDef)
            and node.name == "resolve_manuscript_ti"
        ):
            keep.append(node)
    module = ast.Module(body=keep, type_ignores=[])
    ast.fix_missing_locations(module)
    namespace: dict[str, object] = {}
    exec(compile(module, str(GUI_PATH), "exec"), namespace)
    return namespace["resolve_manuscript_ti"]


resolve_manuscript_ti = _load_resolve_manuscript_ti()


def test_leftover_custom_duration_does_not_change_manuscript_ti() -> None:
    """A custom builder duration must not clip a different manuscript protocol."""
    native_ti = 3091.0
    assert (
        resolve_manuscript_ti(
            native_ti,
            override=560.0,
            override_label="test",
            run_label="440-440",
        )
        == native_ti
    )


def test_matching_manuscript_label_with_longer_override_extends_ti() -> None:
    """Trailing dark in the builder may extend the same manuscript protocol."""
    native_ti = 3091.0
    assert (
        resolve_manuscript_ti(
            native_ti,
            override=3200.0,
            override_label="440-440",
            run_label="440-440",
        )
        == 3200.0
    )


def test_matching_label_with_shorter_override_leaves_native_ti() -> None:
    """Override must not shrink ti; manuscript timings are not truncated."""
    native_ti = 3091.0
    assert (
        resolve_manuscript_ti(
            native_ti,
            override=560.0,
            override_label="440-440",
            run_label="440-440",
        )
        == native_ti
    )


def test_no_override_leaves_native_ti() -> None:
    native_ti = 3091.0
    assert (
        resolve_manuscript_ti(
            native_ti,
            override=None,
            override_label=None,
            run_label="440-440",
        )
        == native_ti
    )


def test_custom_run_label_ignores_override() -> None:
    """Custom stimuli keep their own spec duration even if an override is set."""
    native_ti = 560.0
    assert (
        resolve_manuscript_ti(
            native_ti,
            override=3200.0,
            override_label="440-440",
            run_label="test",
        )
        == native_ti
    )
