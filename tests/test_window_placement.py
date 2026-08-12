"""Unit tests for pure window-placement geometry (no display required)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

# Load the module by path so ``myutils.__init__`` (heavy scipy/pandas) is not
# imported — same pattern as other lightweight geometry/static tests.
_MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "myutils" / "window_placement.py"
)
_spec = importlib.util.spec_from_file_location("window_placement", _MODULE_PATH)
assert _spec is not None and _spec.loader is not None
_wp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_wp)
choose_window_position = _wp.choose_window_position


def test_edge_places_child_to_the_right_when_it_fits() -> None:
    """Child docks just outside the parent's right edge on the same work area."""
    work = (0, 0, 1920, 1080)
    x, y = choose_window_position(
        parent_x=100,
        parent_y=80,
        parent_w=800,
        parent_h=600,
        child_w=400,
        child_h=300,
        work_area=work,
        mode="edge",
    )
    assert x == 100 + 800 + 8  # parent right + gap
    assert y == 80
    assert x + 400 <= work[2]
    assert y + 300 <= work[3]


def test_maximized_parent_falls_back_to_inner_corner() -> None:
    """When no outside slot fits, place inside the parent (same monitor)."""
    # Parent fills the work area; child cannot sit beside it.
    work = (0, 0, 1920, 1080)
    x, y = choose_window_position(
        parent_x=0,
        parent_y=0,
        parent_w=1920,
        parent_h=1080,
        child_w=400,
        child_h=300,
        work_area=work,
        mode="edge",
    )
    # Top-right inner corner: parent_w - child_w - inset, inset
    assert x == 1920 - 400 - 24
    assert y == 24
    assert x >= work[0] and y >= work[1]
    assert x + 400 <= work[2] and y + 300 <= work[3]


def test_two_tool_windows_avoid_same_slot() -> None:
    """Second window picks a different edge when the preferred slot is occupied."""
    work = (0, 0, 1920, 1080)
    first = choose_window_position(
        parent_x=200,
        parent_y=100,
        parent_w=700,
        parent_h=500,
        child_w=350,
        child_h=280,
        work_area=work,
        mode="edge",
    )
    # Mark first slot as occupied with the first child's size.
    occupied = [(first[0], first[1], 350, 280)]
    second = choose_window_position(
        parent_x=200,
        parent_y=100,
        parent_w=700,
        parent_h=500,
        child_w=350,
        child_h=280,
        work_area=work,
        mode="edge",
        occupied=occupied,
    )
    assert second != first
    # Second should still be fully on the work area.
    assert second[0] >= work[0] and second[1] >= work[1]
    assert second[0] + 350 <= work[2] and second[1] + 280 <= work[3]


def test_center_is_parent_midpoint() -> None:
    work = (0, 0, 1920, 1080)
    x, y = choose_window_position(
        parent_x=100,
        parent_y=50,
        parent_w=800,
        parent_h=600,
        child_w=400,
        child_h=200,
        work_area=work,
        mode="center",
    )
    assert x == 100 + (800 - 400) // 2
    assert y == 50 + (600 - 200) // 2


def test_center_clamps_when_child_larger_than_work_area() -> None:
    work = (0, 0, 800, 600)
    x, y = choose_window_position(
        parent_x=0,
        parent_y=0,
        parent_w=800,
        parent_h=600,
        child_w=900,
        child_h=700,
        work_area=work,
        mode="center",
    )
    # Cannot fully fit; clamp to work-area origin.
    assert x == 0
    assert y == 0


def test_outside_slot_rejected_when_it_would_spill_off_monitor() -> None:
    """Narrow leftover space on the right forces left (or inner) placement."""
    # Parent near the right edge: right-of-parent would spill past 1920.
    work = (0, 0, 1920, 1080)
    x, y = choose_window_position(
        parent_x=1600,
        parent_y=100,
        parent_w=300,
        parent_h=700,
        child_w=400,
        child_h=300,
        work_area=work,
        mode="edge",
    )
    # Right slot (1600+300+8=1908) cannot fit 400px. Left might work:
    # 1600 - 400 - 8 = 1192 — fits. Prefer left over inner corners.
    assert x == 1600 - 400 - 8
    assert y == 100
    assert x + 400 <= work[2]
