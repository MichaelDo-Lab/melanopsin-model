"""Place Tk Toplevel windows near their parent on the same monitor.

Tk does not set Toplevel geometry by default; on multi-monitor Windows setups
new windows often open at the virtual-desktop origin. This module pins dialogs
to the parent's monitor, tucks tool windows along the parent's edges, and
centers modal prompts on their parent.
"""

from __future__ import annotations

import sys
from typing import Iterable, Literal, Sequence

PlacementMode = Literal["edge", "center"]

# Gap between parent outer edge and child when docking beside/above/below.
_EDGE_GAP = 8
# Inset from parent inner corners when the child must sit inside the parent.
_INNER_INSET = 24
# Cascade step when every preferred slot overlaps existing windows.
_CASCADE_STEP = 28


def _rect_overlap_area(
    ax: int, ay: int, aw: int, ah: int,
    bx: int, by: int, bw: int, bh: int,
) -> int:
    """Return intersection area of two axis-aligned rectangles (0 if none)."""
    left = max(ax, bx)
    top = max(ay, by)
    right = min(ax + aw, bx + bw)
    bottom = min(ay + ah, by + bh)
    if right <= left or bottom <= top:
        return 0
    return (right - left) * (bottom - top)


def _clamp_to_work_area(
    x: int, y: int, w: int, h: int,
    work: tuple[int, int, int, int],
) -> tuple[int, int]:
    """Clamp top-left ``(x, y)`` so the ``w``×``h`` window stays in ``work``.

    ``work`` is ``(left, top, right, bottom)`` in screen coordinates (right/bottom
    exclusive in the usual Win32 RECT sense: usable width = right - left).
    """
    left, top, right, bottom = work
    max_x = right - w
    max_y = bottom - h
    if max_x < left:
        x = left
    else:
        x = min(max(x, left), max_x)
    if max_y < top:
        y = top
    else:
        y = min(max(y, top), max_y)
    return x, y


def _fits_in_work_area(
    x: int, y: int, w: int, h: int,
    work: tuple[int, int, int, int],
) -> bool:
    left, top, right, bottom = work
    return x >= left and y >= top and x + w <= right and y + h <= bottom


def _total_overlap(
    x: int, y: int, w: int, h: int,
    occupied: Sequence[tuple[int, int, int, int]],
) -> int:
    total = 0
    for ox, oy, ow, oh in occupied:
        total += _rect_overlap_area(x, y, w, h, ox, oy, ow, oh)
    return total


def choose_window_position(
    *,
    parent_x: int,
    parent_y: int,
    parent_w: int,
    parent_h: int,
    child_w: int,
    child_h: int,
    work_area: tuple[int, int, int, int],
    mode: PlacementMode = "edge",
    occupied: Sequence[tuple[int, int, int, int]] = (),
    edge_gap: int = _EDGE_GAP,
    inner_inset: int = _INNER_INSET,
    cascade_step: int = _CASCADE_STEP,
) -> tuple[int, int]:
    """Return ``(x, y)`` for a child window relative to its parent.

    Pure geometry: no Tk calls. ``work_area`` is ``(left, top, right, bottom)``.
    ``occupied`` lists other mapped windows as ``(x, y, w, h)`` to avoid.
    """
    if child_w < 1:
        child_w = 1
    if child_h < 1:
        child_h = 1

    if mode == "center":
        x = parent_x + (parent_w - child_w) // 2
        y = parent_y + (parent_h - child_h) // 2
        return _clamp_to_work_area(x, y, child_w, child_h, work_area)

    # Preferred outside-edge slots (fully on the same monitor).
    outside: list[tuple[int, int]] = [
        (parent_x + parent_w + edge_gap, parent_y),  # right
        (parent_x - child_w - edge_gap, parent_y),  # left
        (parent_x, parent_y + parent_h + edge_gap),  # below
        (parent_x, parent_y - child_h - edge_gap),  # above
    ]

    # Inner corners of the parent (always same monitor if parent is).
    inside: list[tuple[int, int]] = [
        (parent_x + parent_w - child_w - inner_inset, parent_y + inner_inset),  # top-right
        (parent_x + inner_inset, parent_y + inner_inset),  # top-left
        (
            parent_x + parent_w - child_w - inner_inset,
            parent_y + parent_h - child_h - inner_inset,
        ),  # bottom-right
        (parent_x + inner_inset, parent_y + parent_h - child_h - inner_inset),  # bottom-left
    ]

    candidates: list[tuple[int, int]] = []
    for x, y in outside:
        if _fits_in_work_area(x, y, child_w, child_h, work_area):
            candidates.append((x, y))
    for x, y in inside:
        cx, cy = _clamp_to_work_area(x, y, child_w, child_h, work_area)
        candidates.append((cx, cy))

    # De-duplicate while preserving order.
    seen: set[tuple[int, int]] = set()
    unique: list[tuple[int, int]] = []
    for pos in candidates:
        if pos not in seen:
            seen.add(pos)
            unique.append(pos)

    best = unique[0]
    best_overlap = _total_overlap(best[0], best[1], child_w, child_h, occupied)
    for pos in unique[1:]:
        overlap = _total_overlap(pos[0], pos[1], child_w, child_h, occupied)
        if overlap < best_overlap:
            best = pos
            best_overlap = overlap
            if best_overlap == 0:
                break

    if best_overlap == 0:
        return best

    # Cascade from the best slot until overlap drops or we run out of steps.
    bx, by = best
    for step in range(1, 12):
        cx = bx + step * cascade_step
        cy = by + step * cascade_step
        cx, cy = _clamp_to_work_area(cx, cy, child_w, child_h, work_area)
        overlap = _total_overlap(cx, cy, child_w, child_h, occupied)
        if overlap < best_overlap:
            return cx, cy
    return best


def _work_area_win32(hwnd: int) -> tuple[int, int, int, int] | None:
    """Return the work area of the monitor containing ``hwnd``, or None."""
    if sys.platform != "win32":
        return None
    try:
        import ctypes
        from ctypes import wintypes

        user32 = ctypes.windll.user32

        class RECT(ctypes.Structure):
            _fields_ = [
                ("left", wintypes.LONG),
                ("top", wintypes.LONG),
                ("right", wintypes.LONG),
                ("bottom", wintypes.LONG),
            ]

        class MONITORINFO(ctypes.Structure):
            _fields_ = [
                ("cbSize", wintypes.DWORD),
                ("rcMonitor", RECT),
                ("rcWork", RECT),
                ("dwFlags", wintypes.DWORD),
            ]

        MONITOR_DEFAULTTONEAREST = 2
        hmon = user32.MonitorFromWindow(wintypes.HWND(hwnd), MONITOR_DEFAULTTONEAREST)
        if not hmon:
            return None
        info = MONITORINFO()
        info.cbSize = ctypes.sizeof(MONITORINFO)
        if not user32.GetMonitorInfoW(hmon, ctypes.byref(info)):
            return None
        r = info.rcWork
        return int(r.left), int(r.top), int(r.right), int(r.bottom)
    except Exception:
        return None


def _work_area_from_widget(widget) -> tuple[int, int, int, int]:
    """Work area for the monitor containing ``widget`` (fallback: full screen)."""
    try:
        hwnd = int(widget.winfo_id())
        area = _work_area_win32(hwnd)
        if area is not None:
            return area
    except Exception:
        pass

    try:
        # Virtual root covers the multi-monitor desktop on some platforms.
        left = int(widget.winfo_vrootx())
        top = int(widget.winfo_vrooty())
        width = int(widget.winfo_vrootwidth())
        height = int(widget.winfo_vrootheight())
        if width > 0 and height > 0:
            return left, top, left + width, top + height
    except Exception:
        pass

    try:
        sw = int(widget.winfo_screenwidth())
        sh = int(widget.winfo_screenheight())
        return 0, 0, sw, sh
    except Exception:
        return 0, 0, 1920, 1080


def _mapped_toplevel_rects(root, exclude) -> list[tuple[int, int, int, int]]:
    """Collect ``(x, y, w, h)`` for mapped Toplevels under ``root``, skipping ``exclude``."""
    occupied: list[tuple[int, int, int, int]] = []
    try:
        children: Iterable = root.winfo_children()
    except Exception:
        return occupied
    for child in children:
        if child is exclude:
            continue
        try:
            if not child.winfo_exists():
                continue
            # Only consider Toplevel-like windows (not Frames).
            if child.winfo_toplevel() is not child:
                continue
            if str(child.state()) != "normal":
                continue
            x = int(child.winfo_rootx())
            y = int(child.winfo_rooty())
            w = max(int(child.winfo_width()), 1)
            h = max(int(child.winfo_height()), 1)
            occupied.append((x, y, w, h))
        except Exception:
            continue
    return occupied


def place_toplevel(
    win,
    parent,
    *,
    mode: PlacementMode = "edge",
) -> tuple[int, int]:
    """Position ``win`` relative to ``parent`` and return the chosen ``(x, y)``.

    Call after the dialog's widgets are packed and before ``grab_set()`` /
    ``wait_window()``. Does not reposition windows that the user has already
    moved (callers should only invoke this on first open, not on ``lift()``).
    """
    try:
        win.update_idletasks()
        parent.update_idletasks()
    except Exception:
        pass

    try:
        parent_x = int(parent.winfo_rootx())
        parent_y = int(parent.winfo_rooty())
        parent_w = max(int(parent.winfo_width()), 1)
        parent_h = max(int(parent.winfo_height()), 1)
    except Exception:
        parent_x, parent_y, parent_w, parent_h = 100, 100, 800, 600

    try:
        child_w = max(int(win.winfo_reqwidth()), int(win.winfo_width()), 1)
        child_h = max(int(win.winfo_reqheight()), int(win.winfo_height()), 1)
    except Exception:
        child_w, child_h = 400, 300

    work_area = _work_area_from_widget(parent)

    root = parent.winfo_toplevel() if hasattr(parent, "winfo_toplevel") else parent
    occupied = _mapped_toplevel_rects(root, exclude=win) if mode == "edge" else []

    x, y = choose_window_position(
        parent_x=parent_x,
        parent_y=parent_y,
        parent_w=parent_w,
        parent_h=parent_h,
        child_w=child_w,
        child_h=child_h,
        work_area=work_area,
        mode=mode,
        occupied=occupied,
    )
    try:
        win.geometry(f"+{x}+{y}")
    except Exception:
        pass
    return x, y
