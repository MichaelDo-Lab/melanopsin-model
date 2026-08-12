"""Unit tests for Data Comparator dataset import helpers (Tk-free)."""

from __future__ import annotations

import numpy as np
import pytest

from myutils.melanopsin_gui import (
    _parse_time_unit_scale,
    _read_data_trace,
)


def test_parse_time_unit_scale_seconds_default() -> None:
    scale, label = _parse_time_unit_scale("time")
    assert scale == 1.0
    assert label == "s"


def test_parse_time_unit_scale_ms_in_parens() -> None:
    scale, label = _parse_time_unit_scale("Time (ms)")
    assert scale == pytest.approx(1e-3)
    assert label == "ms"


def test_parse_time_unit_scale_underscore_suffix() -> None:
    scale, label = _parse_time_unit_scale("time_ms")
    assert scale == pytest.approx(1e-3)
    assert label == "ms"


def test_parse_time_unit_scale_minutes() -> None:
    scale, label = _parse_time_unit_scale("elapsed [min]")
    assert scale == pytest.approx(60.0)
    assert label == "min"


def test_headerless_single_column(tmp_path) -> None:
    path = tmp_path / "trace.csv"
    path.write_text("1.0\n2.0\n3.0\n", encoding="utf-8")
    loaded = _read_data_trace(str(path))
    assert loaded["headerless"] is True
    assert loaded["time_s"] is None
    assert loaded["time_column"] is None
    assert list(loaded["columns"].keys()) == ["value"]
    np.testing.assert_allclose(loaded["columns"]["value"], [1.0, 2.0, 3.0])


def test_time_current_seconds(tmp_path) -> None:
    path = tmp_path / "with_time.csv"
    path.write_text("time,current\n0.0,1.0\n0.1,2.0\n0.2,3.0\n", encoding="utf-8")
    loaded = _read_data_trace(str(path))
    assert loaded["headerless"] is False
    assert loaded["time_column"] == "time"
    assert loaded["time_unit_label"] == "s"
    assert loaded["time_offset_s"] == pytest.approx(0.0)
    assert list(loaded["columns"].keys()) == ["current"]
    np.testing.assert_allclose(loaded["time_s"], [0.0, 0.1, 0.2])
    np.testing.assert_allclose(loaded["columns"]["current"], [1.0, 2.0, 3.0])


def test_time_ms_scaled_and_zero_based(tmp_path) -> None:
    path = tmp_path / "time_ms.csv"
    path.write_text(
        "Time (ms),signal\n1000,10\n2000,20\n3000,30\n",
        encoding="utf-8",
    )
    loaded = _read_data_trace(str(path))
    assert loaded["time_column"] == "Time (ms)"
    assert loaded["time_unit_label"] == "ms"
    assert loaded["time_offset_s"] == pytest.approx(1.0)
    # Zero-based: original 1 s, 2 s, 3 s -> 0, 1, 2
    np.testing.assert_allclose(loaded["time_s"], [0.0, 1.0, 2.0])
    np.testing.assert_allclose(loaded["columns"]["signal"], [10.0, 20.0, 30.0])


def test_multicolumn_headered_no_time(tmp_path) -> None:
    path = tmp_path / "cells.csv"
    path.write_text(
        "Intensity,cell_a,cell_b\n1.0,0.1,0.2\n2.0,0.3,0.4\n",
        encoding="utf-8",
    )
    loaded = _read_data_trace(str(path))
    assert loaded["time_s"] is None
    assert loaded["headerless"] is False
    assert list(loaded["columns"].keys()) == ["Intensity", "cell_a", "cell_b"]


def test_nan_rows_retained_until_selection(tmp_path) -> None:
    """NaNs stay in arrays so time/value stay aligned; filtering is at selection."""
    path = tmp_path / "with_nan.csv"
    path.write_text(
        "time,current\n0.0,1.0\n0.1,nan\n0.2,3.0\n",
        encoding="utf-8",
    )
    loaded = _read_data_trace(str(path))
    assert loaded["time_s"] is not None
    assert loaded["time_s"].shape == loaded["columns"]["current"].shape
    assert loaded["time_s"].size == 3
    assert np.isnan(loaded["columns"]["current"][1])
    # Joint finite mask (what _selected_dataset_trace applies)
    mask = np.isfinite(loaded["time_s"]) & np.isfinite(loaded["columns"]["current"])
    t = loaded["time_s"][mask]
    v = loaded["columns"]["current"][mask]
    np.testing.assert_allclose(t, [0.0, 0.2])
    np.testing.assert_allclose(v, [1.0, 3.0])
