"""Unit tests for spectrum-delete reference helpers (Tk-free)."""

from __future__ import annotations

from myutils import library_storage


def _spec(*refs: str) -> dict:
    blocks = [
        library_storage.make_interval_block(ref, 1e10 if ref != "Dark" else 0.0, 10.0)
        for ref in refs
    ]
    return {
        "version": library_storage.STIMULUS_SCHEMA_VERSION,
        "name": "demo",
        "total_duration": 10.0 * len(blocks),
        "blocks": blocks,
    }


def test_spec_references_spectrum_hit() -> None:
    spec = _spec("440 nm", "custom_A", "Dark")
    assert library_storage.spec_references_spectrum(spec, "custom_A") is True


def test_spec_references_spectrum_miss() -> None:
    spec = _spec("440 nm", "Dark")
    assert library_storage.spec_references_spectrum(spec, "custom_A") is False


def test_spec_references_spectrum_multi_block() -> None:
    spec = _spec("custom_A", "560 nm", "custom_A")
    assert library_storage.spec_references_spectrum(spec, "custom_A") is True
    assert library_storage.spec_references_spectrum(spec, "560 nm") is True
    assert library_storage.spec_references_spectrum(spec, "Xenon") is False


def test_replace_spectrum_ref_rewrites_matches() -> None:
    spec = _spec("custom_A", "Dark", "custom_A")
    original_ref = spec["blocks"][0]["spectrum_ref"]
    out = library_storage.replace_spectrum_ref(spec, "custom_A", "560 nm")
    assert out["blocks"][0]["spectrum_ref"] == "560 nm"
    assert out["blocks"][1]["spectrum_ref"] == "Dark"
    assert out["blocks"][2]["spectrum_ref"] == "560 nm"
    assert out["blocks"][0]["intensity"] == 1e10
    # Original unchanged
    assert spec["blocks"][0]["spectrum_ref"] == original_ref
    assert spec["blocks"][0]["spectrum_ref"] == "custom_A"


def test_replace_spectrum_ref_dark_zeros_intensity() -> None:
    spec = _spec("custom_A", "560 nm")
    out = library_storage.replace_spectrum_ref(spec, "custom_A", "Dark")
    assert out["blocks"][0]["spectrum_ref"] == "Dark"
    assert out["blocks"][0]["intensity"] == 0.0
    assert out["blocks"][1]["spectrum_ref"] == "560 nm"
    assert out["blocks"][1]["intensity"] == 1e10
    assert spec["blocks"][0]["intensity"] == 1e10


def test_replace_spectrum_ref_preserves_metadata() -> None:
    spec = _spec("custom_A")
    out = library_storage.replace_spectrum_ref(spec, "custom_A", "Xenon")
    assert out["name"] == spec["name"]
    assert out["version"] == spec["version"]
    assert out["total_duration"] == spec["total_duration"]
    assert len(out["blocks"]) == len(spec["blocks"])
