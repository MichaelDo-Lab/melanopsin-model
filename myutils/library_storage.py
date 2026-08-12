"""User-library storage backend for spectra, stimuli, and configurations.

Stores are rooted at ``data/user_library/``:

* ``spectra/``  - one ``<safe_name>.npz`` per user-facing spectrum name plus a
  ``registry.json`` mapping names to file/hash metadata. Older installs may
  still have content-hashed ``<sha256>.npz`` filenames; those keep working.
* ``stimuli/``  - one tiny ``<stim_id>.json`` per stimulus, written in the v2
  block-based schema. ``index.json`` maps user names to file ids.
* ``configurations/`` - user-saved model-parameter presets (``.json``).

Manuscript-shipped spectra/stimuli are baked into the model code and not
managed here. On first launch after the upgrade, legacy paths
(``data/custom_spectra/``, ``data/longform_uploads/``, ``configs/config.json``)
are removed entirely and the new tree is created from scratch. The reset is
guarded by ``data/user_library/.reset_done`` so it only runs once.
"""

from __future__ import annotations

import copy
import hashlib
import json
import numbers
import re
import shutil
import uuid
from pathlib import Path
from typing import Iterable

import numpy as np

from .app_paths import configs_dir, data_dir

# ---------------------------------------------------------------------------
# Layout / constants
# ---------------------------------------------------------------------------

# Block-based schema introduced as part of the storage overhaul.
STIMULUS_SCHEMA_VERSION = 2

# Reserved names that map onto code-level (manuscript) spectra and are never
# stored by this module. ``Dark`` is a synthetic zero-intensity placeholder.
BUILTIN_SPECTRUM_NAMES = ("440 nm", "560 nm", "Xenon", "Xenon (eye)", "Dark")

_LEGACY_CUSTOM_SPECTRA_DIR = data_dir() / "custom_spectra"
_LEGACY_LONGFORM_DIR = data_dir() / "longform_uploads"
_LEGACY_STIMULUS_LIBRARY = configs_dir() / "config.json"


def user_library_root() -> Path:
    return data_dir() / "user_library"


def spectra_dir() -> Path:
    return user_library_root() / "spectra"


def stimuli_dir() -> Path:
    return user_library_root() / "stimuli"


def configurations_dir() -> Path:
    """Path to user-saved model-parameter presets under ``user_library``."""
    return user_library_root() / "configurations"


def _spectra_registry_path() -> Path:
    return spectra_dir() / "registry.json"


def _stimuli_index_path() -> Path:
    return stimuli_dir() / "index.json"


def _reset_marker_path() -> Path:
    return user_library_root() / ".reset_done"


# ---------------------------------------------------------------------------
# Reset / migration
# ---------------------------------------------------------------------------

def reset_legacy_paths_if_needed() -> bool:
    """Remove pre-overhaul storage and create the new ``user_library/`` tree.

    Runs at most once per install (guarded by ``.reset_done``). Returns True
    when a reset actually happened, False if the tree already exists.
    """
    if _reset_marker_path().exists():
        ensure_library_tree()
        return False

    for path in (_LEGACY_CUSTOM_SPECTRA_DIR, _LEGACY_LONGFORM_DIR):
        if path.exists():
            shutil.rmtree(path, ignore_errors=True)
    if _LEGACY_STIMULUS_LIBRARY.exists():
        try:
            _LEGACY_STIMULUS_LIBRARY.unlink()
        except OSError:
            pass
    # Also drop the .bak file left by the old in-place embed migration.
    legacy_bak = _LEGACY_STIMULUS_LIBRARY.with_suffix(".json.bak")
    if legacy_bak.exists():
        try:
            legacy_bak.unlink()
        except OSError:
            pass

    ensure_library_tree()
    try:
        _reset_marker_path().write_text("1", encoding="utf-8")
    except OSError:
        pass
    return True


def ensure_library_tree() -> None:
    """Create the empty ``user_library/`` skeleton if it does not yet exist."""
    for d in (spectra_dir(), stimuli_dir(), configurations_dir()):
        d.mkdir(parents=True, exist_ok=True)
    if not _spectra_registry_path().exists():
        _atomic_write_json(_spectra_registry_path(), {"version": 1, "spectra": {}})
    if not _stimuli_index_path().exists():
        _atomic_write_json(_stimuli_index_path(), {"version": 1, "stimuli": {}})


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def _atomic_write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    tmp.replace(path)


def _read_json(path: Path, default):
    if not path.exists():
        return default
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return default


def _safe_slug(name: str, *, fallback: str = "item") -> str:
    slug = re.sub(r"[^\w\-]+", "_", name).strip("_")
    return (slug[:60] or fallback)


def _new_id(prefix: str, used: Iterable[str]) -> str:
    used_set = set(used)
    while True:
        candidate = f"{prefix}_{uuid.uuid4().hex[:10]}"
        if candidate not in used_set:
            return candidate


# ---------------------------------------------------------------------------
# Spectrum registry
# ---------------------------------------------------------------------------

def _hash_spectrum(wlen: np.ndarray, intensity: np.ndarray) -> str:
    """Return a stable content hash for ``(wlen, intensity)``."""
    wlen_arr = np.ascontiguousarray(np.asarray(wlen, dtype=np.float64))
    int_arr = np.ascontiguousarray(np.asarray(intensity, dtype=np.float64))
    h = hashlib.sha256()
    h.update(b"wlen:")
    h.update(wlen_arr.tobytes())
    h.update(b"|int:")
    h.update(int_arr.tobytes())
    return h.hexdigest()[:32]


def _spectrum_entry_filename(entry: dict | None) -> str:
    if not isinstance(entry, dict):
        return ""
    fname = entry.get("file")
    if isinstance(fname, str) and fname:
        return fname
    content_hash = entry.get("hash")
    if isinstance(content_hash, str) and content_hash:
        return f"{content_hash}.npz"
    return ""


def _spectrum_filename_for_name(
    name: str, spectra: dict, content_hash: str
) -> str:
    """Choose an on-disk NPZ name derived from the user-facing spectrum name."""
    existing = spectra.get(name)
    if isinstance(existing, dict):
        current = _spectrum_entry_filename(existing)
        if current:
            return current

    base = _safe_slug(name, fallback="spectrum")
    fname = f"{base}.npz"
    used = {
        _spectrum_entry_filename(entry)
        for other, entry in spectra.items()
        if other != name and isinstance(entry, dict)
    }
    used.discard("")
    if fname not in used:
        return fname
    return f"{base}_{content_hash[:8]}.npz"


def _load_spectra_registry() -> dict:
    ensure_library_tree()
    data = _read_json(_spectra_registry_path(), {"version": 1, "spectra": {}})
    if not isinstance(data, dict) or "spectra" not in data:
        data = {"version": 1, "spectra": {}}
    if not isinstance(data["spectra"], dict):
        data["spectra"] = {}
    return data


def _write_spectra_registry(registry: dict) -> None:
    _atomic_write_json(_spectra_registry_path(), registry)


def list_spectra_names() -> list[str]:
    return sorted(_load_spectra_registry().get("spectra", {}).keys())


def has_spectrum(name: str) -> bool:
    return name in _load_spectra_registry().get("spectra", {})


def save_spectrum(name: str, wlen: np.ndarray, intensity: np.ndarray) -> str:
    """Persist a spectrum under ``name`` as ``data/user_library/spectra/<name>.npz``.

    The on-disk filename is a filesystem-safe form of the user-facing name.
    Returns the user-facing name (unchanged). Raises ``ValueError`` if the
    name is already taken by a different content hash.
    """
    name = str(name).strip()
    if not name:
        raise ValueError("Spectrum name cannot be empty.")
    if name in BUILTIN_SPECTRUM_NAMES:
        raise ValueError(f"{name!r} is reserved for a manuscript spectrum.")

    wlen_arr = np.asarray(wlen, dtype=float)
    int_arr = np.asarray(intensity, dtype=float)
    if wlen_arr.ndim != 1 or int_arr.ndim != 1 or wlen_arr.shape != int_arr.shape:
        raise ValueError("wlen and intensity must be 1-D arrays of the same length.")

    content_hash = _hash_spectrum(wlen_arr, int_arr)
    registry = _load_spectra_registry()
    spectra = registry["spectra"]

    existing = spectra.get(name)
    if existing is not None and existing.get("hash") != content_hash:
        raise ValueError(
            f"A different spectrum is already saved under {name!r}; choose a different name."
        )

    fname = _spectrum_filename_for_name(name, spectra, content_hash)
    fpath = spectra_dir() / fname
    np.savez_compressed(str(fpath), wlen=wlen_arr, intensity=int_arr)

    spectra[name] = {"hash": content_hash, "file": fname}
    _write_spectra_registry(registry)
    return name


def load_spectrum(name: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Return ``(wlen, intensity)`` for a custom spectrum, or None if absent."""
    registry = _load_spectra_registry()
    entry = registry.get("spectra", {}).get(name)
    if not isinstance(entry, dict):
        return None
    fname = _spectrum_entry_filename(entry)
    if not fname:
        return None
    fpath = spectra_dir() / fname
    if not fpath.exists():
        return None
    try:
        data = np.load(str(fpath))
        return (
            np.asarray(data["wlen"], dtype=float),
            np.asarray(data["intensity"], dtype=float),
        )
    except Exception:
        return None


def load_all_spectra() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Eagerly load every saved spectrum into memory."""
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name in list_spectra_names():
        data = load_spectrum(name)
        if data is not None:
            out[name] = data
    return out


def delete_spectrum(name: str) -> bool:
    """Drop the registry entry for ``name`` and remove its NPZ if orphaned."""
    registry = _load_spectra_registry()
    spectra = registry.get("spectra", {})
    entry = spectra.pop(name, None)
    if entry is None:
        return False
    _write_spectra_registry(registry)
    fname = _spectrum_entry_filename(entry)
    if not fname:
        return True
    # Remove the NPZ only when nothing else still points at this file.
    if not any(_spectrum_entry_filename(e) == fname for e in spectra.values()):
        fpath = spectra_dir() / fname
        if fpath.exists():
            try:
                fpath.unlink()
            except OSError:
                pass
    return True


# ---------------------------------------------------------------------------
# Stimulus store (v2 block-based)
# ---------------------------------------------------------------------------

def _load_stimuli_index() -> dict:
    ensure_library_tree()
    data = _read_json(_stimuli_index_path(), {"version": 1, "stimuli": {}})
    if not isinstance(data, dict) or not isinstance(data.get("stimuli"), dict):
        data = {"version": 1, "stimuli": {}}
    return data


def _write_stimuli_index(index: dict) -> None:
    _atomic_write_json(_stimuli_index_path(), index)


def list_stimulus_names() -> list[str]:
    return sorted(_load_stimuli_index().get("stimuli", {}).keys())


def stimulus_summary(name: str) -> dict | None:
    index = _load_stimuli_index()
    entry = index.get("stimuli", {}).get(name)
    if not isinstance(entry, dict):
        return None
    return {
        "name": name,
        "file": entry.get("file"),
        "total_duration": float(entry.get("total_duration", 0.0)),
        "n_blocks": int(entry.get("n_blocks", 0)),
        "kind": entry.get("kind", "manual"),
    }


def load_stimulus(name: str) -> dict | None:
    index = _load_stimuli_index()
    entry = index.get("stimuli", {}).get(name)
    if not isinstance(entry, dict):
        return None
    fname = entry.get("file")
    if not fname:
        return None
    spec = _read_json(stimuli_dir() / fname, None)
    if spec is None:
        return None
    try:
        return validate_stimulus_spec(spec)
    except ValueError:
        return None


def load_all_stimuli() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for name in list_stimulus_names():
        spec = load_stimulus(name)
        if spec is not None:
            out[name] = spec
    return out


def save_stimulus(spec: dict) -> str:
    """Validate ``spec`` and persist it; returns the canonical name."""
    spec = validate_stimulus_spec(spec)
    name = spec["name"]
    index = _load_stimuli_index()
    existing = index["stimuli"].get(name)
    fname = (existing or {}).get("file") or f"stim_{uuid.uuid4().hex[:10]}_{_safe_slug(name)}.json"
    _atomic_write_json(stimuli_dir() / fname, spec)
    index["stimuli"][name] = {
        "file": fname,
        "total_duration": float(spec["total_duration"]),
        "n_blocks": len(spec["blocks"]),
        "kind": "manual",
    }
    _write_stimuli_index(index)
    return name


def delete_stimulus(name: str) -> bool:
    index = _load_stimuli_index()
    entry = index.get("stimuli", {}).pop(name, None)
    if entry is None:
        return False
    _write_stimuli_index(index)
    fname = entry.get("file")
    if fname:
        fpath = stimuli_dir() / fname
        if fpath.exists():
            try:
                fpath.unlink()
            except OSError:
                pass
    return True


# ---------------------------------------------------------------------------
# v2 block-based spec validation / expansion
# ---------------------------------------------------------------------------

def _is_positive_number(value) -> bool:
    return isinstance(value, numbers.Real) and float(value) > 0 and np.isfinite(float(value))


def _is_nonneg_number(value) -> bool:
    return isinstance(value, numbers.Real) and float(value) >= 0 and np.isfinite(float(value))


def validate_stimulus_spec(data: object) -> dict:
    """Return a normalized v2 stimulus spec. Raises ``ValueError`` on problems.

    The v2 format describes a stimulus as an ordered list of ``blocks``. Each
    block is:

    - ``{"type": "interval", "spectrum_ref": str, "intensity": float, "duration": float}``
      where ``spectrum_ref`` is a built-in name or a custom spectrum name.
    """
    if not isinstance(data, dict):
        raise ValueError("Stimulus spec must be a JSON object.")
    version = data.get("version")
    if version != STIMULUS_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported stimulus schema version {version!r} "
            f"(expected {STIMULUS_SCHEMA_VERSION})."
        )
    name = data.get("name")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("Stimulus spec is missing a non-empty 'name'.")

    blocks_raw = data.get("blocks")
    if not isinstance(blocks_raw, list) or not blocks_raw:
        raise ValueError("Stimulus spec 'blocks' must be a non-empty list.")

    blocks: list[dict] = []
    total = 0.0
    for i, raw in enumerate(blocks_raw):
        if not isinstance(raw, dict):
            raise ValueError(f"Block {i + 1}: must be an object.")
        kind = raw.get("type")
        if kind == "interval":
            spectrum_ref = raw.get("spectrum_ref")
            if not isinstance(spectrum_ref, str) or not spectrum_ref.strip():
                raise ValueError(f"Block {i + 1}: missing 'spectrum_ref'.")
            intensity = raw.get("intensity", 0.0)
            if not _is_nonneg_number(intensity):
                raise ValueError(
                    f"Block {i + 1}: 'intensity' must be a non-negative number."
                )
            duration = raw.get("duration")
            if not _is_positive_number(duration):
                raise ValueError(
                    f"Block {i + 1}: 'duration' must be a positive number."
                )
            normalized = {
                "type": "interval",
                "spectrum_ref": spectrum_ref.strip(),
                "intensity": float(intensity),
                "duration": float(duration),
            }
            total += float(duration)
            blocks.append(normalized)
        else:
            raise ValueError(f"Block {i + 1}: unknown type {kind!r}.")

    declared_total = data.get("total_duration")
    if declared_total is None:
        total_duration = total
    else:
        if not _is_positive_number(declared_total):
            raise ValueError("'total_duration' must be a positive number.")
        if total - float(declared_total) > 1e-6 * max(1.0, total):
            raise ValueError(
                "'total_duration' is shorter than the sum of block durations "
                f"({total:g} s)."
            )
        total_duration = float(declared_total)

    return {
        "version": STIMULUS_SCHEMA_VERSION,
        "name": name.strip(),
        "total_duration": float(total_duration),
        "blocks": blocks,
    }

def n_intervals_in_block(block: dict) -> int:
    """Return the number of expanded intervals represented by ``block``."""
    if block.get("type") == "interval":
        return 1
    return 0


def block_total_duration(block: dict) -> float:
    if block.get("type") == "interval":
        return float(block["duration"])
    return 0.0


def total_intervals(spec: dict) -> int:
    return sum(n_intervals_in_block(b) for b in spec.get("blocks", []))


def iter_block_durations(block: dict):
    """Yield per-interval ``duration`` values for ``block`` (no spectra slice)."""
    if block.get("type") == "interval":
        yield float(block["duration"])


def iter_block_intervals(block: dict, upload_loader=None):
    """Yield ``(spectrum_ref, intensity, wlen, intensity_arr, duration)`` per interval.

    For ``interval`` blocks, ``wlen`` and ``intensity_arr`` are None (the caller
    resolves the spectrum by ref). ``upload_loader`` is unused and kept only for
    call-site compatibility.
    """
    kind = block.get("type")
    if kind == "interval":
        yield (
            block["spectrum_ref"],
            float(block["intensity"]),
            None,
            None,
            float(block["duration"]),
        )
        return
    raise ValueError(f"Unknown block type {kind!r}.")

# ---------------------------------------------------------------------------
# Block helpers used by the GUI
# ---------------------------------------------------------------------------

def make_interval_block(spectrum_ref: str, intensity: float, duration: float) -> dict:
    return {
        "type": "interval",
        "spectrum_ref": str(spectrum_ref).strip(),
        "intensity": float(intensity),
        "duration": float(duration),
    }


def spec_references_spectrum(spec: dict, name: str) -> bool:
    """Return True if any interval block in ``spec`` references ``name``."""
    target = str(name).strip()
    for block in spec.get("blocks", []):
        if not isinstance(block, dict):
            continue
        if block.get("type") != "interval":
            continue
        if str(block.get("spectrum_ref", "")).strip() == target:
            return True
    return False


def replace_spectrum_ref(spec: dict, old: str, new: str) -> dict:
    """Return a deep-copied spec with matching ``spectrum_ref`` values rewritten.

    When ``new`` is ``"Dark"``, matching blocks also have their intensity forced
    to ``0.0`` (matching the Stimulus Builder's Dark lock). Other fields
    (``name``, ``version``, ``total_duration``, block order) are left untouched.
    """
    old_name = str(old).strip()
    new_name = str(new).strip()
    out = copy.deepcopy(spec)
    for block in out.get("blocks", []):
        if not isinstance(block, dict):
            continue
        if block.get("type") != "interval":
            continue
        if str(block.get("spectrum_ref", "")).strip() != old_name:
            continue
        block["spectrum_ref"] = new_name
        if new_name == "Dark":
            block["intensity"] = 0.0
    return out


def coerce_blocks_from_legacy(intervals: list[dict]) -> list[dict]:
    """Convert old v1 ``intervals`` (flat list) to v2 interval blocks.

    Kept for callers that still build a flat list and want to round-trip
    through the new validator without restructuring (e.g. the Stimulus Builder's
    manual-grid save path).
    """
    blocks: list[dict] = []
    for iv in intervals:
        blocks.append(
            make_interval_block(
                spectrum_ref=str(iv.get("spectrum", iv.get("spectrum_ref", ""))).strip(),
                intensity=float(iv.get("intensity", 0.0)),
                duration=float(iv.get("duration", 0.0)),
            )
        )
    return blocks
