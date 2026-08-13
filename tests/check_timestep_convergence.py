"""Compare integrated activity at the default timestep vs half that timestep.

Runs all six manuscript intensity-response protocols at the default model
``rate`` (0.01 s) and at ``rate / 2``, then reports the percent change in the
trapezoid integral of ``currentGlobalGain``.

From the repository root::

    python tests/check_timestep_convergence.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from myutils.init import loadInSpectra
from myutils.model_functions import (
    get_predict_melanopsin_defaults,
    predictMelanopsin,
    setupModelRun,
)

MANUSCRIPT_LABELS = (
    "440-440",
    "560-560",
    "440-560",
    "560-440",
    "Xenon",
    "Xenon with ocular filtering",
)


def _spectra_env() -> dict:
    (
        _spec440,
        _spec560,
        intensity440,
        intensity560,
        wlenshared,
        xenon,
        xenon_eye,
    ) = loadInSpectra()
    return {
        "wlenshared": wlenshared,
        "intensity440": intensity440,
        "intensity560": intensity560,
        "xenon": xenon,
        "xenon_eye": xenon_eye,
    }


def _build_stimulus(label: str, env: dict) -> dict:
    """Match GUI ``_build_stimulus`` / tutorial manuscript cells."""
    wlenshared = env["wlenshared"]
    intensity440 = env["intensity440"]
    intensity560 = env["intensity560"]
    xph = np.asarray(env["xenon"]["photons"], dtype=float)
    xeph = np.asarray(env["xenon_eye"]["photons"], dtype=float)

    if label == "440-440":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 440",
            specwlen=wlenshared,
            specint=intensity440,
        )
    if label == "560-560":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 560",
            specwlen=wlenshared,
            specint=intensity560,
        )
    if label == "440-560":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 440 560",
            specwlen=[wlenshared, wlenshared],
            specint=[intensity440, intensity560],
        )
    if label == "560-440":
        return setupModelRun(
            "intensity response",
            "spectrum matched to 560 440",
            specwlen=[wlenshared, wlenshared],
            specint=[intensity560, intensity440],
        )
    if label == "Xenon":
        return setupModelRun(
            "intensity response",
            "xenon",
            specwlen=wlenshared,
            specint=xph,
        )
    if label == "Xenon with ocular filtering":
        return setupModelRun(
            "intensity response",
            "xenon eye",
            specwlen=wlenshared,
            specint=xeph,
        )
    raise ValueError(f"Unknown stimulus label: {label!r}")


def _integrated_activity(model: dict) -> float:
    return float(np.trapezoid(model["currentGlobalGain"], model["xaxis"]))


def main() -> None:
    default_rate = float(get_predict_melanopsin_defaults()["rate"])
    half_rate = default_rate / 2.0
    env = _spectra_env()

    print(
        f"Default rate = {default_rate:g} s; halved rate = {half_rate:g} s\n"
        "Integrated activity = trapezoid of currentGlobalGain over time.\n"
    )
    header = (
        f"{'stimulus':<32} {'I(default)':>14} {'I(half)':>14} {'% change':>12}"
    )
    print(header)
    print("-" * len(header))

    rows: list[tuple[str, float, float, float]] = []
    for label in MANUSCRIPT_LABELS:
        stim = _build_stimulus(label, env)
        print(f"\n=== {label}: rate={default_rate:g} ===")
        model_default = predictMelanopsin(stim, config={"rate": default_rate})
        print(f"=== {label}: rate={half_rate:g} ===")
        model_half = predictMelanopsin(stim, config={"rate": half_rate})

        i_default = _integrated_activity(model_default)
        i_half = _integrated_activity(model_half)
        pct = 100.0 * (i_half - i_default) / i_default
        rows.append((label, i_default, i_half, pct))

    print("\n")
    print(header)
    print("-" * len(header))
    for label, i_default, i_half, pct in rows:
        print(f"{label:<32} {i_default:14.6g} {i_half:14.6g} {pct:11.4f}%")


if __name__ == "__main__":
    main()
