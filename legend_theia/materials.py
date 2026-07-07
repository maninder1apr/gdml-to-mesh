"""
legend_theia/materials.py

Builds Theia Medium objects from metadata/materials.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from theia.material import Medium
from theia.lookup import Table
from theia.model import RayleighScatteringPhaseFunction
from theia.property import FloatProperty, TableProperty
from theia.volume import Attenuating, Transparent, Fluorescent
import theia.units as u

from .registry import canonical

# ─────────────────────────────────────────────────────────────────────────────
# Global wavelength grid — locked per base design §4
# ─────────────────────────────────────────────────────────────────────────────

WAVELENGTH_MIN: float = 110.0   # nm
WAVELENGTH_MAX: float = 700.0   # nm
WAVELENGTH_N:   int   = 512


def _grid() -> np.ndarray:
    return np.linspace(WAVELENGTH_MIN, WAVELENGTH_MAX, WAVELENGTH_N)


def _interp(wavelengths: list[float], values: list[float]) -> np.ndarray:
    """Interpolate onto the global wavelength grid."""
    return np.interp(_grid(), wavelengths, values)


def _mfp_to_coef(mfp_mm: np.ndarray) -> np.ndarray:
    """Convert mean free path in mm to absorption/scattering coefficient in 1/m."""
    return 1000.0 / np.clip(mfp_mm, 1e-6, None)

def _compute_cdf(wavelengths, values, num_values = 1024):
    wl_support = np.linspace(wavelengths[0], wavelengths[-1], num_values)
    values_interp = np.interp(wl_support, wavelengths, values)

    cdf = cumulative_trapezoid(values_interp, wl_support, initial=0)
    # normalize
    cdf /= cdf[-1]

    x_values = np.linspace(0, 1, num_values)
    return np.interp(x_values, wl_support, cdf)

def _compute_ppf(wavelengths, values, num_values = 1024):
    wl_support = np.linspace(wavelengths[0], wavelengths[-1], num_values)
    values_interp = np.interp(wl_support, wavelengths, values)

    cdf = cumulative_trapezoid(values_interp, wl_support, initial=0)
    # normalize
    cdf /= cdf[-1]

    x_values = np.linspace(0, 1, num_values)
    # inverse cdf -> swap x-axis and y-axis
    return np.interp(x_values, cdf, wl_support)


def build_media(
    materials_json: str | Path,
    *,
    wavelength_min: float = WAVELENGTH_MIN,
    wavelength_max: float = WAVELENGTH_MAX,
    n_samples: int = WAVELENGTH_N,
) -> dict[str, Medium]:
    """
    Build a dict of canonical_name → Medium from materials.json.

    Parameters
    ----------
    materials_json
        Path to metadata/materials.json produced by gdml-to-mesh.
    wavelength_min
        Minimum wavelength in nm for the simulation grid.
    wavelength_max
        Maximum wavelength in nm for the simulation grid.
    n_samples
        Number of uniform samples on the wavelength grid.

    Returns
    -------
    Dict mapping canonical material name to a Theia Medium.
    """
    with open(materials_json) as f:
        raw = json.load(f)

    entries = raw.get("materials", raw) if isinstance(raw, dict) else raw

    media: dict[str, Medium] = {}

    for entry in entries:
        if not isinstance(entry, dict):
            continue

        canon = canonical(entry.get("canonical_name") or entry.get("name", ""))
        props = entry.get("properties", {})

        kwargs: dict = {}

        def _prop(arr):
            return TableProperty(Table(arr.astype("float32"), (wavelength_min, wavelength_max)))

        def _prop01(arr):
            return TableProperty(Table(arr.astype("float32"), (0, 1)))

        if "RINDEX" in props:
            p = props["RINDEX"]
            kwargs["refractive_index"] = _prop(_interp(p["wavelength_nm"], p["values"]))

        if "GROUPVEL" in props:
            p = props["GROUPVEL"]
            kwargs["group_velocity"] = _prop(_interp(p["wavelength_nm"], p["values"]) / 1000.0)

        if "ABSLENGTH" in props:
            p = props["ABSLENGTH"]
            kwargs["absorption_coef"] = _prop(_mfp_to_coef(_interp(p["wavelength_nm"], p["values"])))

        if "RAYLEIGH" in props:
            p = props["RAYLEIGH"]
            scatlen = _interp(p["wavelength_nm"], p["values"])
            kwargs["scattering_coef"] = _prop(_mfp_to_coef(scatlen))
            cosines = np.linspace(0, 1, len(scatlen))
            kwargs["log_phase_function"] = _prop01(
                RayleighScatteringPhaseFunction().log_phase_function(cosines)
            )
            kwargs["phase_sampling"] = _prop01(
                RayleighScatteringPhaseFunction().phase_sampling(cosines)
            )

        for n in (1, 2, 3):
            key = f"SCINTILLATIONCOMPONENT{n}"
            if key in props:
                p = props[key]
                kwargs[f"scintillation_spectrum_{n}_sampler"] = _prop01(
                    _compute_ppf(p["wavelength_nm"], p["values"])
                )

        for particle, suffix in (
            ("", ""),
            ("ELECTRON", "_electron"),
            ("PROTON", "_proton"),
            ("ALPHA", "_alpha"),
            ("ION", "_ion"),
            ("DEUTERON", "_deuteron"),
            ("TRITON", "_triton"),
        ):
            for n in (1, 2, 3):
                key = f"{particle}SCINTILLATIONTIMECONSTANT{n}"
                if key in props:
                    kwargs[f"scintillation_time_constant_{n}{suffix}"] = FloatProperty(
                        float(props[key]["value"])
                    )

            for n in (1, 2, 3):
                key = f"{particle}SCINTILLATIONYIELD{n}"
                if key in props:
                    kwargs[f"scintillation_branching_ratio_{n}{suffix}"] = FloatProperty(
                        float(props[key]["value"])
                    )

            key = f"{particle}SCINTILLATIONYIELD"
            if key in props:
                p = props[key]
                scint_yield = p["values"][0] * 806554.394 * p["wavelength_nm"][0] * 1e-9 / u.eV
                kwargs[f"scintillation_yield{suffix}"] = FloatProperty(float(scint_yield))

        if "RESOLUTIONSCALE" in props:
            kwargs["scintillation_resolutionscale"] = FloatProperty(
                float(props["RESOLUTIONSCALE"]["value"])
            )

        for n in (1, 2, 3):
            key = f"SCINTILLATIONRISETIME{n}"
            if key in props:
                kwargs[f"scintillation_rise_time_{n}"] = FloatProperty(float(props[key]["value"]))

        if "WLSABSLENGTH" in props:
            p = props["WLSABSLENGTH"]
            kwargs["fluorescence_coef"] = _prop(_mfp_to_coef(_interp(p["wavelength_nm"], p["values"])))

        if "WLSMEANNUMBERPHOTONS" in props:
            kwargs["fluorescence_efficiency"] = FloatProperty(
                float(props["WLSMEANNUMBERPHOTONS"]["value"])
            )

        if "WLSTIMECONSTANT" in props:
            kwargs["fluorescence_time_shift"] = FloatProperty(
                float(props["WLSTIMECONSTANT"]["value"])
            )

        if "WLSCOMPONENT" in props:
            p = props["WLSCOMPONENT"]
            kwargs["fluorescence_emission_sampling"] = _prop01(
                _compute_ppf(p["wavelength_nm"], p["values"])
            )
            kwargs["fluorescence_emission_quantile"] = _prop01(
                _compute_cdf(p["wavelength_nm"], p["values"])
            )

        if ("WLSABSLENGTH" in props or "WLSMEANNUMBERPHOTONS" in props
                or "WLSTIMECONSTANT" in props or "WLSCOMPONENT" in props):
            volume_model = Fluorescent(timeModel="exponential")
        elif("ABSLENGTH" in props or "RAYLEIGH" in props):
            volume_model = Attenuating()
        else:
            volume_model = Transparent()

        medium = Medium(canon, (wavelength_min, wavelength_max), kwargs, volume_model)
        media[canon] = medium
        print(f"  medium: {canon} ({list(kwargs.keys())})")

    return media
