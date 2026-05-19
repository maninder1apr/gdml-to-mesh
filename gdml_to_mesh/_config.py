"""
_config.py — auto-detect paths for Geant4, OCC, HDF5, Qt5.

Priority order for each dependency:
  1. Environment variable (e.g. GEANT4_INSTALL)
  2. Known default locations (from the build cache)
  3. Homebrew / system paths
"""

import os
import shutil
from pathlib import Path

# ============================================================
# known default install paths (from CMakeCache.txt)
# ============================================================

_DEFAULTS = {
    "geant4": [
        os.environ.get("GEANT4_INSTALL", ""),
        "/Users/maninder/Desktop/Programs/geant4-install",
        "/usr/local/geant4",
        "/opt/geant4",
    ],
    "occ": [
        os.environ.get("OCC_INSTALL", ""),
        "/opt/homebrew/Cellar/opencascade/7.9.3",
        "/opt/homebrew/opt/opencascade",
        "/usr/local/opencascade",
    ],
    "hdf5": [
        os.environ.get("HDF5_INSTALL", ""),
        "/Users/maninder/Desktop/Programs/hdf5-install",
        "/opt/homebrew/opt/hdf5",
        "/usr/local/hdf5",
    ],
    "qt5": [
        os.environ.get("Qt5_DIR", ""),
        "/opt/homebrew/opt/qt@5",
        "/opt/homebrew/opt/qt5",
    ],
}


def find_path(name: str, marker: str) -> Path | None:
    """Find a dependency by checking default paths for a marker file."""
    for candidate in _DEFAULTS[name]:
        if not candidate:
            continue
        p = Path(candidate)
        if (p / marker).exists():
            return p
    return None


def geant4_dir() -> Path:
    p = find_path("geant4", "bin/geant4-config")
    if p is None:
        raise RuntimeError(
            "Geant4 not found. Set GEANT4_INSTALL=/path/to/geant4-install "
            "or install Geant4 to a standard location."
        )
    return p


def occ_dir() -> Path:
    p = find_path("occ", "include/opencascade/TopoDS_Shape.hxx")
    if p is None:
        raise RuntimeError(
            "OpenCASCADE not found. Set OCC_INSTALL=/path/to/occ "
            "or install via: brew install opencascade"
        )
    return p


def hdf5_dir() -> Path:
    p = find_path("hdf5", "cmake/hdf5-config.cmake")
    if p is None:
        # try cmake subdir
        p = find_path("hdf5", "lib/cmake/hdf5/hdf5-config.cmake")
    if p is None:
        raise RuntimeError(
            "HDF5 not found. Set HDF5_INSTALL=/path/to/hdf5-install "
            "or install via: brew install hdf5"
        )
    return p


def qt5_dir() -> Path:
    p = find_path("qt5", "lib/cmake/Qt5Core")
    if p is None:
        raise RuntimeError(
            "Qt5 not found. Set Qt5_DIR=/path/to/qt5 "
            "or install via: brew install qt@5"
        )
    return p


def cmake_bin() -> str:
    c = shutil.which("cmake")
    if c is None:
        raise RuntimeError("cmake not found. Install via: brew install cmake")
    return c


def print_config():
    """Print detected dependency paths."""
    print("gdml-to-mesh dependency detection:")
    for name, finder in [
        ("Geant4", geant4_dir),
        ("OCC",    occ_dir),
        ("HDF5",   hdf5_dir),
        ("Qt5",    qt5_dir),
    ]:
        try:
            p = finder()
            print(f"  {name:8s}: {p}")
        except RuntimeError as e:
            print(f"  {name:8s}: NOT FOUND — {e}")
