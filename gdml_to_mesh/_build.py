"""
_build.py — builds the occ_mesher C++ binary via cmake.
Called automatically by run() if the binary is not found.
"""

import os
import subprocess
import sys
from pathlib import Path

from . import _config

# repo root is two levels up from this file
_REPO_ROOT    = Path(__file__).parent.parent
_SOURCE_DIR   = _REPO_ROOT / "geant4_extract"
_BUILD_DIR    = _SOURCE_DIR / "build"
_BINARY       = _BUILD_DIR / "occ_mesher"


def binary_path() -> Path:
    return _BINARY


def is_built() -> bool:
    return _BINARY.exists()


def build(verbose: bool = False, jobs: int = 8) -> None:
    """Configure and build the occ_mesher binary."""

    g4   = _config.geant4_dir()
    occ  = _config.occ_dir()
    hdf5 = _config.hdf5_dir()
    qt5  = _config.qt5_dir()
    cmake = _config.cmake_bin()

    _BUILD_DIR.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # cmake configure
    # --------------------------------------------------------
    configure_cmd = [
        cmake,
        str(_SOURCE_DIR),
        f"-DGeant4_DIR={g4}/lib/cmake/Geant4",
        f"-DOpenCASCADE_DIR={occ}/lib/cmake/opencascade",
        f"-DHDF5_DIR={hdf5}/cmake",
        f"-DQt5_DIR={qt5}/lib/cmake/Qt5",
        "-DCMAKE_BUILD_TYPE=Release",
    ]

    print("gdml-to-mesh: configuring C++ build...")
    if verbose:
        print("  " + " ".join(configure_cmd))

    result = subprocess.run(
        configure_cmd,
        cwd=str(_BUILD_DIR),
        capture_output=not verbose,
        text=True,
    )

    if result.returncode != 0:
        print("CMake configure failed:")
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("CMake configure failed. See output above.")

    # --------------------------------------------------------
    # cmake build
    # --------------------------------------------------------
    build_cmd = [cmake, "--build", ".", f"-j{jobs}", "--target", "occ_mesher"]

    print(f"gdml-to-mesh: building occ_mesher (jobs={jobs})...")

    result = subprocess.run(
        build_cmd,
        cwd=str(_BUILD_DIR),
        capture_output=not verbose,
        text=True,
    )

    if result.returncode != 0:
        print("Build failed:")
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("Build failed. Run with verbose=True for details.")

    print(f"gdml-to-mesh: built successfully → {_BINARY}")


def ensure_built(verbose: bool = False) -> Path:
    """Build if not already built. Returns path to binary."""
    if not is_built():
        print("gdml-to-mesh: occ_mesher binary not found, building...")
        build(verbose=verbose)
    return _BINARY
