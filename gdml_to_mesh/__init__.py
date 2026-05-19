"""
gdml_to_mesh — Python wrapper for the LEGEND-Theia GDML geometry engine.

Usage
-----
    from gdml_to_mesh import run, load_interfaces, load_materials

    # run the geometry engine on a GDML file
    result = run("path/to/detector.gdml")

    # result.interfaces  — list of interface dicts
    # result.materials   — materials dict
    # result.surfaces    — list of surface dicts
    # result.output_dir  — Path to cad/ and metadata/ outputs

    # or load previously generated outputs directly
    interfaces = load_interfaces("path/to/metadata/interfaces.json")
"""

import json
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from . import _build


# ============================================================
# result container
# ============================================================

@dataclass
class GeometryResult:
    output_dir:  Path
    interfaces:  list  = field(default_factory=list)
    materials:   dict  = field(default_factory=dict)
    surfaces:    list  = field(default_factory=list)

    @property
    def n_interfaces(self) -> int:
        return len(self.interfaces)

    @property
    def interface_stls(self) -> list[Path]:
        return [self.output_dir / i["stl"] for i in self.interfaces]

    def summary(self) -> str:
        from collections import Counter
        by_surface = Counter(i["surface"] for i in self.interfaces)
        lines = [
            f"GeometryResult — {self.output_dir}",
            f"  interfaces : {self.n_interfaces}",
        ]
        for s, n in by_surface.most_common():
            lines.append(f"    {s:12s}: {n}")
        lines.append(f"  materials  : {len(self.materials.get('materials', []))}")
        lines.append(f"  surfaces   : {len(self.surfaces)}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()


# ============================================================
# public API
# ============================================================

def run(
    gdml_file: str | Path,
    output_dir: Optional[str | Path] = None,
    rebuild: bool = False,
    verbose: bool = False,
    jobs: int = 8,
) -> GeometryResult:
    """
    Run the occ_mesher on a GDML file and return a GeometryResult.

    Parameters
    ----------
    gdml_file   : path to the .gdml file
    output_dir  : where to write cad/ and metadata/ (default: cwd)
    rebuild     : force rebuild of the C++ binary
    verbose     : show cmake/make output
    jobs        : parallel build jobs

    Returns
    -------
    GeometryResult with loaded interfaces, materials, surfaces
    """

    gdml_path = Path(gdml_file).resolve()
    if not gdml_path.exists():
        raise FileNotFoundError(f"GDML file not found: {gdml_path}")

    # build binary if needed
    if rebuild:
        _build.build(verbose=verbose, jobs=jobs)
    binary = _build.ensure_built(verbose=verbose)

    # set working directory for outputs
    work_dir = Path(output_dir).resolve() if output_dir else Path.cwd()
    work_dir.mkdir(parents=True, exist_ok=True)

    # run the mesher
    print(f"gdml-to-mesh: running on {gdml_path.name}...")
    result = subprocess.run(
        [str(binary), str(gdml_path)],
        cwd=str(work_dir),
        capture_output=not verbose,
        text=True,
    )

    if result.returncode != 0:
        print(result.stdout[-3000:] if result.stdout else "")
        print(result.stderr[-1000:] if result.stderr else "")
        raise RuntimeError(f"occ_mesher failed with exit code {result.returncode}")

    # load outputs
    return _load_result(work_dir)


def load_result(output_dir: str | Path) -> GeometryResult:
    """Load a previously generated geometry result from output_dir."""
    return _load_result(Path(output_dir).resolve())


def load_interfaces(interfaces_json: str | Path) -> list:
    """Load interfaces.json directly."""
    with open(interfaces_json) as f:
        return json.load(f)


def load_materials(materials_json: str | Path) -> dict:
    """Load materials.json directly."""
    with open(materials_json) as f:
        return json.load(f)


def load_surfaces(surfaces_json: str | Path) -> list:
    """Load surfaces.json directly."""
    with open(surfaces_json) as f:
        return json.load(f)


def check_dependencies() -> None:
    """Print detected dependency paths. Useful for debugging installs."""
    from . import _config
    _config.print_config()


# ============================================================
# internal helpers
# ============================================================

def _load_result(work_dir: Path) -> GeometryResult:
    meta = work_dir / "metadata"

    interfaces = []
    materials  = {}
    surfaces   = []

    if (meta / "interfaces.json").exists():
        with open(meta / "interfaces.json") as f:
            interfaces = json.load(f)

    if (meta / "materials.json").exists():
        with open(meta / "materials.json") as f:
            materials = json.load(f)

    if (meta / "surfaces.json").exists():
        with open(meta / "surfaces.json") as f:
            surfaces = json.load(f)

    return GeometryResult(
        output_dir=work_dir,
        interfaces=interfaces,
        materials=materials,
        surfaces=surfaces,
    )
