# gdml-to-mesh

Converts GDML detector geometry files into triangulated mesh outputs for optical simulation and visualization. Built for the [LEGEND-Theia](https://legend-exp.org) optical photon transport pipeline.

---

## Requirements

### System dependencies (macOS)

Install via Homebrew:

```bash
brew install cmake opencascade qt@5
```

Install Geant4 and HDF5 manually (or from your experiment software stack):

| Dependency | Version | Default path |
|-----------|---------|-------------|
| Geant4 | ≥ 11 (with GDML) | `~/Desktop/Programs/geant4-install` |
| OpenCASCADE | ≥ 7.6 | `/opt/homebrew/Cellar/opencascade/7.9.3` |
| HDF5 | any | `~/Desktop/Programs/hdf5-install` |
| Qt5 | 5.x | `/opt/homebrew/opt/qt@5` |
| CMake | ≥ 3.16 | via Homebrew |

If your paths differ from the defaults, set environment variables:

```bash
export GEANT4_INSTALL=/path/to/geant4-install
export OCC_INSTALL=/path/to/opencascade
export HDF5_INSTALL=/path/to/hdf5-install
export Qt5_DIR=/path/to/qt5
```

### Python dependencies

- Python 3.9+
- `trimesh`, `matplotlib`, `numpy`, `scipy` (installed automatically by pip)

---

## Installation

```bash
git clone https://github.com/maninder1apr/gdml-to-mesh.git
cd gdml-to-mesh

# create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# install the Python package (builds the C++ binary automatically)
pip install -e .
```

The first `pip install` will:
1. Detect your Geant4, OCC, HDF5, Qt5 paths
2. Run `cmake` and `make` to build the `occ_mesher` binary
3. Install the `gdml-to-mesh` CLI and Python package

---

## Verify installation

```bash
gdml-to-mesh check
```

Expected output:
```
gdml-to-mesh dependency detection:
  Geant4  : /path/to/geant4-install
  OCC     : /path/to/opencascade
  HDF5    : /path/to/hdf5-install
  Qt5     : /path/to/qt5

  binary built : True
  binary path  : .../geant4_extract/build/occ_mesher
```

---

## Usage

### Command line

```bash
# run the geometry engine on a GDML file
gdml-to-mesh run gdml/scarf_pen.gdml

# with verbose output
gdml-to-mesh run gdml/scarf_pen.gdml --verbose

# write outputs to a specific directory
gdml-to-mesh run gdml/scarf_pen.gdml --output-dir ./output

# rebuild the C++ binary
gdml-to-mesh build

# launch the visualizer
gdml-to-mesh visualize
gdml-to-mesh visualize --surface detector
gdml-to-mesh visualize --surface blackbody
gdml-to-mesh visualize --surface specular
```

### Python API

```python
from gdml_to_mesh import run, load_interfaces, load_materials

# run the geometry engine
result = run("gdml/scarf_pen.gdml")

print(result.summary())
# GeometryResult
#   interfaces : 298
#     specular    : 272
#     blackbody   : 14
#     detector    : 12
#   materials  : 9
#   surfaces   : 554

# access the data
for iface in result.interfaces:
    print(iface["lv_inside"], iface["surface"], iface["area_mm2"])

# load previously generated outputs
result = load_interfaces("output/metadata/interfaces.json")
```

---

## Output structure

```
cad/
  interfaces/
    interface_0.stl       # LAr ↔ BeGe boundary (blackbody)
    interface_1.stl       # LAr ↔ PEN boundary (specular)
    ...                   # 298 total for scarf_pen.gdml
  volumes/
    lar_pv.stl            # key detector volumes for visualization
    bege_pv.stl
    pen_bege_pv.stl
    icpc_pv.stl
    sipm_top_0.stl
    ...
  detector.brep

metadata/
  interfaces.json         # one entry per STL
  materials.json          # optical MPT tables with canonical LEGEND names
  surfaces.json           # 554 G4LogicalBorderSurface entries with MPT
```

### Interface surface types

| Type | Volumes | Count |
|------|---------|-------|
| `blackbody` | HPGe (BeGe, ICPC) + SiPM Cu wraps | 14 |
| `specular` | PEN encapsulations + TPB fiber coatings | 272 |
| `detector` | SiPM silicon photocathodes | 12 |

---

## Visualizer

```bash
# all interfaces + detector volumes (color-coded by surface type)
python3 geant4_extract/visualize.py

# filter by surface type
python3 geant4_extract/visualize.py --surface detector
python3 geant4_extract/visualize.py --surface blackbody
python3 geant4_extract/visualize.py --surface specular

# inspect a single interface with boundary validation
python3 geant4_extract/visualize_single.py --id 0
python3 geant4_extract/visualize_single.py --id 1
```

Colors: 🔴 blackbody · 🔵 specular · 🟢 detector

---

## Metadata format

### `interfaces.json`

```json
[
  {
    "id": 0,
    "stl": "cad/interfaces/interface_0.stl",
    "lv_inside": "bege_pv",
    "lv_outside": "lar_pv",
    "material_inside": "EnrichedGermanium0.076",
    "material_outside": "liquid_argon",
    "surface": "blackbody",
    "detector_id": null,
    "n_triangles": 2688,
    "area_mm2": 15216.9
  }
]
```

### `materials.json`

Exports all Geant4 optical material property tables (RINDEX, RAYLEIGH, ABSLENGTH, GROUPVEL, scintillation components) with a `canonical_name` field mapping G4 material names to LEGEND registry keys.

### `surfaces.json`

Exports all `G4LogicalBorderSurface` and `G4LogicalSkinSurface` entries with model, finish, type, and full wavelength-indexed MPT tables. 554 surfaces for `scarf_pen.gdml`.

---

## Surface classification

| Material | Surface type | Notes |
|----------|-------------|-------|
| `EnrichedGermanium0.076`, `NaturalGermanium` | `blackbody` | HPGe approximation |
| `PEN` | `specular` | WLS not yet implemented |
| `tpb_on_fibers` | `specular` | Diffuse/WLS not yet implemented |
| `metal_silicon` | `detector` | SiPM photocathode |
| `metal_copper` | `blackbody` | SiPM copper wrap |

---

## License

See [LICENSE](LICENSE).
