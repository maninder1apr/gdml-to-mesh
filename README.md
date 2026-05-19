# gdml-to-mesh

Converts GDML detector geometry files into triangulated mesh outputs for optical simulation and visualization. Built for the [LEGEND-Theia](https://legend-exp.org) optical photon transport pipeline.

---

## What it does

Given a GDML file (e.g. `scarf_pen.gdml`), the tool:

1. Parses the full detector geometry using Geant4's GDML parser
2. Converts every logical volume solid into an OCC CAD shape with global transforms baked in
3. Identifies all touching volume pairs (LAr ↔ germanium, LAr ↔ PEN, LAr ↔ fibers, LAr ↔ SiPMs)
4. Computes the shared boundary surface for each pair via OCC Boolean intersection
5. Exports each boundary as a triangulated STL mesh
6. Emits three metadata JSONs that fully describe the geometry handoff to the scene builder

The output is a self-contained directory that the downstream Python scene builder (`legend_theia.geometry`) reads without needing Geant4 or OCC.

---

## Output structure
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
detector.brep           # full OCC assembly
metadata/
interfaces.json         # one entry per STL: lv_inside, lv_outside,
# material, surface type, n_triangles, area_mm2
materials.json          # optical MPT tables with canonical LEGEND names
surfaces.json           # 554 G4LogicalBorderSurface entries with full MPT

### Interface surface types (v0)

| Type | Volumes | Count |
|------|---------|-------|
| `blackbody` | HPGe (BeGe, ICPC) + SiPM Cu wraps | 14 |
| `specular` | PEN encapsulations + TPB fiber coatings | 272 |
| `detector` | SiPM silicon photocathodes | 12 |

WLS (PEN/TPB) and diffuse (PTFE/Tetratex) surfaces are deferred to v1 — treated as specular in v0. The `surfaces.json` already exports all the MPT data needed for v1.

---

## Prerequisites

- **Geant4** ≥ 11 with GDML support (`-DGEANT4_USE_GDML=ON`)
- **OpenCASCADE** ≥ 7.6
- **CMake** ≥ 3.16
- **C++17** compiler
- **nlohmann/json** (fetched automatically by CMake)

For the Python visualizer:
- Python 3.9+
- `trimesh`, `matplotlib`, `numpy`, `scipy`

---

## Build

```bash
git clone https://github.com/maninder1apr/gdml-to-mesh.git
cd gdml-to-mesh/geant4_extract

mkdir build && cd build
cmake ..
make -j8
```

---

## Run

```bash
cd geant4_extract
./build/occ_mesher ../gdml/scarf_pen.gdml
```

Output is written to `cad/` and `metadata/` in the current directory.

---

## Visualizer

### Setup

```bash
cd geant4_extract
python3 -m venv .venv
source .venv/bin/activate
pip install trimesh matplotlib numpy scipy
```

### Usage

```bash
# all 298 interfaces + detector volumes
python3 visualize.py

# stats only
python3 visualize.py --list

# filter by surface type
python3 visualize.py --surface detector
python3 visualize.py --surface blackbody
python3 visualize.py --surface specular

# filter by material
python3 visualize.py --material tpb_on_fibers --no-volumes

# inspect a single interface with boundary validation
python3 visualize_single.py --id 0
python3 visualize_single.py --id 1
```

Interfaces are color-coded by surface type: 🔴 blackbody · 🔵 specular · 🟢 detector

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

Exports all Geant4 optical material property tables (RINDEX, RAYLEIGH, ABSLENGTH, GROUPVEL, scintillation components) with a `canonical_name` field mapping G4 material names to the LEGEND registry keys used by `legend_theia`.

### `surfaces.json`

Exports all `G4LogicalBorderSurface` and `G4LogicalSkinSurface` entries with model, finish, type, and full wavelength-indexed MPT tables. 554 surfaces for `scarf_pen.gdml`.

---

## Surface classification

| Material | Surface type | Notes |
|----------|-------------|-------|
| `EnrichedGermanium0.076`, `NaturalGermanium` | `blackbody` | v0 approximation |
| `PEN` | `specular` | WLS deferred to v1 |
| `tpb_on_fibers` | `specular` | Diffuse/WLS deferred to v1 |
| `metal_silicon` | `detector` | SiPM photocathode |
| `metal_copper` | `blackbody` | SiPM copper wrap |

---

## Next steps

- `legend_theia.geometry` — Python scene builder that consumes this output
- `legend_theia.scint` — LAr scintillation generator from LH5 energy deposits
- v1: Add WLS and diffuse surface kinds to Theia using `surfaces.json` MPT data

---

## License

See [LICENSE](LICENSE).You said: we can remove this v1 and v0 classification in readmewe can remove this v1 and v0 classification in readme3:11 PMClaude responded: Here it is without the v0/v1 references:Here it is without the v0/v1 references:
markdown# gdml-to-mesh

Converts GDML detector geometry files into triangulated mesh outputs for optical simulation and visualization. Built for the [LEGEND-Theia](https://legend-exp.org) optical photon transport pipeline.

---

## What it does

Given a GDML file (e.g. `scarf_pen.gdml`), the tool:

1. Parses the full detector geometry using Geant4's GDML parser
2. Converts every logical volume solid into an OCC CAD shape with global transforms baked in
3. Identifies all touching volume pairs (LAr ↔ germanium, LAr ↔ PEN, LAr ↔ fibers, LAr ↔ SiPMs)
4. Computes the shared boundary surface for each pair via OCC Boolean intersection
5. Exports each boundary as a triangulated STL mesh
6. Emits three metadata JSONs that fully describe the geometry handoff to the scene builder

The output is a self-contained directory that the downstream Python scene builder (`legend_theia.geometry`) reads without needing Geant4 or OCC.

---

## Output structure
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
detector.brep           # full OCC assembly
metadata/
interfaces.json         # one entry per STL: lv_inside, lv_outside,
# material, surface type, n_triangles, area_mm2
materials.json          # optical MPT tables with canonical LEGEND names
surfaces.json           # 554 G4LogicalBorderSurface entries with full MPT

### Interface surface types

| Type | Volumes | Count |
|------|---------|-------|
| `blackbody` | HPGe (BeGe, ICPC) + SiPM Cu wraps | 14 |
| `specular` | PEN encapsulations + TPB fiber coatings | 272 |
| `detector` | SiPM silicon photocathodes | 12 |

---

## Prerequisites

- **Geant4** ≥ 11 with GDML support (`-DGEANT4_USE_GDML=ON`)
- **OpenCASCADE** ≥ 7.6
- **CMake** ≥ 3.16
- **C++17** compiler
- **nlohmann/json** (fetched automatically by CMake)

For the Python visualizer:
- Python 3.9+
- `trimesh`, `matplotlib`, `numpy`, `scipy`

---

## Build

```bash
git clone https://github.com/maninder1apr/gdml-to-mesh.git
cd gdml-to-mesh/geant4_extract

mkdir build && cd build
cmake ..
make -j8
```

---

## Run

```bash
cd geant4_extract
./build/occ_mesher ../gdml/scarf_pen.gdml
```

Output is written to `cad/` and `metadata/` in the current directory.

---

## Visualizer

### Setup

```bash
cd geant4_extract
python3 -m venv .venv
source .venv/bin/activate
pip install trimesh matplotlib numpy scipy
```

### Usage

```bash
# all interfaces + detector volumes
python3 visualize.py

# stats only
python3 visualize.py --list

# filter by surface type
python3 visualize.py --surface detector
python3 visualize.py --surface blackbody
python3 visualize.py --surface specular

# filter by material
python3 visualize.py --material tpb_on_fibers --no-volumes

# inspect a single interface with boundary validation
python3 visualize_single.py --id 0
python3 visualize_single.py --id 1
```

Interfaces are color-coded by surface type: 🔴 blackbody · 🔵 specular · 🟢 detector

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

Exports all Geant4 optical material property tables (RINDEX, RAYLEIGH, ABSLENGTH, GROUPVEL, scintillation components) with a `canonical_name` field mapping G4 material names to the LEGEND registry keys used by `legend_theia`.

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
