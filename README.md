# gdml-to-mesh

A reproducible pipeline to convert **Geant4 GDML detector geometry** into
**GPU-ready triangle meshes and placement data**, suitable for
**Theia / Hephaistos** ray tracing and optical simulations.

The pipeline preserves:
- exact geometry
- instancing
- placements (translations + rotations)
- material classification

and avoids baking transforms into meshes.

---

## Overview

This project implements the following workflow:

1. Load detector geometry from **Geant4 GDML**
2. Extract **triangle meshes** for each unique solid
3. Extract **placements** (translation + rotation) for each physical volume
4. Export data to **NumPy `.npz` files**
5. Visually validate geometry + placements in Python
6. Prepare geometry for GPU loading with **Theia**

---

## Repository Structure
