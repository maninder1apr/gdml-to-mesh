import numpy as np
import vtk
from pathlib import Path

from pyg4ometry.gdml import Reader
from pyg4ometry.visualisation import VtkViewer


# ============================================================
# CONFIG
# ============================================================
GDML_FILE = "/Users/maninder/Desktop/Simulations/pygeom-scarf/tests/scarf_pen.gdml"
OUTDIR = Path("meshes_npz")

vtk.vtkPolyDataMapper.SetResolveCoincidentTopologyToPolygonOffset()


# ============================================================
# 1. Load geometry from GDML
# ============================================================
reader = Reader(GDML_FILE)
reg = reader.getRegistry()

print("\n========== LOGICAL VOLUMES IN REGISTRY ==========")
for i, lv in enumerate(reg.logicalVolumeDict.values()):
    solid = lv.solid
    print(
        f"[{i:02d}] LV='{lv.name}' | "
        f"Solid='{solid.name}' | "
        f"Type={type(solid).__name__}"
    )
print("================================================\n")


# ============================================================
# 2. Create VTK viewer (for meshes)
# ============================================================
viewer = VtkViewer()
viewer.addLogicalVolume(reg.getWorldVolume())


# ============================================================
# 3. Prepare output directory
# ============================================================
OUTDIR.mkdir(exist_ok=True)


# ============================================================
# 4. Export meshes (logical volumes)
# ============================================================
mesh_names = set()
written = 0

print("\nActors found:")
for i, actor in enumerate(viewer.actors):
    print(f"  {i}: {actor.GetObjectName()}")

for i, actor in enumerate(viewer.actors):

    mapper = actor.GetMapper()
    if mapper is None:
        continue

    polydata = mapper.GetInput()
    if polydata is None:
        continue

    points = polydata.GetPoints()
    if points is None or points.GetNumberOfPoints() == 0:
        continue

    vertices = np.array(
        [points.GetPoint(j) for j in range(points.GetNumberOfPoints())],
        dtype=np.float32,
    )

    polys = polydata.GetPolys()
    polys.InitTraversal()

    faces = []
    ids = vtk.vtkIdList()

    while polys.GetNextCell(ids):
        if ids.GetNumberOfIds() == 3:
            faces.append([ids.GetId(0), ids.GetId(1), ids.GetId(2)])

    if not faces:
        continue

    indices = np.array(faces, dtype=np.int32)

    name = actor.GetObjectName()
    if not name:
        name = f"mesh_{i}"

    np.savez(
        OUTDIR / f"{name}.npz",
        vertices=vertices,
        indices=indices,
    )

    mesh_names.add(name)

    print(f"✅ wrote {name}.npz "
          f"({vertices.shape[0]} verts, {indices.shape[0]} tris)")
    written += 1


# ============================================================
# 5. Extract placements (physical volumes)
# ============================================================
translations = []
rotations = []
lv_names = []

for pv in reg.physicalVolumeDict.values():
    lv = pv.logicalVolume
    lv_name = lv.name

    # only keep placements for which we have a mesh
    if lv_name not in mesh_names:
        continue

    t = pv.translation
    r = pv.rotation

    translations.append([t[0], t[1], t[2]])

    if r is None:
        rotations.append([1, 0, 0,
                          0, 1, 0,
                          0, 0, 1])
    else:
        rotations.append([
            r[0][0], r[0][1], r[0][2],
            r[1][0], r[1][1], r[1][2],
            r[2][0], r[2][1], r[2][2],
        ])

    lv_names.append(lv_name)


np.savez(
    OUTDIR / "placements.npz",
    translations=np.array(translations, dtype=np.float32),
    rotations=np.array(rotations, dtype=np.float32),
    lv_names=np.array(lv_names),
)


# ============================================================
# SUMMARY
# ============================================================
print("\n🎉 DONE")
print(f"  Meshes written       : {written}")
print(f"  Placement instances  : {len(lv_names)}")
print(f"  Output directory     : {OUTDIR.resolve()}")
