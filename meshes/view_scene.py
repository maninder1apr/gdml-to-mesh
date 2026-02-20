import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ------------------------------------------------------------
# Resolve paths robustly
# ------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
MESH_DIR = PROJECT_ROOT / "meshes"

# ------------------------------------------------------------
# Mesh table (must match extract_placements.cc)
# ------------------------------------------------------------
MESH_TABLE = {
    0: MESH_DIR / "lar/lar_s.npz",
    1: MESH_DIR / "bege/bege_lv.npz",
    2: MESH_DIR / "icpc/icpc_lv.npz",
}

COLORS = {
    0: "lightblue",   # LAr
    1: "red",         # BEGe
    2: "green",       # ICPC
}

# ------------------------------------------------------------
# Load placements
# ------------------------------------------------------------
p = np.load(MESH_DIR / "placements.npz")

translations = p["translations"]      # (N,3)
rotations = p["rotations"]            # (N,9)
mesh_ids = p["mesh_ids"]

print(f"Loaded {len(mesh_ids)} instances")

# ------------------------------------------------------------
# Setup matplotlib scene
# ------------------------------------------------------------
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection="3d")

all_points = []

# ------------------------------------------------------------
# Draw each instance
# ------------------------------------------------------------
for i, mid in enumerate(mesh_ids):
    mesh_path = MESH_TABLE[mid]
    mesh = np.load(mesh_path)

    V = mesh["vertices"]     # (Nv,3)
    F = mesh["indices"]      # (Nt,3)

    R = rotations[i].reshape(3, 3)
    T = translations[i]

    # Apply transform: v' = R v + T
    Vt = (V @ R.T) + T

    all_points.append(Vt)

    for tri in F:
        pts = Vt[tri]
        ax.plot(
            pts[[0, 1, 2, 0], 0],
            pts[[0, 1, 2, 0], 1],
            pts[[0, 1, 2, 0], 2],
            color=COLORS.get(mid, "black"),
            linewidth=0.6,
            alpha=0.85
        )

# ------------------------------------------------------------
# Equal aspect ratio
# ------------------------------------------------------------
all_points = np.vstack(all_points)
mins = all_points.min(axis=0)
maxs = all_points.max(axis=0)

center = 0.5 * (mins + maxs)
extent = 0.5 * np.max(maxs - mins)

ax.set_xlim(center[0] - extent, center[0] + extent)
ax.set_ylim(center[1] - extent, center[1] + extent)
ax.set_zlim(center[2] - extent, center[2] + extent)

ax.set_box_aspect([1, 1, 1])

# ------------------------------------------------------------
# Labels & title
# ------------------------------------------------------------
ax.set_xlabel("X [mm]")
ax.set_ylabel("Y [mm]")
ax.set_zlabel("Z [mm]")
ax.set_title("SCARF geometry (meshes + placements)")

plt.tight_layout()
plt.show()