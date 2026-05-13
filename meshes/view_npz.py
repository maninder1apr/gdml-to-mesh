import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def visualize_npz(path):
    data = np.load(path)
    vertices = data["vertices"]
    indices = data["indices"]

    print(f"Loaded {path}")
    print(f"  vertices: {vertices.shape}")
    print(f"  triangles: {indices.shape}")

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")

    # draw triangles
    for tri in indices:
        pts = vertices[tri]
        ax.plot_trisurf(
            pts[:, 0],
            pts[:, 1],
            pts[:, 2],
            triangles=[[0, 1, 2]],
            color="lightblue",
            edgecolor="k",
            linewidth=0.3,
            alpha=0.8,
        )

    # axis labels
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_zlabel("Z [mm]")

    # equal aspect ratio
    mins = vertices.min(axis=0)
    maxs = vertices.max(axis=0)
    center = 0.5 * (mins + maxs)
    extent = (maxs - mins).max() * 0.5

    ax.set_xlim(center[0] - extent, center[0] + extent)
    ax.set_ylim(center[1] - extent, center[1] + extent)
    ax.set_zlim(center[2] - extent, center[2] + extent)

    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python view_npz.py <mesh.npz>")
        sys.exit(1)

    visualize_npz(sys.argv[1])
