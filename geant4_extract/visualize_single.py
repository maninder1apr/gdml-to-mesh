#!/usr/bin/env python3
"""
visualize_single.py — inspect one or all interfaces with bounding volumes.

Usage:
    python3 visualize_single.py           # show all volumes + all interfaces
    python3 visualize_single.py --id 0   # highlight interface 0
    python3 visualize_single.py --id 1   # highlight interface 1
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import trimesh

# ============================================================
# styling
# ============================================================

INTERFACE_COLORS = ["#E74C3C", "#3498DB", "#2ECC71", "#F39C12"]

VOLUME_STYLE = {
    "lar_pv":       dict(color="#AED6F1", alpha=0.05, edge="#AED6F133", lw=0.2, label="LAr (liquid argon)"),
    "bege_pv":      dict(color="#A9DFBF", alpha=0.35, edge="#A9DFBF99", lw=0.3, label="BeGe (Germanium)"),
    "pen_bege_pv":  dict(color="#D7BDE2", alpha=0.20, edge="#D7BDE266", lw=0.3, label="PEN encapsulation"),
    "icpc_pv":      dict(color="#A9DFBF", alpha=0.35, edge="#A9DFBF99", lw=0.3, label="ICPC (Germanium)"),
    "pen_icpc_pv":  dict(color="#D7BDE2", alpha=0.20, edge="#D7BDE266", lw=0.3, label="PEN (ICPC)"),
    "sipm_top_0":   dict(color="#F9E79F", alpha=0.50, edge="#F9E79FBB", lw=0.3, label="SiPM top"),
    "sipm_bot_0":   dict(color="#F9E79F", alpha=0.50, edge="#F9E79FBB", lw=0.3, label="SiPM bot"),
}

SURFACE_MARKER = {"blackbody": "●", "specular": "◆", "detector": "★"}

# ============================================================
# helpers
# ============================================================

def load_interfaces(metadata_dir="metadata"):
    path = Path(metadata_dir) / "interfaces.json"
    if not path.exists():
        print(f"ERROR: {path} not found.")
        sys.exit(1)
    with open(path) as f:
        return json.load(f)


def load_mesh(stl_path, max_faces=5000):
    p = Path(stl_path)
    if not p.exists():
        return None, None
    mesh = trimesh.load(str(p))
    if isinstance(mesh, trimesh.Scene):
        mesh = trimesh.util.concatenate(mesh.dump())
    faces = mesh.faces
    if len(faces) > max_faces:
        idx = np.random.choice(len(faces), max_faces, replace=False)
        tris = mesh.vertices[faces[idx]]
    else:
        tris = mesh.vertices[faces]
    return mesh, tris


def add_mesh_to_ax(ax3d, ax_top, ax_sid, tris, color, alpha, edge, lw,
                   step3d=1, step2d=3):
    ax3d.add_collection3d(Poly3DCollection(
        tris[::step3d], alpha=alpha,
        facecolor=color, edgecolor=edge, linewidth=lw, zsort='average'
    ))
    for tri in tris[::step2d]:
        ax_top.fill(tri[:,0], tri[:,1], color=color, alpha=min(alpha*1.5, 1.0), linewidth=0)
        ax_top.plot(np.append(tri[:,0], tri[0,0]),
                   np.append(tri[:,1], tri[0,1]),
                   color=edge, lw=0.15, alpha=0.3)
        ax_sid.fill(tri[:,0], tri[:,2], color=color, alpha=min(alpha*1.5, 1.0), linewidth=0)
        ax_sid.plot(np.append(tri[:,0], tri[0,0]),
                   np.append(tri[:,2], tri[0,2]),
                   color=edge, lw=0.15, alpha=0.3)


def boundary_check(mesh_vol, mesh_iface):
    try:
        dists = trimesh.proximity.closest_point(mesh_vol, mesh_iface.vertices)[1]
        return dists.max(), dists.mean()
    except Exception:
        return None, None

# ============================================================
# render
# ============================================================

def render(interfaces, volumes_dir="cad/volumes", highlight_id=None):

    fig = plt.figure(figsize=(17, 9), facecolor="#0F1923")
    ax3d  = fig.add_subplot(131, projection='3d', facecolor="#0F1923")
    ax_top = fig.add_subplot(132, facecolor="#0F1923")
    ax_sid = fig.add_subplot(133, facecolor="#0F1923")

    legend_patches = []
    all_verts = []
    val_lines = []

    # --------------------------------------------------------
    # 1. volumes — draw in back-to-front order (LAr last = outermost)
    # --------------------------------------------------------
    volume_order = [
        "sipm_bot_0", "sipm_top_0",
        "bege_pv", "icpc_pv",
        "pen_bege_pv", "pen_icpc_pv",
        "lar_pv",
    ]

    for vol_name in volume_order:
        stl_path = Path(volumes_dir) / f"{vol_name}.stl"
        if not stl_path.exists():
            continue
        style = VOLUME_STYLE.get(vol_name, dict(color="#CCCCCC", alpha=0.1, edge="#CCCCCC33", lw=0.2, label=vol_name))
        mesh, tris = load_mesh(str(stl_path), max_faces=2500)
        if tris is None:
            continue
        all_verts.append(mesh.vertices)
        add_mesh_to_ax(ax3d, ax_top, ax_sid, tris,
                       style['color'], style['alpha'], style['edge'], style['lw'],
                       step3d=1, step2d=4)
        legend_patches.append(mpatches.Patch(
            facecolor=style['color'], edgecolor="#ffffff22",
            alpha=0.7, label=f"vol: {style['label']}"
        ))

    # --------------------------------------------------------
    # 2. interfaces — solid, color-coded, on top
    # --------------------------------------------------------
    for iface in interfaces:
        color = INTERFACE_COLORS[iface['id'] % len(INTERFACE_COLORS)]
        marker = SURFACE_MARKER.get(iface['surface'], "?")

        # dim non-highlighted interfaces if a specific one is selected
        alpha = 0.90 if (highlight_id is None or iface['id'] == highlight_id) else 0.25

        mesh_if, tris_if = load_mesh(iface['stl'], max_faces=6000)
        if tris_if is None:
            continue
        all_verts.append(mesh_if.vertices)

        add_mesh_to_ax(ax3d, ax_top, ax_sid, tris_if,
                       color, alpha, "#00000033", 0.15,
                       step3d=1, step2d=2)

        legend_patches.append(mpatches.Patch(
            facecolor=color, edgecolor="#ffffff44",
            label=f"[{iface['id']}] {marker} {iface['pv_inside']} ({iface['surface']})"
        ))

        # boundary check for highlighted or all
        if highlight_id is None or iface['id'] == highlight_id:
            inside_stl = Path(volumes_dir) / f"{iface['pv_inside']}.stl"
            if inside_stl.exists():
                mesh_vol, _ = load_mesh(str(inside_stl), max_faces=99999)
                if mesh_vol:
                    mx, mn = boundary_check(mesh_vol, mesh_if)
                    if mx is not None:
                        status = "PASS" if mx < 0.15 else "FAIL"
                        val_lines.append(
                            (f"[{iface['id']}] {iface['pv_inside']} boundary: "
                             f"max={mx:.4f}mm  mean={mn:.4f}mm  [{status}]",
                             "#2ECC71" if status == "PASS" else "#E74C3C")
                        )

    # --------------------------------------------------------
    # auto-scale
    # --------------------------------------------------------
    if all_verts:
        vv = np.vstack(all_verts)
        mn_v, mx_v = vv.min(axis=0), vv.max(axis=0)
        mid = (mn_v + mx_v) / 2
        rng = (mx_v - mn_v).max() / 2 * 1.25

        ax3d.set_xlim(mid[0]-rng, mid[0]+rng)
        ax3d.set_ylim(mid[1]-rng, mid[1]+rng)
        ax3d.set_zlim(mid[2]-rng, mid[2]+rng)

        for ax2, xi, yi, xl, yl in [
            (ax_top, 0, 1, "X (mm)", "Y (mm)"),
            (ax_sid, 0, 2, "X (mm)", "Z (mm)"),
        ]:
            ax2.set_xlim(mid[xi]-rng, mid[xi]+rng)
            ax2.set_ylim(mid[yi]-rng, mid[yi]+rng)
            ax2.set_xlabel(xl, color="#8899AA", fontsize=9)
            ax2.set_ylabel(yl, color="#8899AA", fontsize=9)
            ax2.tick_params(colors="#8899AA", labelsize=7)
            ax2.set_aspect('equal')
            for spine in ax2.spines.values():
                spine.set_edgecolor("#1E2D3D")

    # 3D axis style
    for axis in [ax3d.xaxis, ax3d.yaxis, ax3d.zaxis]:
        axis.pane.fill = False
        axis.pane.set_edgecolor("#1E2D3D")
    ax3d.tick_params(colors="#8899AA", labelsize=7)
    ax3d.set_xlabel("X (mm)", color="#8899AA", labelpad=6)
    ax3d.set_ylabel("Y (mm)", color="#8899AA", labelpad=6)
    ax3d.set_zlabel("Z (mm)", color="#8899AA", labelpad=6)
    ax3d.set_title("3D View", color="#E8EDF2", pad=10, fontsize=10)
    ax_top.set_title("Top view (XY)", color="#E8EDF2", pad=10, fontsize=10)
    ax_sid.set_title("Side view (XZ)", color="#E8EDF2", pad=10, fontsize=10)

    ax3d.legend(
        handles=legend_patches,
        loc="upper left",
        facecolor="#131F2B",
        edgecolor="#2E5FA3",
        labelcolor="#E8EDF2",
        fontsize=7.5,
        framealpha=0.92,
    )

    # boundary validation footer
    if val_lines:
        text = "\n".join(t for t, _ in val_lines)
        color = val_lines[0][1]
        fig.text(0.01, 0.01, text, color=color,
                fontsize=8, fontfamily="monospace",
                verticalalignment="bottom")

    title = "LEGEND-Theia v0  ·  scarf_pen.gdml  ·  LAr / PEN / BeGe + Interfaces"
    if highlight_id is not None:
        iface = next(i for i in interfaces if i['id'] == highlight_id)
        marker = SURFACE_MARKER.get(iface['surface'], "?")
        title = (f"Interface {highlight_id} {marker}  {iface['pv_inside']} | {iface['pv_outside']}"
                 f"  —  {iface['surface'].upper()}  ·  {iface['n_triangles']:,} tris")

    plt.suptitle(title, color="#E8EDF2", fontsize=10, y=1.01)
    plt.tight_layout()
    plt.show()


# ============================================================
# main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Geometry + interface inspector")
    parser.add_argument("--id",      type=int, default=None, help="Highlight one interface (default: show all)")
    parser.add_argument("--metadata",default="metadata",    help="Metadata directory")
    parser.add_argument("--volumes", default="cad/volumes", help="Volume STL directory")
    args = parser.parse_args()

    interfaces = load_interfaces(args.metadata)
    render(interfaces, args.volumes, highlight_id=args.id)


if __name__ == "__main__":
    main()
