#!/usr/bin/env python3
"""
visualize.py — LEGEND-Theia geometry visualizer
Shows all 298 interfaces color-coded by surface type.

Usage:
    python3 visualize.py                        # all interfaces + volumes
    python3 visualize.py --list                 # stats only
    python3 visualize.py --no-volumes           # interfaces only
    python3 visualize.py --surface blackbody    # filter by surface type
    python3 visualize.py --surface detector
    python3 visualize.py --surface specular
    python3 visualize.py --material tpb_on_fibers
"""

import argparse
import json
import sys
from pathlib import Path
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import trimesh

# ============================================================
# styling — color by surface type
# ============================================================

SURFACE_COLORS = {
    "blackbody": ("#E74C3C", 0.85),   # red   — germanium + Cu wraps
    "specular":  ("#3498DB", 0.55),   # blue  — PEN + TPB fibers
    "detector":  ("#2ECC71", 0.90),   # green — SiPM photocathodes
}

VOLUME_STYLE = {
    "lar_pv":      dict(color="#AED6F1", alpha=0.04, edge="#AED6F122", lw=0.2, label="LAr"),
    "bege_pv":     dict(color="#A9DFBF", alpha=0.30, edge="#A9DFBF88", lw=0.3, label="BeGe"),
    "pen_bege_pv": dict(color="#D7BDE2", alpha=0.18, edge="#D7BDE255", lw=0.3, label="PEN (BeGe)"),
    "icpc_pv":     dict(color="#A9DFBF", alpha=0.30, edge="#A9DFBF88", lw=0.3, label="ICPC"),
    "pen_icpc_pv": dict(color="#D7BDE2", alpha=0.18, edge="#D7BDE255", lw=0.3, label="PEN (ICPC)"),
    "sipm_top_0":  dict(color="#F9E79F", alpha=0.45, edge="#F9E79FBB", lw=0.3, label="SiPM top"),
    "sipm_bot_0":  dict(color="#F9E79F", alpha=0.45, edge="#F9E79FBB", lw=0.3, label="SiPM bot"),
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


def load_stl(stl_path, max_faces=800):
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


def discover_volumes(volumes_dir="cad/volumes"):
    d = Path(volumes_dir)
    if not d.exists():
        return []
    return sorted(d.glob("*.stl"))

# ============================================================
# print summary
# ============================================================

def print_summary(interfaces, volume_stls):
    print("\n" + "="*62)
    print("  LEGEND-Theia — Full Geometry Summary")
    print("="*62)
    print(f"\n  Total interfaces : {len(interfaces)}")
    by_surface = Counter(i['surface'] for i in interfaces)
    by_mat     = Counter(i['material_inside'] for i in interfaces)
    for s, n in by_surface.most_common():
        marker = SURFACE_MARKER.get(s, "?")
        print(f"  {marker}  {s:12s} : {n}")
    print()
    for m, n in by_mat.most_common():
        print(f"    {m:35s} : {n}")
    print(f"\n  Volume STLs      : {len(volume_stls)}")
    print("="*62 + "\n")

# ============================================================
# render
# ============================================================

def render(interfaces, volume_stls,
           surface_filter=None, material_filter=None,
           show_volumes=True):

    fig = plt.figure(figsize=(16, 9), facecolor="#0F1923")
    ax3d  = fig.add_subplot(121, projection='3d', facecolor="#0F1923")
    ax_sid = fig.add_subplot(122, facecolor="#0F1923")  # side XZ view

    all_verts = []
    legend_patches = []
    added_surface_labels = set()

    # --------------------------------------------------------
    # 1. volumes
    # --------------------------------------------------------
    if show_volumes:
        vol_order = ["lar_pv", "pen_bege_pv", "pen_icpc_pv",
                     "bege_pv", "icpc_pv", "sipm_top_0", "sipm_bot_0"]
        for vol_name in vol_order:
            stl_path = next((v for v in volume_stls if v.stem == vol_name), None)
            if not stl_path:
                continue
            style = VOLUME_STYLE.get(vol_name,
                dict(color="#CCCCCC", alpha=0.1, edge="#CCCCCC33", lw=0.2, label=vol_name))
            mesh, tris = load_stl(str(stl_path), max_faces=1500)
            if tris is None:
                continue
            all_verts.append(mesh.vertices)

            ax3d.add_collection3d(Poly3DCollection(
                tris, alpha=style['alpha'],
                facecolor=style['color'],
                edgecolor=style['edge'],
                linewidth=style['lw'], zsort='average'
            ))
            for tri in tris[::4]:
                ax_sid.fill(tri[:,0], tri[:,2],
                           color=style['color'], alpha=style['alpha']*2, linewidth=0)

            legend_patches.append(mpatches.Patch(
                facecolor=style['color'], edgecolor="#ffffff22",
                alpha=0.7, label=f"vol: {style['label']}"
            ))

    # --------------------------------------------------------
    # 2. interfaces — grouped and color-coded by surface type
    # --------------------------------------------------------
    filtered = interfaces
    if surface_filter:
        filtered = [i for i in interfaces if i['surface'] == surface_filter]
    if material_filter:
        filtered = [i for i in filtered if i['material_inside'] == material_filter]

    print(f"  Rendering {len(filtered)} interfaces...")

    for iface in filtered:
        surf  = iface['surface']
        color, alpha = SURFACE_COLORS.get(surf, ("#AAAAAA", 0.5))

        mesh, tris = load_stl(iface['stl'], max_faces=300)
        if tris is None:
            continue
        all_verts.append(mesh.vertices)

        ax3d.add_collection3d(Poly3DCollection(
            tris, alpha=alpha,
            facecolor=color,
            edgecolor="none",
            linewidth=0, zsort='average'
        ))
        for tri in tris[::2]:
            ax_sid.fill(tri[:,0], tri[:,2],
                       color=color, alpha=alpha*0.7, linewidth=0)

        if surf not in added_surface_labels:
            marker = SURFACE_MARKER.get(surf, "?")
            cnt = Counter(i['surface'] for i in filtered)[surf]
            legend_patches.append(mpatches.Patch(
                facecolor=color, edgecolor="#ffffff44",
                label=f"{marker} {surf} ({cnt} interfaces)"
            ))
            added_surface_labels.add(surf)

    # --------------------------------------------------------
    # auto-scale
    # --------------------------------------------------------
    if not all_verts:
        print("Nothing to display.")
        return

    vv = np.vstack(all_verts)
    mn, mx = vv.min(axis=0), vv.max(axis=0)
    mid = (mn + mx) / 2
    rng = (mx - mn).max() / 2 * 1.25

    ax3d.set_xlim(mid[0]-rng, mid[0]+rng)
    ax3d.set_ylim(mid[1]-rng, mid[1]+rng)
    ax3d.set_zlim(mid[2]-rng, mid[2]+rng)

    ax_sid.set_xlim(mid[0]-rng, mid[0]+rng)
    ax_sid.set_ylim(mid[2]-rng, mid[2]+rng)
    ax_sid.set_xlabel("X (mm)", color="#8899AA", fontsize=9)
    ax_sid.set_ylabel("Z (mm)", color="#8899AA", fontsize=9)
    ax_sid.tick_params(colors="#8899AA", labelsize=7)
    ax_sid.set_aspect('equal')
    for spine in ax_sid.spines.values():
        spine.set_edgecolor("#1E2D3D")

    for axis in [ax3d.xaxis, ax3d.yaxis, ax3d.zaxis]:
        axis.pane.fill = False
        axis.pane.set_edgecolor("#1E2D3D")
    ax3d.tick_params(colors="#8899AA", labelsize=7)
    ax3d.set_xlabel("X (mm)", color="#8899AA", labelpad=6)
    ax3d.set_ylabel("Y (mm)", color="#8899AA", labelpad=6)
    ax3d.set_zlabel("Z (mm)", color="#8899AA", labelpad=6)
    ax3d.set_title("3D View  (drag to rotate)", color="#E8EDF2", pad=10, fontsize=10)
    ax_sid.set_title("Side view XZ", color="#E8EDF2", pad=10, fontsize=10)

    ax3d.legend(
        handles=legend_patches,
        loc="upper left",
        facecolor="#131F2B",
        edgecolor="#2E5FA3",
        labelcolor="#E8EDF2",
        fontsize=8,
        framealpha=0.92,
    )

    title = f"LEGEND-Theia v0  ·  {len(filtered)} / {len(interfaces)} interfaces"
    if surface_filter:
        title += f"  ·  filter: {surface_filter}"
    if material_filter:
        title += f"  ·  mat: {material_filter}"

    plt.suptitle(title, color="#E8EDF2", fontsize=10, y=1.01)
    plt.tight_layout()
    plt.show()

# ============================================================
# main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="LEGEND-Theia geometry visualizer")
    parser.add_argument("--list",       action="store_true")
    parser.add_argument("--no-volumes", action="store_true")
    parser.add_argument("--surface",    default=None,
                        choices=["blackbody", "specular", "detector", "default"],
                        help="Filter by surface type")
    parser.add_argument("--material",   default=None,
                        help="Filter by material_inside")
    parser.add_argument("--metadata",   default="metadata")
    parser.add_argument("--volumes",    default="cad/volumes")
    args = parser.parse_args()

    interfaces  = load_interfaces(args.metadata)
    volume_stls = [] if args.no_volumes else discover_volumes(args.volumes)

    print_summary(interfaces, volume_stls)

    if args.list:
        return

    render(interfaces, volume_stls,
           surface_filter=args.surface,
           material_filter=args.material,
           show_volumes=not args.no_volumes)


if __name__ == "__main__":
    main()
