#!/usr/bin/env python3
"""
visualize_pen.py — focused view of PEN encapsulation border surface with cross-section.

Usage:
    python3 visualize_pen.py            # BeGe PEN encapsulation
    python3 visualize_pen.py --icpc     # ICPC PEN encapsulation
    python3 visualize_pen.py --clip     # show cross-section (half cut)
    python3 visualize_pen.py --icpc --clip
"""

import argparse
import json
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import trimesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

BORDER_COLOR = "#3498DB"   # blue   — LAr/PEN border surface
SKIN_COLOR   = "#F39C12"   # orange — PEN skin surface

CONFIGS = {
    "bege": {
        "det_vol":   "bege_pv",
        "pen_vol":   "pen_bege_pv",
        "pen_lv":    "pen_bege_lv",
        "det_label": "BeGe (Germanium)",
        "pen_label": "PEN (BeGe encapsulation)",
        "det_color": "#A9DFBF",
        "pen_color": "#D7BDE2",
    },
    "icpc": {
        "det_vol":   "icpc_pv",
        "pen_vol":   "pen_icpc_pv",
        "pen_lv":    "pen_icpc_lv",
        "det_label": "ICPC (Germanium)",
        "pen_label": "PEN (ICPC encapsulation)",
        "det_color": "#A9DFBF",
        "pen_color": "#D7BDE2",
    },
}


def load_mesh(stl_path, max_faces=5000):
    p = Path(stl_path)
    if not p.exists():
        print(f"  WARNING: {p} not found")
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


def cap_mesh(mesh, plane_normal=(1, 0, 0), plane_origin=(0, 0, 0)):
    """Slice mesh at a plane and cap the cut face for a clean cross-section."""
    try:
        sliced = trimesh.intersections.slice_mesh_plane(
            mesh,
            plane_normal=np.array(plane_normal),
            plane_origin=np.array(plane_origin),
            cap=True,
        )
        return sliced
    except Exception as e:
        print(f"  WARNING: mesh slicing failed: {e}")
        return mesh


def mesh_to_tris(mesh, max_faces=6000):
    faces = mesh.faces
    if len(faces) > max_faces:
        idx = np.random.choice(len(faces), max_faces, replace=False)
        return mesh.vertices[faces[idx]]
    return mesh.vertices[faces]


def add_solid(ax3d, ax_top, ax_sid, tris, color, alpha,
              edge="none", lw=0, step3d=1, step2d=3):
    if tris is None or len(tris) == 0:
        return
    ax3d.add_collection3d(Poly3DCollection(
        tris[::step3d], alpha=alpha,
        facecolor=color, edgecolor=edge,
        linewidth=lw, zsort='average'
    ))
    for tri in tris[::step2d]:
        ax_top.fill(tri[:, 0], tri[:, 1],
                    color=color, alpha=min(alpha * 1.8, 1.0), linewidth=0)
        ax_sid.fill(tri[:, 0], tri[:, 2],
                    color=color, alpha=min(alpha * 1.8, 1.0), linewidth=0)


def add_skin(ax3d, ax_top, ax_sid, tris, color, step3d=1, step2d=4):
    if tris is None or len(tris) == 0:
        return
    ax3d.add_collection3d(Poly3DCollection(
        tris[::step3d], alpha=0.15,
        facecolor=color, edgecolor=color,
        linewidth=0.3, zsort='average'
    ))
    for tri in tris[::step2d]:
        ax_top.plot(np.append(tri[:, 0], tri[0, 0]),
                    np.append(tri[:, 1], tri[0, 1]),
                    color=color, lw=0.5, alpha=0.5)
        ax_sid.plot(np.append(tri[:, 0], tri[0, 0]),
                    np.append(tri[:, 2], tri[0, 2]),
                    color=color, lw=0.5, alpha=0.5)


def render(detector="bege", clip=False):
    cfg = CONFIGS[detector]

    interfaces = json.load(open("metadata/interfaces.json"))
    surfaces   = json.load(open("metadata/surfaces.json"))

    pen_border = next((i for i in interfaces
                       if i['material_inside'] == 'PEN'
                       and cfg['pen_vol'] in i['lv_inside']), None)

    pen_skins = [s for s in surfaces
                 if s['type'] == 'skin' and cfg['pen_lv'] in s['lv_skin']]

    if clip:
        fig = plt.figure(figsize=(10, 9), facecolor="#0F1923")
        ax3d   = fig.add_subplot(111, projection='3d', facecolor="#0F1923")
        ax_top = None
        ax_sid = None
    else:
        fig = plt.figure(figsize=(17, 9), facecolor="#0F1923")
        ax3d   = fig.add_subplot(131, projection='3d', facecolor="#0F1923")
        ax_top = fig.add_subplot(132, facecolor="#0F1923")
        ax_sid = fig.add_subplot(133, facecolor="#0F1923")

    if clip:
        ax3d.view_init(elev=0, azim=90)
    else:
        ax3d.view_init(elev=20, azim=-60)

    legend_patches = []
    all_verts = []

    # 1. detector volume — commented out
    # mesh_d, tris_d = load_mesh(f"cad/volumes/{cfg['det_vol']}.stl", max_faces=3000)
    # if tris_d is not None and len(tris_d) > 0:
    #     all_verts.append(mesh_d.vertices)
    #     add_solid(ax3d, ax_top, ax_sid, tris_d,
    #               cfg['det_color'], 0.40,
    #               edge="#A9DFBF88", lw=0.3, step3d=1, step2d=4)
    #     legend_patches.append(mpatches.Patch(
    #         facecolor=cfg['det_color'], edgecolor="#ffffff22",
    #         alpha=0.7, label=f"vol: {cfg['det_label']}"
    #     ))

    # 2. PEN volume — commented out
    # mesh_p, tris_p = load_mesh(f"cad/volumes/{cfg['pen_vol']}.stl", max_faces=3000)
    # if tris_p is not None and len(tris_p) > 0:
    #     all_verts.append(mesh_p.vertices)
    #     add_solid(ax3d, ax_top, ax_sid, tris_p,
    #               cfg['pen_color'], 0.06,
    #               edge="none", lw=0, step3d=1, step2d=4)
    #     legend_patches.append(mpatches.Patch(
    #         facecolor=cfg['pen_color'], edgecolor="#ffffff22",
    #         alpha=0.4, label=f"vol: {cfg['pen_label']} (faint)"
    #     ))

    # 3. PEN border surface — solid blue
    if pen_border:
        mesh_b, _ = load_mesh(pen_border['stl'], max_faces=6000)
        if mesh_b is not None:
            if clip:
                mesh_b = cap_mesh(mesh_b)
            tris_b = mesh_to_tris(mesh_b, max_faces=6000)
            if tris_b is not None and len(tris_b) > 0:
                all_verts.append(mesh_b.vertices)
                ax3d.add_collection3d(Poly3DCollection(
                    tris_b[::1], alpha=0.85,
                    facecolor=BORDER_COLOR, edgecolor="#00000022",
                    linewidth=0.1, zsort='average'
                ))
                if not clip and ax_top is not None and ax_sid is not None:
                    for tri in tris_b[::2]:
                        ax_top.fill(tri[:, 0], tri[:, 1],
                                    color=BORDER_COLOR, alpha=0.85, linewidth=0)
                        ax_sid.fill(tri[:, 0], tri[:, 2],
                                    color=BORDER_COLOR, alpha=0.85, linewidth=0)
                legend_patches.append(mpatches.Patch(
                    facecolor=BORDER_COLOR, edgecolor="#ffffff44",
                    label=f"border: {pen_border['lv_inside']} | {pen_border['lv_outside']}\n"
                          f"  {pen_border['surface']}  ·  "
                          f"{pen_border['n_triangles']:,} tris  ·  "
                          f"{pen_border['area_mm2']:.0f} mm2"
                ))

    # 4. PEN skin surface — only in non-clip mode
    if pen_skins and not clip:
        mesh_s, tris_s = load_mesh(f"cad/volumes/{cfg['pen_vol']}.stl", max_faces=2000)
        if tris_s is not None and len(tris_s) > 0:
            skin = pen_skins[0]
            add_skin(ax3d, ax_top, ax_sid, tris_s, SKIN_COLOR,
                     step3d=1, step2d=4)
            legend_patches.append(mpatches.Patch(
                facecolor=SKIN_COLOR, edgecolor=SKIN_COLOR,
                alpha=0.5,
                label=f"skin: {skin['lv_skin']}\n"
                      f"  {skin['finish']}  ·  {skin['surf_type']}"
            ))

    # auto-scale
    if all_verts:
        vv = np.vstack(all_verts)
        mn_v, mx_v = vv.min(axis=0), vv.max(axis=0)
        mid = (mn_v + mx_v) / 2
        rng = (mx_v - mn_v).max() / 2 * 1.3

        ax3d.set_xlim(0 if clip else mid[0] - rng, mid[0] + rng)
        ax3d.set_ylim(mid[1] - rng, mid[1] + rng)
        ax3d.set_zlim(mid[2] - rng, mid[2] + rng)

        if not clip and ax_top is not None and ax_sid is not None:
            for ax2, xi, yi, xl, yl in [
                (ax_top, 0, 1, "X (mm)", "Y (mm)"),
                (ax_sid, 0, 2, "X (mm)", "Z (mm)"),
            ]:
                ax2.set_xlim(mid[xi] - rng, mid[xi] + rng)
                ax2.set_ylim(mid[yi] - rng, mid[yi] + rng)
                ax2.set_xlabel(xl, color="#8899AA", fontsize=9)
                ax2.set_ylabel(yl, color="#8899AA", fontsize=9)
                ax2.tick_params(colors="#8899AA", labelsize=7)
                ax2.set_aspect('equal')
                for spine in ax2.spines.values():
                    spine.set_edgecolor("#1E2D3D")

    for axis in [ax3d.xaxis, ax3d.yaxis, ax3d.zaxis]:
        axis.pane.fill = False
        axis.pane.set_edgecolor("#1E2D3D")
    ax3d.tick_params(colors="#8899AA", labelsize=7)
    ax3d.set_xlabel("X (mm)", color="#8899AA", labelpad=6)
    ax3d.set_ylabel("Y (mm)", color="#8899AA", labelpad=6)
    ax3d.set_zlabel("Z (mm)", color="#8899AA", labelpad=6)
    ax3d.set_title("3D View (clipped)" if clip else "3D View",
                   color="#E8EDF2", pad=10, fontsize=10)

    if not clip and ax_top is not None and ax_sid is not None:
        ax_top.set_title("Top view (XY)", color="#E8EDF2", pad=10, fontsize=10)
        ax_sid.set_title("Side view (XZ)", color="#E8EDF2", pad=10, fontsize=10)

    ax3d.legend(
        handles=legend_patches,
        loc="upper left",
        facecolor="#131F2B",
        edgecolor="#2E5FA3",
        labelcolor="#E8EDF2",
        fontsize=8,
        framealpha=0.92,
    )

    det_name = "BeGe" if detector == "bege" else "ICPC"
    clip_note = "  [cross-section x>=0]" if clip else ""
    plt.suptitle(
        f"PEN Encapsulation ({det_name})  -  "
        f"Border surface (solid blue)  |  Skin surface (orange){clip_note}",
        color="#E8EDF2", fontsize=10, y=1.01
    )
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="PEN surface visualizer")
    parser.add_argument("--icpc", action="store_true",
                        help="Show ICPC instead of BeGe")
    parser.add_argument("--clip", action="store_true",
                        help="Show cross-section (half cut at x=0)")
    args = parser.parse_args()
    render("icpc" if args.icpc else "bege", clip=args.clip)


if __name__ == "__main__":
    main()