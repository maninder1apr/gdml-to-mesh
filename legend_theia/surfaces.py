"""
legend_theia/surfaces.py

Builds surface lookup from metadata/surfaces.json.

Matching mirrors Geant4 semantics:
  * border surfaces (G4LogicalBorderSurface) are defined between two *physical*
    volumes  -> matched via pv_pair  <-> interface pv_inside/pv_outside
  * skin surfaces (G4LogicalSkinSurface) wrap a whole *logical* volume
    -> matched via lv_skin  <-> interface lv_inside/lv_outside
"""

from __future__ import annotations

import json
from pathlib import Path


def load_surfaces(surfaces_json: str | Path) -> dict[str, dict]:
    """
    Load surfaces.json and build lookup indices.

    Parameters
    ----------
    surfaces_json
        Path to metadata/surfaces.json produced by gdml-to-mesh.

    Returns
    -------
    Dict with two sub-indices:
      ``{"border": {(pv_a, pv_b): surf}, "skin": {lv: surf}}``
    """
    with open(surfaces_json) as f:
        surfaces = json.load(f)

    border: dict[tuple[str, str], dict] = {}
    skin:   dict[str, dict]             = {}

    for surf in surfaces:
        if surf.get("type") == "border":
            pair = surf.get("pv_pair", [])
            if len(pair) == 2:
                border[(pair[0], pair[1])] = surf

        elif surf.get("type") == "skin":
            lv = surf.get("lv_skin")
            if isinstance(lv, str) and lv:
                skin[lv] = surf

    return {"border": border, "skin": skin}


def get_surface_for_interface(
    interface: dict,
    surface_index: dict[str, dict],
) -> dict | None:
    """
    Look up the G4 surface definition for a given interface entry.

    Border surfaces take precedence over skin surfaces, matching Geant4's
    surface resolution order.

    Parameters
    ----------
    interface
        One entry from interfaces.json.
    surface_index
        Index returned by load_surfaces().

    Returns
    -------
    Surface dict or None if no surface found.
    """
    border = surface_index.get("border", {})
    skin   = surface_index.get("skin", {})

    # 1. border surface on physical volumes
    pv_in  = interface.get("pv_inside", "")
    pv_out = interface.get("pv_outside", "")
    surf = border.get((pv_in, pv_out)) or border.get((pv_out, pv_in))
    if surf is not None:
        return surf

    # 2. skin surface, matched on either logical volume
    lv_in  = interface.get("lv_inside", "")
    lv_out = interface.get("lv_outside", "")
    return skin.get(lv_in) or skin.get(lv_out)
