#pragma once

#include <string>

#include <TopoDS_Shape.hxx>

struct OpticalInterface {

  int id;

  int volumeA;
  int volumeB;

  std::string nameA;
  std::string nameB;

  std::string materialA;
  std::string materialB;

  // filled by InterfaceExtractor::Extract
  std::string pv_inside;
  std::string pv_outside;
  std::string lv_inside;
  std::string lv_outside;
  std::string mat_inside;
  std::string mat_outside;
  bool is_detector = false;
  int detector_channel = -1; // GDML channel id, -1 if not a detector

  // normals point from pv_inside toward pv_outside
  TopoDS_Shape boundary;

  // filled by SurfaceMesher
  int n_triangles = 0;
  double area_mm2 = 0.0; // now taken from the triangulation (see below)

  // mesh-quality audit (filled by SurfaceMesher)
  double mesh_area_mm2 = 0.0;  // summed from triangles; source of truth
  double occ_area_mm2 = 0.0;   // analytic GProp area, kept as cross-check
  int n_open_edges = 0;        // edges used by exactly 1 triangle
  int n_nonmanifold_edges = 0; // edges used by >2 triangles
  int n_degenerate_tris = 0;   // triangles with ~zero area
  bool mesh_empty = false;     // tessellation produced no triangles
};
