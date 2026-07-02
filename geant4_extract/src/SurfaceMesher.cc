#include "SurfaceMesher.hh"

#include "DetectorAssembly.hh"
#include "InterfaceExtractor.hh"
#include "OpticalInterface.hh"

#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <StlAPI_Writer.hxx>

#include <Poly_Triangulation.hxx>

#include <TopExp_Explorer.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <TopLoc_Location.hxx>

#include <TopAbs_ShapeEnum.hxx>

#include <filesystem>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

// Set to true to also dump a .brep per interface (debug only). These
// are not needed for plotting or for Theia and roughly double the file
// count, so the default run skips them.
static constexpr bool kWriteBrep = false;

// ============================================================
// mesh all optical interfaces
// ============================================================

void SurfaceMesher::MeshInterfaces(

    DetectorAssembly &assembly) {

  std::cout << "\nMeshing optical interfaces...\n" << std::endl;

  // --------------------------------------------------------
  // create output directories
  // --------------------------------------------------------

  fs::create_directories("cad/interfaces");
  fs::create_directories("metadata");

  // --------------------------------------------------------
  // loop over interfaces (parallel).
  //
  // Each iteration meshes a DISTINCT iface.boundary and writes
  // its own file — no cross-iteration shared state except the
  // report stream, which is guarded. iface.n_triangles/area_mm2
  // are written per-element, so those are race-free too.
  //
  // schedule(dynamic): interface sizes are very uneven (the
  // LAr-facing boundaries dwarf the small copper patches), so
  // dynamic hand-out keeps all threads busy instead of leaving
  // one thread stuck on the heavy ones under static chunking.
  // --------------------------------------------------------

  const int n = static_cast<int>(assembly.interfaces.size());

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {

    auto &iface = assembly.interfaces[i];

    // ----------------------------------------------------
    // verify boundary topology
    // ----------------------------------------------------

    TopAbs_ShapeEnum shape_type = iface.boundary.ShapeType();

    // ----------------------------------------------------
    // reject volumetric intersections
    // ----------------------------------------------------

    if (shape_type == TopAbs_SOLID || shape_type == TopAbs_COMPSOLID) {

#pragma omp critical(mesher_cout)
      std::cout << "Skipping volumetric intersection:\n"
                << "  interface " << iface.id << std::endl;

      continue;
    }

    // ----------------------------------------------------
    // tessellate OCC boundary surface
    //
    // Standard_False for the parallel flag: the outer OpenMP loop
    // owns the cores. Leaving BRepMesh's own in-parallel flag on
    // would nest parallelism and oversubscribe.
    // ----------------------------------------------------

    BRepMesh_IncrementalMesh mesher(
        iface.boundary,
        0.1,            // linear deflection (relative fraction when isRelative)
        Standard_True,  // isRelative — scale tolerance to each face's size
        0.5,            // angular deflection in radians (~28°)
        Standard_False  // NOT in parallel — outer OpenMP loop owns the cores
    );

    mesher.Perform();

    // ----------------------------------------------------
    // export raw OCC BREP (debug only)
    // ----------------------------------------------------

    std::string brep_name =
        "cad/interfaces/interface_" + std::to_string(iface.id) + ".brep";

    if (kWriteBrep) {
      BRepTools::Write(iface.boundary, brep_name.c_str());
    }

    // ----------------------------------------------------
    // export triangulated STL (binary)
    // ----------------------------------------------------

    std::string stl_name =
        "cad/interfaces/interface_" + std::to_string(iface.id) + ".stl";

    StlAPI_Writer stl_writer;
    stl_writer.ASCIIMode() = Standard_False; // write binary STL

    stl_writer.Write(iface.boundary, stl_name.c_str());

    size_t total_vertices = 0;
    size_t total_triangles = 0;
    size_t total_faces = 0;

    // ----------------------------------------------------
    // iterate surface faces
    // ----------------------------------------------------

    for (TopExp_Explorer explorer(iface.boundary, TopAbs_FACE);
         explorer.More(); explorer.Next()) {

      total_faces++;

      TopoDS_Face face = TopoDS::Face(explorer.Current());

      TopLoc_Location location;

      Handle(Poly_Triangulation) triangulation =
          BRep_Tool::Triangulation(face, location);

      if (triangulation.IsNull())
        continue;

      total_vertices += triangulation->NbNodes();
      total_triangles += triangulation->NbTriangles();
    }

    // ----------------------------------------------------
    // compute surface area (mm²)
    // ----------------------------------------------------

    GProp_GProps props;
    BRepGProp::SurfaceProperties(iface.boundary, props);
    double area_mm2 = props.Mass();

    // ----------------------------------------------------
    // store back onto the interface struct (this element only)
    // ----------------------------------------------------

    iface.n_triangles = (int)total_triangles;
    iface.area_mm2 = area_mm2;

    // ----------------------------------------------------
    // report
    // ----------------------------------------------------

#pragma omp critical(mesher_cout)
    std::cout << "Interface " << iface.id << " : " << iface.materialA << " ↔ "
              << iface.materialB << "\n"
              << "  shape type = " << shape_type << "\n"
              << "  faces      = " << total_faces << "\n"
              << "  vertices   = " << total_vertices << "\n"
              << "  triangles  = " << total_triangles << "\n"
              << "  area_mm2   = " << area_mm2 << "\n"
              << "  brep file  = " << (kWriteBrep ? brep_name : "(skipped)")
              << "\n"
              << "  stl file   = " << stl_name << "\n"
              << std::endl;
  }

  std::cout << "Finished interface meshing.\n" << std::endl;

  // --------------------------------------------------------
  // write interfaces.json now that all stats are populated
  // --------------------------------------------------------

  InterfaceExtractor extractor;
  extractor.WriteInterfacesJSON(assembly, ".");
}