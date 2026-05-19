#include "SurfaceMesher.hh"

#include "DetectorAssembly.hh"
#include "OpticalInterface.hh"
#include "InterfaceExtractor.hh"

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

namespace fs = std::filesystem;

// ============================================================
// mesh all optical interfaces
// ============================================================

void SurfaceMesher::MeshInterfaces(

    DetectorAssembly& assembly
) {

    std::cout
        << "\nMeshing optical interfaces...\n"
        << std::endl;

    // --------------------------------------------------------
    // create output directories
    // --------------------------------------------------------

    fs::create_directories("cad/interfaces");
    fs::create_directories("metadata");

    // --------------------------------------------------------
    // loop over interfaces
    // --------------------------------------------------------

    for (auto& iface : assembly.interfaces)
    {

        // ----------------------------------------------------
        // verify boundary topology
        // ----------------------------------------------------

        TopAbs_ShapeEnum shape_type =

            iface.boundary.ShapeType();

        // ----------------------------------------------------
        // reject volumetric intersections
        // ----------------------------------------------------

        if (
            shape_type == TopAbs_SOLID ||
            shape_type == TopAbs_COMPSOLID
        )
        {

            std::cout
                << "Skipping volumetric intersection:\n"
                << "  interface "
                << iface.id
                << std::endl;

            continue;
        }

        // ----------------------------------------------------
        // tessellate OCC boundary surface
        // ----------------------------------------------------

        BRepMesh_IncrementalMesh mesher(

            iface.boundary,

            0.01   // mesh resolution in mm
        );

        mesher.Perform();

        // ----------------------------------------------------
        // export raw OCC BREP
        // ----------------------------------------------------

        std::string brep_name =

            "cad/interfaces/interface_" +
            std::to_string(iface.id) +
            ".brep";

        BRepTools::Write(

            iface.boundary,

            brep_name.c_str()
        );

        // ----------------------------------------------------
        // export triangulated STL
        // ----------------------------------------------------

        std::string stl_name =

            "cad/interfaces/interface_" +
            std::to_string(iface.id) +
            ".stl";

        StlAPI_Writer stl_writer;

        stl_writer.Write(

            iface.boundary,

            stl_name.c_str()
        );

        size_t total_vertices  = 0;
        size_t total_triangles = 0;
        size_t total_faces     = 0;

        // ----------------------------------------------------
        // iterate surface faces
        // ----------------------------------------------------

        for (
            TopExp_Explorer explorer(

                iface.boundary,
                TopAbs_FACE
            );

            explorer.More();

            explorer.Next()
        )
        {

            total_faces++;

            TopoDS_Face face =

                TopoDS::Face(
                    explorer.Current()
                );

            TopLoc_Location location;

            Handle(Poly_Triangulation) triangulation =

                BRep_Tool::Triangulation(

                    face,
                    location
                );

            if (triangulation.IsNull())
                continue;

            total_vertices +=
                triangulation->NbNodes();

            total_triangles +=
                triangulation->NbTriangles();
        }

        // ----------------------------------------------------
        // compute surface area (mm²)
        // ----------------------------------------------------

        GProp_GProps props;
        BRepGProp::SurfaceProperties(iface.boundary, props);
        double area_mm2 = props.Mass();

        // ----------------------------------------------------
        // store back onto the interface struct
        // ----------------------------------------------------

        iface.n_triangles = (int)total_triangles;
        iface.area_mm2    = area_mm2;

        // ----------------------------------------------------
        // report
        // ----------------------------------------------------

        std::cout
            << "Interface "
            << iface.id
            << " : "
            << iface.materialA
            << " ↔ "
            << iface.materialB
            << "\n"
            << "  shape type = "
            << shape_type
            << "\n"
            << "  faces      = "
            << total_faces
            << "\n"
            << "  vertices   = "
            << total_vertices
            << "\n"
            << "  triangles  = "
            << total_triangles
            << "\n"
            << "  area_mm2   = "
            << area_mm2
            << "\n"
            << "  brep file  = "
            << brep_name
            << "\n"
            << "  stl file   = "
            << stl_name
            << "\n"
            << std::endl;
    }

    std::cout
        << "Finished interface meshing.\n"
        << std::endl;

    // --------------------------------------------------------
    // write interfaces.json now that all stats are populated
    // --------------------------------------------------------

    InterfaceExtractor extractor;
    extractor.WriteInterfacesJSON(assembly, ".");
}
