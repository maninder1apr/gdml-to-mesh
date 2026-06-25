#include "occ_mesher.hh"

#include "SolidConverter.hh"
#include "AssemblyBuilder.hh"
#include "MaterialExporter.hh"
#include "SurfaceExporter.hh"
#include "InterfaceExtractor.hh"
#include "SurfaceMesher.hh"
#include "VolumeSTLExporter.hh"

#include <G4GDMLParser.hh>
#include <G4GDMLAuxStructType.hh>
#include <G4GeometryTolerance.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SystemOfUnits.hh>
#include <G4VSolid.hh>

#include <algorithm>

#include <map>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>

#include <BRepTools.hxx>

#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

// ============================================================
// load GDML + convert solids
// ============================================================

void OCCMesher::run(
    const std::string& gdml_file
) {

    fs::path gdml_path =
        fs::absolute(gdml_file);

    if (!fs::exists(gdml_path)) {

        std::cerr
            << "GDML file not found:\n"
            << gdml_path
            << std::endl;

        return;
    }

    // ------------------------------------------------------------
    // IMPORTANT:
    // switch cwd so relative includes work
    // ------------------------------------------------------------

    fs::path old_cwd  = fs::current_path();
    fs::path gdml_dir = gdml_path.parent_path();

    fs::current_path(gdml_dir);

    G4GDMLParser parser;

    parser.Read(

        gdml_path.filename().string(),

        false
    );

    fs::current_path(old_cwd);

    // ------------------------------------------------------------
    // read optical detector map from GDML <userinfo>
    // <auxiliary auxtype="RMG_detector" auxvalue="optical">
    //   <auxiliary auxtype="sipm_top_0" auxvalue="1000"/>
    //   ...
    // </auxiliary>
    // ------------------------------------------------------------

    std::map<std::string, int> optical_detectors;

    const G4GDMLAuxListType* auxList = parser.GetAuxList();

    if (auxList) {
        for (const auto& aux : *auxList) {
            if (aux.type == "RMG_detector" && aux.value == "optical") {
                if (aux.auxList) {
                    for (const auto& child : *aux.auxList)
                        optical_detectors[child.type] =
                            std::stoi(child.value);
                }
            }
        }
    }

    std::cout
        << "Optical detectors from GDML: "
        << optical_detectors.size()
        << std::endl;

    // ------------------------------------------------------------
    // create output directories
    // ------------------------------------------------------------

    fs::create_directories("metadata");
    fs::create_directories("cad");

    // ------------------------------------------------------------
    // export materials
    // ------------------------------------------------------------

    MaterialExporter materialExporter;

    materialExporter.Export(
        "metadata/materials.json"
    );

    // ------------------------------------------------------------
    // export optical surfaces
    // (must happen after GDML parse, while G4 surface tables live)
    // ------------------------------------------------------------

    SurfaceExporter surfaceExporter;

    surfaceExporter.Export(".");

    // ------------------------------------------------------------
    // build full detector assembly
    // ------------------------------------------------------------

    auto* world =
        parser.GetWorldVolume();

    if (!world) {

        std::cerr
            << "No world volume"
            << std::endl;

        return;
    }

    AssemblyBuilder builder;

    DetectorAssembly detector =

        builder.Build(world);

    // ------------------------------------------------------------
    // export key volumes as STL for the Python visualizer
    // ------------------------------------------------------------

    VolumeSTLExporter volumeExporter;

    volumeExporter.Export(detector, "cad/volumes");

    // ------------------------------------------------------------
    // extract touching interfaces
    // ------------------------------------------------------------

    // Fuzzy tolerance for coincident-face detection between touching
    // sibling volumes. Geant4's surface tolerance is world-extent
    // relative and can be unrealistically small/large, so floor it at
    // 0.1 µm
    double surf_tol_mm =
        G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / CLHEP::mm;

    double fuzzy_mm = std::max(surf_tol_mm, 1e-4);

    std::cout << "\nG4 surface tolerance = " << surf_tol_mm
              << " mm; using fuzzy = " << fuzzy_mm << " mm\n";

    InterfaceExtractor extractor;

    extractor.Extract(
        detector,
        optical_detectors,
        fuzzy_mm
    );

    // ------------------------------------------------------------
    // mesh optical interfaces
    // ------------------------------------------------------------

    SurfaceMesher mesher;

    mesher.MeshInterfaces(
        detector
    );

    // ------------------------------------------------------------
    // write OCC assembly
    // ------------------------------------------------------------

    BRepTools::Write(

        detector.assembly,

        "cad/detector.brep"
    );

    std::cout
        << "\nWrote detector assembly:\n"
        << "cad/detector.brep"
        << std::endl;

    // ------------------------------------------------------------
    // statistics
    // ------------------------------------------------------------

    std::cout
        << "\nDetector summary\n"
        << "-----------------------------\n"
        << "Volumes: "
        << detector.volumes.size()
        << "\nInterfaces: "
        << detector.interfaces.size()
        << std::endl;
}
