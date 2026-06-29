#include "VolumeSTLExporter.hh"

#include "DetectorAssembly.hh"
#include "VolumeInstance.hh"

#include <BRepMesh_IncrementalMesh.hxx>
#include <StlAPI_Writer.hxx>
#include <TopoDS_Shape.hxx>

#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ============================================================
// volumes to export by default (key detector volumes)
// add more names here as the geometry grows
// ============================================================

static const std::vector<std::string> kDefaultVolumes = {
    "lar_pv",
    "bege_pv",
    "pen_bege_pv",
    "icpc_pv",
    "pen_icpc_pv",
    "sipm_top_0",
    "sipm_bot_0",
};

// ============================================================
// export
// ============================================================

void VolumeSTLExporter::Export(
    const DetectorAssembly& assembly,
    const std::string& outDir,
    const std::vector<std::string>& names
) {
    fs::create_directories(outDir);

    // build filter set
    const std::vector<std::string>& filter =
        names.empty() ? kDefaultVolumes : names;

    std::set<std::string> wanted(filter.begin(), filter.end());

    int exported = 0;

    for (const auto& vol : assembly.volumes) {

        if (wanted.find(vol.name) == wanted.end())
            continue;

        if (vol.shape.IsNull()) {
            std::cout
                << "VolumeSTLExporter: skipping null shape for "
                << vol.name
                << std::endl;
            continue;
        }

        // ----------------------------------------------------
        // tessellate
        // ----------------------------------------------------

        BRepMesh_IncrementalMesh mesher(vol.shape, 0.1); // 0.1mm deflection
        mesher.Perform();

        // ----------------------------------------------------
        // write STL
        // ----------------------------------------------------

        std::string path =
            outDir + "/" + vol.name + "_" + std::to_string(exported) + ".stl";

        StlAPI_Writer writer;
        writer.Write(vol.shape, path.c_str());

        std::cout
            << "VolumeSTLExporter: wrote "
            << path
            << "  (material: "
            << vol.material
            << ")"
            << std::endl;

        ++exported;
    }

    std::cout
        << "VolumeSTLExporter: exported "
        << exported
        << " / "
        << wanted.size()
        << " requested volumes to "
        << outDir
        << std::endl;
}
