#include "VolumeSTLExporter.hh"

#include "DetectorAssembly.hh"
#include "VolumeInstance.hh"

#include <BRepMesh_IncrementalMesh.hxx>
#include <StlAPI_Writer.hxx>
#include <TopoDS_Shape.hxx>

#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <filesystem>
#include <iostream>
#include <map>
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

// A fiber system can place thousands of individually-tiny volumes (every
// fiber's own core/cl1/cl2/coating segment -- LEGEND-1000 alone places
// ~12000 of them), each of which would otherwise get its own STL file.
// Every volume of the same fiber-layer material is bundled into ONE
// combined STL instead; every other component keeps its own per-volume
// file as before.
static bool IsFiberLayerMaterial(const std::string &mat) {
  return mat == "tpb_on_fibers" || mat == "pmma_cl2" || mat == "pmma" ||
         mat == "ps_fibers";
}

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

    // every fiber-layer volume (core/cl1/cl2/tpb, any length/instance)
    // accumulates into this ONE compound instead of being tessellated/
    // written individually below.
    TopoDS_Compound fiber_bundle;
    BRep_Builder fiber_builder;
    fiber_builder.MakeCompound(fiber_bundle);
    bool has_fiber = false;

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

        if (IsFiberLayerMaterial(vol.material)) {
            fiber_builder.Add(fiber_bundle, vol.shape);
            has_fiber = true;
            ++exported;
            continue;
        }

        // ----------------------------------------------------
        // tessellate
        //
        // Relative deflection (isRelative=true) scales the tolerance to
        // each shape's own size, matching SurfaceMesher's interface
        // meshing. A fixed ABSOLUTE 0.1mm deflection is fine for small
        // parts but wildly over-tessellates meter-scale volumes (a
        // legend-1000 cryostat/LAr shell was producing 30MB+ STLs from
        // this alone) -- relative deflection keeps visual fidelity
        // proportional to actual feature size instead of every shape
        // paying the same absolute facet budget regardless of scale.
        // ----------------------------------------------------

        BRepMesh_IncrementalMesh mesher(vol.shape, 0.1, Standard_True, 0.5);
        mesher.Perform();

        // ----------------------------------------------------
        // write STL
        // ----------------------------------------------------

        std::string path =
            outDir + "/" + vol.name + "_" + std::to_string(exported) + ".stl";

        // Binary STL instead of OCC's default ASCII: identical geometry,
        // typically 5-10x smaller (ASCII writes every vertex/normal as
        // human-readable text; binary uses fixed-size records). Matches
        // what SurfaceMesher already does for interface STLs.
        StlAPI_Writer writer;
        writer.ASCIIMode() = Standard_False;
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

    // write the single fiber bundle once, after every matching volume has
    // been accumulated above.
    if (has_fiber) {
        BRepMesh_IncrementalMesh mesher(fiber_bundle, 0.1, Standard_True, 0.5);
        mesher.Perform();

        std::string path = outDir + "/fiber_bundle.stl";

        StlAPI_Writer writer;
        writer.ASCIIMode() = Standard_False;
        writer.Write(fiber_bundle, path.c_str());

        std::cout
            << "VolumeSTLExporter: wrote "
            << path
            << "  (bundled all fiber-layer volumes)"
            << std::endl;
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
