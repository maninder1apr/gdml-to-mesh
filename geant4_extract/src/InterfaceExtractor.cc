#include "InterfaceExtractor.hh"

#include <VolumeInstance.hh>
#include <OpticalInterface.hh>

#include <BRepAlgoAPI_Common.hxx>

#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>

#include <TopoDS_Shape.hxx>

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

using json = nlohmann::json;

// ============================================================
// bbox overlap helper
// ============================================================

static bool BoxesOverlap(

    const TopoDS_Shape& a,
    const TopoDS_Shape& b
) {

    Bnd_Box boxA;
    Bnd_Box boxB;

    BRepBndLib::Add(a, boxA);
    BRepBndLib::Add(b, boxB);

    return !boxA.IsOut(boxB);
}

// ============================================================
// classify a volume that touches LAr
// returns "blackbody" | "specular" | "detector" | "" (skip)
// ============================================================

static std::string ClassifyVolume(const VolumeInstance& vol)
{
    const std::string& mat  = vol.material;
    const std::string& name = vol.name;

    // germanium detectors — blackbody in v0
    if (mat == "EnrichedGermanium0.076" ||
        mat == "NaturalGermanium")
        return "blackbody";

    // PEN encapsulations — specular in v0 (WLS deferred to v1)
    if (mat == "PEN")
        return "specular";

    // TPB fiber coatings — specular in v0 (diffuse/WLS deferred to v1)
    if (mat == "tpb_on_fibers" ||
        name.find("fiber_IB") != std::string::npos ||
        name.find("fiber_coating") != std::string::npos)
        return "specular";

// SiPM photocathodes — detector (silicon only)
    if (mat == "metal_silicon")
        return "detector";

    // SiPM copper wraps — blackbody
    if (mat == "metal_copper")
        return "blackbody";

    // skip everything else (pmma, ps_fibers, copper, world, etc.)
    return "";
}

// ============================================================
// extract interfaces
// ============================================================

void InterfaceExtractor::Extract(

    DetectorAssembly& assembly
) {

    auto& volumes =
        assembly.volumes;

    std::cout
        << "\nChecking interfaces...\n"
        << std::endl;

    for (size_t i = 0;
         i < volumes.size();
         ++i)
    {

        for (size_t j = i + 1;
             j < volumes.size();
             ++j)
        {

            const auto& A = volumes[i];
            const auto& B = volumes[j];

            // ------------------------------------------------
            // one side must be LAr, the other a classified partner
            // ------------------------------------------------

            bool A_is_lar = (A.material == "liquid_argon");
            bool B_is_lar = (B.material == "liquid_argon");

            if (!A_is_lar && !B_is_lar)
                continue;

            const VolumeInstance& other =
                A_is_lar ? B : A;

            std::string hint = ClassifyVolume(other);

            if (hint.empty())
                continue;

            // ------------------------------------------------
            // bbox reject — fast cull
            // ------------------------------------------------

            if (!BoxesOverlap(A.shape, B.shape))
                continue;

            // ------------------------------------------------
            // OCC Boolean intersection
            // ------------------------------------------------

            BRepAlgoAPI_Common common(A.shape, B.shape);

            common.Build();

            if (!common.IsDone())
                continue;

            TopoDS_Shape result = common.Shape();

            if (result.IsNull())
                continue;

            // ------------------------------------------------
            // construct interface object
            // ------------------------------------------------

            OpticalInterface iface;

            iface.id           = assembly.interfaces.size();
            iface.volumeA      = A.id;
            iface.volumeB      = B.id;
            iface.nameA        = A.name;
            iface.nameB        = B.name;
            iface.materialA    = A.material;
            iface.materialB    = B.material;
            iface.boundary     = result;
            iface.surface_hint = hint;

            assembly.interfaces.push_back(iface);

            // ------------------------------------------------
            // debug output
            // ------------------------------------------------

            std::cout
                << "INTERFACE " << iface.id << ":\n  "
                << A.name << " (" << A.material << ")\n  <-> "
                << B.name << " (" << B.material << ")\n"
                << "  surface: " << hint << "\n"
                << std::endl;
        }
    }

    std::cout
        << "\nTotal interfaces found: "
        << assembly.interfaces.size()
        << std::endl;
}

// ============================================================
// write metadata/interfaces.json
// called after SurfaceMesher has populated n_triangles/area_mm2
// ============================================================

void InterfaceExtractor::WriteInterfacesJSON(
    const DetectorAssembly& assembly,
    const std::string& outDir)
{
    json arr = json::array();

    for (const auto& iface : assembly.interfaces) {

        // orient: lv_inside = non-LAr, lv_outside = LAr
        bool A_is_lar = (iface.materialA == "liquid_argon");

        std::string lv_inside    = A_is_lar ? iface.nameB     : iface.nameA;
        std::string lv_outside   = A_is_lar ? iface.nameA     : iface.nameB;
        std::string mat_inside   = A_is_lar ? iface.materialB : iface.materialA;
        std::string mat_outside  = A_is_lar ? iface.materialA : iface.materialB;

        // use hint computed during extraction
        std::string surface = iface.surface_hint;
        json        det_id  = nullptr;

        if (surface == "detector")
            det_id = iface.id;

        json entry;
        entry["id"]               = iface.id;
        entry["stl"]              = "cad/interfaces/interface_"
                                    + std::to_string(iface.id) + ".stl";
        entry["lv_inside"]        = lv_inside;
        entry["lv_outside"]       = lv_outside;
        entry["material_inside"]  = mat_inside;
        entry["material_outside"] = mat_outside;
        entry["surface"]          = surface;
        entry["detector_id"]      = det_id;
        entry["n_triangles"]      = iface.n_triangles;
        entry["area_mm2"]         = iface.area_mm2;

        arr.push_back(entry);
    }

    std::string path = outDir + "/metadata/interfaces.json";
    std::ofstream f(path);

    if (!f) {
        std::cerr << "ERROR: could not write " << path << std::endl;
        return;
    }

    f << std::setw(2) << arr << std::endl;
    std::cout << "Wrote " << path << std::endl;
}
