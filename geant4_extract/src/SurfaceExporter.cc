#include "SurfaceExporter.hh"

#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4VPhysicalVolume.hh>
#include <G4LogicalVolume.hh>

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using json = nlohmann::json;

// ============================================================
// helper: dump an optical MPT to JSON
// returns { "EFFICIENCY": { "wavelength_nm": [...], "values": [...] }, ... }
// ============================================================

static json ExportMPT(G4MaterialPropertiesTable* mpt)
{
    json props = json::object();

    if (!mpt)
        return props;

    // --------------------------------------------------------
    // vector properties (wavelength-indexed)
    // --------------------------------------------------------

    for (auto& name : mpt->GetMaterialPropertyNames()) {

        auto* vec = mpt->GetProperty(name);

        if (!vec)
            continue;

        std::vector<double> wls;
        std::vector<double> vals;

        for (size_t k = 0; k < vec->GetVectorLength(); ++k) {

            double energy = vec->Energy(k);
            double wl_nm  = (h_Planck * c_light) / energy / nm;

            wls.push_back(wl_nm);
            vals.push_back((*vec)[k]);
        }

        props[std::string(name)] = {
            {"wavelength_nm", wls},
            {"values",        vals}
        };
    }

    // --------------------------------------------------------
    // constant properties
    // --------------------------------------------------------

    for (auto& name : mpt->GetMaterialConstPropertyNames()) {

        if (!mpt->ConstPropertyExists(name))
            continue;

        double value = mpt->GetConstProperty(name);

        props[std::string(name)] = {{"value", value}};
    }

    return props;
}

// ============================================================
// helper: G4 surface model/finish/type → string
// ============================================================

static std::string ModelStr(G4OpticalSurface* s)
{
    if (!s) return "unknown";
    switch (s->GetModel()) {
        case glisur:   return "glisur";
        case unified:  return "unified";
        case LUT:      return "LUT";
        default:       return "unknown";
    }
}

static std::string FinishStr(G4OpticalSurface* s)
{
    if (!s) return "unknown";
    switch (s->GetFinish()) {
        case polished:              return "polished";
        case polishedfrontpainted:  return "polishedfrontpainted";
        case polishedbackpainted:   return "polishedbackpainted";
        case ground:                return "ground";
        case groundfrontpainted:    return "groundfrontpainted";
        case groundbackpainted:     return "groundbackpainted";
        default:                    return "unknown";
    }
}

static std::string TypeStr(G4OpticalSurface* s)
{
    if (!s) return "unknown";
    switch (s->GetType()) {
        case dielectric_metal:      return "dielectric_metal";
        case dielectric_dielectric: return "dielectric_dielectric";
        case dielectric_LUT:        return "dielectric_LUT";
        case firsov:                return "firsov";
        case x_ray:                 return "x_ray";
        default:                    return "unknown";
    }
}

// ============================================================
// export metadata/surfaces.json
// ============================================================

void SurfaceExporter::Export(const std::string& outDir)
{
    json arr = json::array();

    // --------------------------------------------------------
    // G4LogicalBorderSurface — between two physical volumes
    // --------------------------------------------------------

    auto* borderTable =
        G4LogicalBorderSurface::GetSurfaceTable();

    if (borderTable) {

        for (auto& kv : *borderTable) {

            auto* surf = kv.second;

            if (!surf)
                continue;

            auto* vol1 = surf->GetVolume1();
            auto* vol2 = surf->GetVolume2();

            if (!vol1 || !vol2)
                continue;

            auto* osurf = dynamic_cast<G4OpticalSurface*>(
                surf->GetSurfaceProperty()
            );

            json entry;

            entry["type"]    = "border";
            entry["g4_name"] = std::string(surf->GetName());

            entry["pv_pair"] = {
                std::string(vol1->GetName()),
                std::string(vol2->GetName())
            };

            entry["lv_pair"] = {
                std::string(vol1->GetLogicalVolume()->GetName()),
                std::string(vol2->GetLogicalVolume()->GetName())
            };

            entry["model"]  = ModelStr(osurf);
            entry["finish"] = FinishStr(osurf);
            entry["surf_type"] = TypeStr(osurf);

            entry["MPT"] = osurf
                ? ExportMPT(osurf->GetMaterialPropertiesTable())
                : json::object();

            arr.push_back(entry);
        }
    }

    // --------------------------------------------------------
    // G4LogicalSkinSurface — wraps an entire logical volume
    // --------------------------------------------------------

    auto* skinTable =
        G4LogicalSkinSurface::GetSurfaceTable();

    if (skinTable) {

        for (auto& kv : *skinTable) {

            auto* surf = kv.second;

            if (!surf)
                continue;

            auto* lv = surf->GetLogicalVolume();

            if (!lv)
                continue;

            auto* osurf = dynamic_cast<G4OpticalSurface*>(
                surf->GetSurfaceProperty()
            );

            json entry;

            entry["type"]    = "skin";
            entry["g4_name"] = std::string(surf->GetName());
            entry["lv_skin"] = std::string(lv->GetName());

            entry["model"]     = ModelStr(osurf);
            entry["finish"]    = FinishStr(osurf);
            entry["surf_type"] = TypeStr(osurf);

            entry["MPT"] = osurf
                ? ExportMPT(osurf->GetMaterialPropertiesTable())
                : json::object();

            arr.push_back(entry);
        }
    }

    // --------------------------------------------------------
    // write file
    // --------------------------------------------------------

    std::string path = outDir + "/metadata/surfaces.json";
    std::ofstream f(path);

    if (!f) {
        std::cerr
            << "ERROR: could not write "
            << path
            << std::endl;
        return;
    }

    f << std::setw(2) << arr << std::endl;

    std::cout
        << "Wrote "
        << path
        << " ("
        << arr.size()
        << " surfaces)"
        << std::endl;
}
