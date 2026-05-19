#include "MaterialExporter.hh"

#include <G4Material.hh>
#include <G4MaterialTable.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>
#include <G4Version.hh>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>

// ============================================================
// canonical name registry
// G4 material name → legend_theia registry key
// ============================================================

static const std::unordered_map<std::string, std::string> kCanonicalNames = {
    {"liquid_argon",           "lar"},
    {"EnrichedGermanium0.076", "germanium_natural"},
    {"NaturalGermanium",       "germanium_natural"},
    {"PEN",                    "pen"},
    {"TPB",                    "tpb"},
    {"BCF91A_core",            "bcf91a_core"},
    {"BCF91A_cladding1",       "bcf91a_cladding_1"},
    {"Copper",                 "copper"},
    {"Tetratex",               "tetratex"},
    {"FusedSilica",            "silica_fused"},
    {"SiPM_photocathode",      "sipm_photocathode"},
};

// ============================================================
// helper structure
// ============================================================

struct SpectralPoint {
    double wavelength_nm;
    double value;
};

// ============================================================
// helper: property units
// ============================================================

static std::string GetUnitForProperty(
    const G4String& name
) {

    if (
        name == "RINDEX"
    ) {
        return "dimensionless";
    }

    if (
        name == "ABSLENGTH" ||
        name == "RAYLEIGH"  ||
        name == "WLSABSLENGTH"
    ) {
        return "mm";
    }

    if (
        name == "FASTCOMPONENT" ||
        name == "SLOWCOMPONENT" ||
        name == "WLSCOMPONENT"
    ) {
        return "relative";
    }

    return "unknown";
}

// ============================================================
// helper: convert units
// ============================================================

static double ConvertPropertyValue(
    const G4String& name,
    double value
) {

    if (
        name == "ABSLENGTH" ||
        name == "RAYLEIGH"  ||
        name == "WLSABSLENGTH"
    ) {
        return value / mm;
    }

    return value;
}

// ============================================================
// export Geant4 optical materials
// ============================================================

void MaterialExporter::Export(
    const char* filename
) {

    std::ofstream out(filename);

    out << std::setprecision(10);

    out << "{\n";

    // --------------------------------------------------------
    // metadata
    // --------------------------------------------------------

    out << "  \"metadata\": {\n";

    out << "    \"geant4_version\": \""
        << G4Version
        << "\",\n";

    out << "    \"wavelength_unit\": \"nm\",\n";

    out << "    \"density_unit\": \"g/cm3\"\n";

    out << "  },\n";

    // --------------------------------------------------------
    // materials
    // --------------------------------------------------------

    out << "  \"materials\": [\n";

    auto* table =
        G4Material::GetMaterialTable();

    if (!table) {

        out << "  ]\n";
        out << "}\n";

        return;
    }

    // --------------------------------------------------------
    // loop over materials
    // --------------------------------------------------------

    for (
        size_t i = 0;
        i < table->size();
        ++i
    ) {

        auto* mat =
            (*table)[i];

        if (!mat)
            continue;

        out << "    {\n";

        // ----------------------------------------------------
        // basic material info
        // ----------------------------------------------------

        out << "      \"name\": \""
            << mat->GetName()
            << "\",\n";

        // ----------------------------------------------------
        // canonical name lookup
        // ----------------------------------------------------

        {
            auto it = kCanonicalNames.find(std::string(mat->GetName()));
            std::string canon = (it != kCanonicalNames.end())
                ? it->second
                : std::string(mat->GetName()); // fallback: use G4 name as-is

            out << "      \"canonical_name\": \""
                << canon
                << "\",\n";
        }

        out << "      \"density_g_cm3\": "
            << mat->GetDensity() / (g / cm3)
            << ",\n";

        out << "      \"state\": "
            << mat->GetState()
            << ",\n";

        // ----------------------------------------------------
        // properties
        // ----------------------------------------------------

        out << "      \"properties\": {\n";

        auto* mpt =
            mat->GetMaterialPropertiesTable();

        bool firstProperty = true;

        if (mpt) {

            // =================================================
            // vector properties
            // =================================================

            std::vector<G4String> names =
                mpt->GetMaterialPropertyNames();

            for (
                size_t j = 0;
                j < names.size();
                ++j
            ) {

                auto name =
                    names[j];

                auto* vec =
                    mpt->GetProperty(name);

                if (!vec)
                    continue;

                // --------------------------------------------
                // collect points
                // --------------------------------------------

                std::vector<SpectralPoint> points;

                for (
                    size_t k = 0;
                    k < vec->GetVectorLength();
                    ++k
                ) {

                    double energy =
                        vec->Energy(k);

                    double wavelength_nm =
                        (h_Planck * c_light)
                        / energy
                        / nm;

                    double value =
                        ConvertPropertyValue(
                            name,
                            (*vec)[k]
                        );

                    points.push_back({
                        wavelength_nm,
                        value
                    });
                }

                // --------------------------------------------
                // sort by wavelength ascending
                // --------------------------------------------

                std::sort(
                    points.begin(),
                    points.end(),
                    [](
                        const SpectralPoint& a,
                        const SpectralPoint& b
                    ) {
                        return
                            a.wavelength_nm <
                            b.wavelength_nm;
                    }
                );

                // --------------------------------------------
                // commas
                // --------------------------------------------

                if (!firstProperty)
                    out << ",\n";

                firstProperty = false;

                // --------------------------------------------
                // export property
                // --------------------------------------------

                out << "        \""
                    << name
                    << "\": {\n";

                out << "          \"unit\": \""
                    << GetUnitForProperty(name)
                    << "\",\n";

                // wavelengths

                out << "          \"wavelength_nm\": [";

                for (
                    size_t k = 0;
                    k < points.size();
                    ++k
                ) {

                    out << points[k].wavelength_nm;

                    if (k + 1 < points.size())
                        out << ", ";
                }

                out << "],\n";

                // values

                out << "          \"values\": [";

                for (
                    size_t k = 0;
                    k < points.size();
                    ++k
                ) {

                    out << points[k].value;

                    if (k + 1 < points.size())
                        out << ", ";
                }

                out << "]\n";

                out << "        }";
            }

            // =================================================
            // constant properties
            // =================================================

            std::vector<G4String> constNames =
                mpt->GetMaterialConstPropertyNames();

            for (
                size_t j = 0;
                j < constNames.size();
                ++j
            ) {

                auto name =
                    constNames[j];

                if (!mpt->ConstPropertyExists(name))
                    continue;

                double value =
                    mpt->GetConstProperty(name);

                if (!firstProperty)
                    out << ",\n";

                firstProperty = false;

                out << "        \""
                    << name
                    << "\": {\n";

                out << "          \"value\": "
                    << value
                    << "\n";

                out << "        }";
            }
        }

        out << "\n";
        out << "      }\n";

        out << "    }";

        if (i + 1 < table->size())
            out << ",";

        out << "\n";
    }

    out << "  ]\n";
    out << "}\n";

    out.close();

    std::cout
        << "Wrote materials:\n"
        << filename
        << std::endl;
}
