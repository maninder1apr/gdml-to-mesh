#include <G4GDMLParser.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4VSolid.hh>
#include <G4Polyhedron.hh>
#include <G4SystemOfUnits.hh>

#include "cnpy/cnpy.h"

#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Usage: extract_mesh <lar|bege|icpc|list>\n";
        return 1;
    }

    std::string mode = argv[1];
    if (mode != "lar" &&
        mode != "bege" &&
        mode != "icpc" &&
        mode != "list") {
        std::cerr << "Invalid mode: " << mode << "\n";
        std::cerr << "Usage: extract_mesh <lar|bege|icpc|list>\n";
        return 1;
    }

    // ------------------------------------------------------------
    // read GDML
    // ------------------------------------------------------------
    G4GDMLParser parser;
    parser.Read("../gdml/scarf_pen.gdml");

    auto lvStore = G4LogicalVolumeStore::GetInstance();
    if (!lvStore) {
        std::cerr << "ERROR: No logical volume store found\n";
        return 1;
    }

    // ------------------------------------------------------------
    // LIST MODE (diagnostic only)
    // ------------------------------------------------------------
    if (mode == "list") {
        std::unordered_set<std::string> seen;
        std::cout << "\n=== Solid names in GDML ===\n";
        for (auto lv : *lvStore) {
            if (!lv || !lv->GetSolid()) continue;
            const std::string& name = lv->GetSolid()->GetName();
            if (seen.insert(name).second) {
                std::cout << name << "\n";
            }
        }
        std::cout << "==========================\n";
        return 0;
    }

    std::string outdir = "../meshes/" + mode + "/";

    std::unordered_set<std::string> written_solids;
    size_t n_written = 0;
    size_t n_skipped_poly = 0;
    size_t n_skipped_filter = 0;

    // ------------------------------------------------------------
    // iterate logical volumes
    // ------------------------------------------------------------
    for (auto lv : *lvStore) {
        if (!lv) continue;

        G4VSolid* solid = lv->GetSolid();
        if (!solid) continue;

        const std::string& name = solid->GetName();

        // --------------------------------------------------------
        // exact solid selection (THIS IS THE FIX)
        // --------------------------------------------------------
        bool accept = false;

        if (mode == "lar") {
            accept = (name == "lar_s");
        }
        else if (mode == "bege") {
            accept = (name == "bege_lv");
        }
        else if (mode == "icpc") {
            accept = (name == "icpc_lv");
        }

        if (!accept) {
            ++n_skipped_filter;
            continue;
        }

        // ---- avoid duplicates
        if (written_solids.count(name)) continue;
        written_solids.insert(name);

        std::cout << "Extracting [" << mode << "] solid: "
                  << name << std::endl;

        // ------------------------------------------------------------
        // create polyhedron
        // ------------------------------------------------------------
        G4Polyhedron* poly = solid->CreatePolyhedron();
        if (!poly) {
            std::cerr << "  -> no polyhedron (boolean or unsupported solid), skipping\n";
            ++n_skipped_poly;
            continue;
        }

        std::vector<float> vertices;
        std::vector<int> indices;

        // ------------------------------------------------------------
        // vertices (Geant4 uses 1-based indexing)
        // ------------------------------------------------------------
        int nVert = poly->GetNoVertices();
        vertices.reserve(nVert * 3);

        for (int i = 1; i <= nVert; ++i) {
            auto v = poly->GetVertex(i);
            vertices.push_back(v.x() / mm);
            vertices.push_back(v.y() / mm);
            vertices.push_back(v.z() / mm);
        }

        // ------------------------------------------------------------
        // facets
        // ------------------------------------------------------------
        int nFace = poly->GetNoFacets();
        indices.reserve(nFace * 3);

        for (int i = 1; i <= nFace; ++i) {
            G4int n;
            G4int idx[4];
            poly->GetFacet(i, n, idx);

            if (n == 3) {
                indices.push_back(idx[0] - 1);
                indices.push_back(idx[1] - 1);
                indices.push_back(idx[2] - 1);
            }
            else if (n == 4) {
                indices.push_back(idx[0] - 1);
                indices.push_back(idx[1] - 1);
                indices.push_back(idx[2] - 1);

                indices.push_back(idx[0] - 1);
                indices.push_back(idx[2] - 1);
                indices.push_back(idx[3] - 1);
            }
        }

        if (vertices.empty() || indices.empty()) {
            std::cerr << "  -> empty mesh, skipping\n";
            continue;
        }

        // ------------------------------------------------------------
        // write NPZ
        // ------------------------------------------------------------
        std::string filename = outdir + name + ".npz";

        cnpy::npz_save(
            filename,
            "vertices",
            vertices.data(),
            { vertices.size() / 3, 3 },
            "w"
        );

        cnpy::npz_save(
            filename,
            "indices",
            indices.data(),
            { indices.size() / 3, 3 },
            "a"
        );

        std::cout << "  -> wrote " << filename
                  << " (" << vertices.size() / 3
                  << " vertices, " << indices.size() / 3
                  << " triangles)\n";

        ++n_written;
    }

    // ------------------------------------------------------------
    // summary
    // ------------------------------------------------------------
    std::cout << "\nSummary for mode = " << mode << "\n";
    std::cout << "  written meshes     : " << n_written << "\n";
    std::cout << "  skipped (filter)   : " << n_skipped_filter << "\n";
    std::cout << "  skipped (no poly)  : " << n_skipped_poly << "\n";
    std::cout << "Done.\n";

    return 0;
}