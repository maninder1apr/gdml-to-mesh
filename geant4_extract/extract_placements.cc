#include <G4GDMLParser.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4VPhysicalVolume.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4SystemOfUnits.hh>

#include "cnpy/cnpy.h"

#include <vector>
#include <string>
#include <iostream>

// ------------------------------------------------------------
// material mapping
// ------------------------------------------------------------
int material_id(const std::string& mat) {
    if (mat.find("LAr") != std::string::npos)  return 0;
    if (mat.find("Germanium") != std::string::npos) return 1;
    if (mat.find("ICPC") != std::string::npos) return 2;
    if (mat.find("PEN") != std::string::npos) return 3;
    return 99;
}

// ------------------------------------------------------------
// mesh mapping
// ------------------------------------------------------------
int mesh_id(const std::string& solid) {
    if (solid == "lar_s")     return 0;
    if (solid == "bege_lv")   return 1;
    if (solid == "icpc_lv")   return 2;
    return -1;
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {

    G4GDMLParser parser;
    parser.Read("../gdml/scarf_pen.gdml");

    auto pvStore = G4PhysicalVolumeStore::GetInstance();
    if (!pvStore) {
        std::cerr << "ERROR: No physical volume store\n";
        return 1;
    }

    std::vector<float> translations;
    std::vector<float> rotations;
    std::vector<int> mesh_ids;
    std::vector<int> material_ids;

    for (auto pv : *pvStore) {
        if (!pv) continue;

        auto lv = pv->GetLogicalVolume();
        if (!lv) continue;

        auto solid = lv->GetSolid();
        if (!solid) continue;

        int mid = mesh_id(solid->GetName());
        if (mid < 0) continue;   // skip world, PEN, etc.

        auto t = pv->GetTranslation();
        auto r = pv->GetRotation();

        translations.insert(translations.end(), {
            (float)(t.x() / mm),
            (float)(t.y() / mm),
            (float)(t.z() / mm)
        });

        if (r) {
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    rotations.push_back((*r)(i,j));
        } else {
            rotations.insert(rotations.end(),
                {1,0,0, 0,1,0, 0,0,1});
        }

        mesh_ids.push_back(mid);

        auto mat = lv->GetMaterial();
        material_ids.push_back(
            mat ? material_id(mat->GetName()) : 99
        );
    }

    cnpy::npz_save("../meshes/placements.npz",
        "translations", translations.data(),
        {translations.size()/3, 3}, "w");

    cnpy::npz_save("../meshes/placements.npz",
        "rotations", rotations.data(),
        {rotations.size()/9, 9}, "a");

    cnpy::npz_save("../meshes/placements.npz",
        "mesh_ids", mesh_ids.data(),
        {mesh_ids.size()}, "a");

    cnpy::npz_save("../meshes/placements.npz",
        "material_ids", material_ids.data(),
        {material_ids.size()}, "a");

    std::cout << "Wrote placements.npz ("
              << mesh_ids.size()
              << " instances)\n";

    return 0;
}
