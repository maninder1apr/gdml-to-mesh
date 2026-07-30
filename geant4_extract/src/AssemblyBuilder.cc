#include "AssemblyBuilder.hh"

#include "SolidConverter.hh"

#include <VolumeInstance.hh>

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4RotationMatrix.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4VPhysicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4VSolid.hh>

#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Transform.hxx>

#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>

#include <gp_Trsf.hxx>

#include <iostream>
#include <set>
#include <string>

// ============================================================
// ctor
// ============================================================

AssemblyBuilder::AssemblyBuilder() {}

// ============================================================
// build full detector assembly
// ============================================================

DetectorAssembly AssemblyBuilder::Build(
    G4VPhysicalVolume* world,
    const std::map<std::string, std::string>& pv_to_lv
) {

    // reset assembly

    assembly_ = DetectorAssembly();

    pv_to_lv_ = pv_to_lv;

    BRep_Builder builder;

    builder.MakeCompound(
        assembly_.assembly
    );

    // --------------------------------------------------------
    // identity transform
    // --------------------------------------------------------

    G4RotationMatrix identity_rot;

    G4ThreeVector identity_trans(
        0,
        0,
        0
    );

    // --------------------------------------------------------
    // recursive traversal
    // --------------------------------------------------------

    Traverse(

        world,

        identity_rot,

        identity_trans,

        VolumeInstance::kNoMother
    );

    return assembly_;
}

// ============================================================
// recursive traversal
// ============================================================

void AssemblyBuilder::Traverse(

    G4VPhysicalVolume* pv,

    const G4RotationMatrix& parent_rot,

    const G4ThreeVector& parent_trans,

    uint64_t parent_id
) {

    if (!pv)
        return;

    auto* lv =
        pv->GetLogicalVolume();

    if (!lv)
        return;

    // id seen by this node's daughters as their mother. Defaults to the
    // inherited parent id so that nodes which produce no VolumeInstance
    // (the World, or any null-shape node) do not orphan their daughters.
    uint64_t this_id = parent_id;

    // --------------------------------------------------------
    // local transform
    // --------------------------------------------------------

    G4RotationMatrix local_rot =
        pv->GetObjectRotationValue();

    G4ThreeVector local_trans =
        pv->GetTranslation();

    // --------------------------------------------------------
    // accumulate hierarchy
    // --------------------------------------------------------

    G4RotationMatrix global_rot =
        parent_rot * local_rot;

    G4ThreeVector global_trans =
        parent_rot * local_trans
        + parent_trans;

    // --------------------------------------------------------
    // skip world
    // --------------------------------------------------------

    bool is_world =
        (pv->GetMotherLogical() == nullptr);

    // Skip volumes whose solid construction is known to be too complex
    // (e.g. 24-deep boolean subtraction chains) or irrelevant for optical sim
    static const std::set<std::string> kSkipNames = {
        "SupportStructure_physical",
        "TubeHolder_logical",
        // fiber<->lar interface extraction reproducibly stalls (burns CPU,
        // zero progress for several minutes+). Root cause found: same
        // per-facet Boolean Cut fragility as mdom's VacuumTube|PMT_0 and
        // legend-200's HPGe<->pen interfaces -- once a fiber accumulates
        // multiple sibling-contact patches (from a genuinely slow but
        // completing sibling-pair phase among ~20 spatially-clustered
        // bend segments), cutting those patches out of its own faceted
        // surface against lar hangs instead of just fragmenting. Fixing
        // SubtractPatches properly would likely resolve all three at
        // once. Skip fibers for now so the rest of legend-200 completes;
        // revisit separately (see memory: project-legend200-fiber-lar-stall).
        "fiber",
    };
    bool skip_this = false;
    for (const auto& s : kSkipNames) {
        if (pv->GetName().find(s) != std::string::npos) {
            skip_this = true;
            break;
        }
    }

    if (!is_world && !skip_this) {

        auto* solid =
            lv->GetSolid();

        if (solid) {

            // ------------------------------------------------
            // convert solid (cached by G4VSolid pointer -- GDML
            // reuses one <volume> definition across many placements
            // via <physvolref>/<volumeref>, e.g. all 24 mDOM PMTs
            // share the identical solid per part. Converting once
            // per unique pointer instead of once per placement skips
            // the repeated tessellation/sewing/RepairAndSolidify
            // cost, which was the actual bottleneck in placement.
            // ------------------------------------------------

            TopoDS_Shape shape;

            auto cached = solid_cache_.find(solid);

            if (cached != solid_cache_.end()) {

                shape = cached->second;
            }
            else {

                SolidConverter converter;

                shape = converter.Convert(solid);

                solid_cache_[solid] = shape;
            }

            if (!shape.IsNull()) {

                // --------------------------------------------
                // translation (mm)
                // --------------------------------------------

                double tx =
                    global_trans.x() / mm;

                double ty =
                    global_trans.y() / mm;

                double tz =
                    global_trans.z() / mm;

                // --------------------------------------------
                // DIRECT affine transform
                //
                // Reproduces:
                //
                // v' = R*v + t
                //
                // exactly like the old working pipeline
                // --------------------------------------------

                gp_Trsf occ_trsf;

                occ_trsf.SetValues(

                    global_rot.xx(),
                    global_rot.xy(),
                    global_rot.xz(),
                    tx,

                    global_rot.yx(),
                    global_rot.yy(),
                    global_rot.yz(),
                    ty,

                    global_rot.zx(),
                    global_rot.zy(),
                    global_rot.zz(),
                    tz
                );

                // --------------------------------------------
                // apply transform
                // --------------------------------------------

                BRepBuilderAPI_Transform transformer(

                    shape,

                    occ_trsf,

                    true
                );

                TopoDS_Shape placed =
                    transformer.Shape();

                // --------------------------------------------
                // add to OCC assembly
                // --------------------------------------------

                BRep_Builder builder;

                builder.Add(

                    assembly_.assembly,

                    placed
                );

                // --------------------------------------------
                // semantic volume info
                // --------------------------------------------

                VolumeInstance inst;

                inst.id =
                    next_volume_id_++;

                inst.mother_id =
                    parent_id;

                // this node's daughters border THIS node as their mother
                this_id = inst.id;

                inst.name =
                    pv->GetName();

                // Use the GDML <volume name=...> via the store,
                // since the GDML parser renames LVs to match physvol names.
                {
                    auto* store = G4LogicalVolumeStore::GetInstance();
                    inst.lv_name = lv->GetName(); // fallback
                    for (auto* slv : *store) {
                        if (slv == lv) {
                            inst.lv_name = slv->GetName();
                            break;
                        }
                    }
                }

                if (lv->GetMaterial()) {

                    inst.material =
                        lv->GetMaterial()->GetName();
                }
                else {

                    inst.material =
                        "UNKNOWN";
                }

                inst.shape =
                    placed;

                assembly_.volumes.push_back(
                    inst
                );

                // --------------------------------------------
                // debug
                // --------------------------------------------

                std::cout
                    << "Placed: "
                    << inst.name
                    << "  material="
                    << inst.material
                    << "  at  "
                    << tx << " "
                    << ty << " "
                    << tz
                    << std::endl;

                std::cout
                    << "Rotation:\n"
                    << global_rot.xx() << " "
                    << global_rot.xy() << " "
                    << global_rot.xz() << "\n"

                    << global_rot.yx() << " "
                    << global_rot.yy() << " "
                    << global_rot.yz() << "\n"

                    << global_rot.zx() << " "
                    << global_rot.zy() << " "
                    << global_rot.zz()
                    << "\n";
            }
            else {

                std::cout
                    << "Skipped null shape: "
                    << pv->GetName()
                    << std::endl;
            }
        }
    }

    // --------------------------------------------------------
    // recurse daughters
    // --------------------------------------------------------

    for (
        int i = 0;
        i < lv->GetNoDaughters();
        ++i
    ) {

        Traverse(

            lv->GetDaughter(i),

            global_rot,

            global_trans,

            this_id
        );
    }
}
