#pragma once

#include "DetectorAssembly.hh"

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>

#include <TopoDS_Shape.hxx>

#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>

class AssemblyBuilder {

public:

    AssemblyBuilder();

    DetectorAssembly Build(
        G4VPhysicalVolume* world,
        const std::map<std::string, std::string>& pv_to_lv = {}
    );

private:

    void Traverse(

        G4VPhysicalVolume* pv,

        const G4RotationMatrix& parent_rot,

        const G4ThreeVector& parent_trans,

        uint64_t parent_id
    );

private:

    DetectorAssembly assembly_;

    uint64_t next_volume_id_ = 0;

    std::map<std::string, std::string> pv_to_lv_;

    // GDML reuses one <volume> definition across many placements
    // (<physvolref>/<volumeref>) -- e.g. all 24 mDOM PMTs share the
    // identical G4VSolid pointer per part. SolidConverter::Convert()
    // (tessellation + sewing + RepairAndSolidify) is expensive and was
    // being re-run once per PLACEMENT even though the input is
    // byte-identical; cache by solid pointer so each unique local
    // shape is converted once and every repeat placement just reuses
    // (then transforms) the cached copy.
    std::unordered_map<const G4VSolid*, TopoDS_Shape> solid_cache_;
};
