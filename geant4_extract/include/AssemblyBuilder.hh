#pragma once

#include "DetectorAssembly.hh"

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4VPhysicalVolume.hh>

#include <cstdint>

class AssemblyBuilder {

public:

    AssemblyBuilder();

    DetectorAssembly Build(
        G4VPhysicalVolume* world
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
};
