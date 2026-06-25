#pragma once

#include <G4VSolid.hh>

#include <G4GenericPolycone.hh>
#include <G4Orb.hh>
#include <G4Sphere.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4DisplacedSolid.hh>

#include <TopoDS_Shape.hxx>

class SolidConverter {

public:

    TopoDS_Shape Convert(
        G4VSolid* solid
    );

    TopoDS_Shape ConvertGenericPolycone(
        G4GenericPolycone* solid
    );

    TopoDS_Shape ConvertUnion(
        G4UnionSolid* solid
    );

    TopoDS_Shape ConvertDisplaced(
        G4DisplacedSolid* solid
    );

    TopoDS_Shape ConvertSubtraction(
        G4SubtractionSolid* solid
    );
};
