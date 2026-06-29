#pragma once

#include <G4VSolid.hh>

#include <G4GenericPolycone.hh>
#include <G4IntersectionSolid.hh>
#include <G4Orb.hh>
#include <G4Polycone.hh>
#include <G4Ellipsoid.hh>
#include <G4Torus.hh>
#include <G4Sphere.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4DisplacedSolid.hh>
#include <G4TessellatedSolid.hh>

#include <TopoDS_Shape.hxx>

class SolidConverter {

public:

    TopoDS_Shape Convert(G4VSolid* solid);

    // Universal fallback: uses G4VSolid::GetPolyhedron() to tessellate
    // any solid Geant4 can render — ellipsoids, polycones, boolean chains,
    // tessellated solids. Returns null only if G4 itself cannot tessellate.
    TopoDS_Shape ConvertViaPolyhedron(G4VSolid* solid);

    // Native OCC converters (fast, exact)
    TopoDS_Shape ConvertGenericPolycone(G4GenericPolycone* solid);
    TopoDS_Shape ConvertPolycone(G4Polycone* solid);
    TopoDS_Shape ConvertEllipsoid(G4Ellipsoid* solid);
    TopoDS_Shape ConvertTorus(G4Torus* solid);
    TopoDS_Shape ConvertUnion(G4UnionSolid* solid);
    TopoDS_Shape ConvertIntersection(G4IntersectionSolid* solid);
    TopoDS_Shape ConvertDisplaced(G4DisplacedSolid* solid);
    TopoDS_Shape ConvertSubtraction(G4SubtractionSolid* solid);
};
