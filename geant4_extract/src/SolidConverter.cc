#include "SolidConverter.hh"

#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4DisplacedSolid.hh>
#include <G4GenericPolycone.hh>
#include <G4Orb.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>

#include <TopLoc_Location.hxx>

#include <cmath>

#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <iostream>

// ============================================================
// helper: Geant4 transform -> OCC transform
// ============================================================

static gp_Trsf BuildTransform(
    const G4RotationMatrix& rot,
    const G4ThreeVector& trans
) {

    gp_Trsf trsf;

    double tx = trans.x() / mm;
    double ty = trans.y() / mm;
    double tz = trans.z() / mm;

    trsf.SetValues(

        rot.xx(),
        rot.xy(),
        rot.xz(),
        tx,

        rot.yx(),
        rot.yy(),
        rot.yz(),
        ty,

        rot.zx(),
        rot.zy(),
        rot.zz(),
        tz
    );

    return trsf;
}

// ============================================================
// dispatcher
// ============================================================

TopoDS_Shape SolidConverter::Convert(
    G4VSolid* solid
) {

    // --------------------------------------------------------
    // box
    // --------------------------------------------------------

    if (auto* box =
        dynamic_cast<G4Box*>(solid)) {

        double dx =
            2.0 * box->GetXHalfLength();

        double dy =
            2.0 * box->GetYHalfLength();

        double dz =
            2.0 * box->GetZHalfLength();

        gp_Pnt corner(
            -dx / 2.0,
            -dy / 2.0,
            -dz / 2.0
        );

        return BRepPrimAPI_MakeBox(
            corner,
            dx,
            dy,
            dz
        );
    }

    // --------------------------------------------------------
    // tubs
    // --------------------------------------------------------

    if (auto* tubs =
        dynamic_cast<G4Tubs*>(solid)) {

        double rmin =
            tubs->GetInnerRadius();

        double rmax =
            tubs->GetOuterRadius();

        double h =
            2.0 * tubs->GetZHalfLength();

        double sphi = tubs->GetStartPhiAngle();
        double dphi = tubs->GetDeltaPhiAngle();

        // A G4Tubs may be a phi-segment (deltaphi < 2pi). OCC's
        // BRepPrimAPI_MakeCylinder sweeps the sector from the ax2 X-direction,
        // so orient X to G4's start-phi; the sector then spans
        // [sphi, sphi + dphi].
        gp_Ax2 ax(
            gp_Pnt(0, 0, -h / 2.0),
            gp_Dir(0, 0, 1),
            gp_Dir(std::cos(sphi), std::sin(sphi), 0.0)
        );

        constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
        bool full = (dphi <= 0.0) || (dphi >= kTwoPi - 1e-9);

        auto make_cyl = [&](double r) -> TopoDS_Shape {
            return full
                ? BRepPrimAPI_MakeCylinder(ax, r, h).Shape()
                : BRepPrimAPI_MakeCylinder(ax, r, h, dphi).Shape();
        };

        // outer (solid sector or full cylinder)
        TopoDS_Shape outer = make_cyl(rmax);

        if (rmin <= 0.0)
            return outer;

        // hollow: subtract the inner sector (same phi range) to leave the
        // annular tube segment with its two radial faces.
        TopoDS_Shape inner = make_cyl(rmin);

        return BRepAlgoAPI_Cut(
            outer,
            inner
        );
    }

    // --------------------------------------------------------
    // cons (conical section / frustum shell, optional phi segment)
    // --------------------------------------------------------

    if (auto* cons =
        dynamic_cast<G4Cons*>(solid)) {

        double rmin1 = cons->GetInnerRadiusMinusZ();   // at z = -h/2
        double rmax1 = cons->GetOuterRadiusMinusZ();
        double rmin2 = cons->GetInnerRadiusPlusZ();    // at z = +h/2
        double rmax2 = cons->GetOuterRadiusPlusZ();
        double h     = 2.0 * cons->GetZHalfLength();

        double sphi = cons->GetStartPhiAngle();
        double dphi = cons->GetDeltaPhiAngle();

        constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
        bool full_phi = (dphi <= 0.0) || (dphi >= kTwoPi - 1e-9);

        // base at z = -h/2, axis +z; X-direction oriented to start-phi so the
        // sector spans [sphi, sphi + dphi] (full cone when dphi >= 2pi).
        gp_Ax2 ax(
            gp_Pnt(0, 0, -h / 2.0),
            gp_Dir(0, 0, 1),
            gp_Dir(std::cos(sphi), std::sin(sphi), 0.0)
        );

        // R1 is the radius at the base (z=-h/2), R2 at the top (z=+h/2).
        // OCC's MakeCone is degenerate for equal radii, so fall back to a
        // cylinder there (e.g. a G4Cons with one constant-radius wall).
        auto make_frustum = [&](double rb, double rt) -> TopoDS_Shape {
            if (std::abs(rb - rt) < 1e-9) {
                return full_phi
                    ? BRepPrimAPI_MakeCylinder(ax, rb, h).Shape()
                    : BRepPrimAPI_MakeCylinder(ax, rb, h, dphi).Shape();
            }
            return full_phi
                ? BRepPrimAPI_MakeCone(ax, rb, rt, h).Shape()
                : BRepPrimAPI_MakeCone(ax, rb, rt, h, dphi).Shape();
        };

        TopoDS_Shape outer = make_frustum(rmax1, rmax2);

        // no inner cavity (solid cone/frustum)
        if (rmin1 <= 0.0 && rmin2 <= 0.0)
            return outer;

        // hollow: subtract the inner frustum (same phi range) to leave the
        // conical shell with its annular end caps and phi-plane faces.
        TopoDS_Shape inner = make_frustum(rmin1, rmin2);

        return BRepAlgoAPI_Cut(
            outer,
            inner
        );
    }

    // --------------------------------------------------------
    // orb (full sphere)
    // --------------------------------------------------------

    if (auto* orb =
        dynamic_cast<G4Orb*>(solid)) {

        double r = orb->GetRadius() / mm;

        return BRepPrimAPI_MakeSphere(r).Shape();
    }

    // --------------------------------------------------------
    // sphere (spherical shell / sector)
    // --------------------------------------------------------

    if (auto* sph =
        dynamic_cast<G4Sphere*>(solid)) {

        double rmin = sph->GetInnerRadius() / mm;
        double rmax = sph->GetOuterRadius() / mm;

        double sphi   = sph->GetStartPhiAngle();
        double dphi   = sph->GetDeltaPhiAngle();
        double stheta = sph->GetStartThetaAngle();
        double dtheta = sph->GetDeltaThetaAngle();

        constexpr double kPi    = 3.14159265358979323846;
        constexpr double kTwoPi = 2.0 * kPi;

        bool full_phi   = (dphi   <= 0.0) || (dphi   >= kTwoPi - 1e-9);
        bool full_theta = (dtheta <= 0.0) || (dtheta >= kPi    - 1e-9);
        bool full = full_phi && full_theta;

        // OCC builds a sphere as a revolution: latitude angles a1<a2 in
        // [-pi/2, +pi/2] (a3 = longitude span). G4 theta is the polar angle
        // from +z, so latitude = pi/2 - theta. Orient the ax2 X-direction to
        // G4's start-phi so the longitude sector spans [sphi, sphi + dphi].
        // Ignoring the angle ranges (the old behaviour) turned every spherical
        // sector/segment into a full shell.
        gp_Ax2 ax(
            gp_Pnt(0, 0, 0),
            gp_Dir(0, 0, 1),
            gp_Dir(std::cos(sphi), std::sin(sphi), 0.0)
        );

        double a1 = kPi / 2.0 - (stheta + dtheta);   // latitude at theta_max
        double a2 = kPi / 2.0 - stheta;              // latitude at theta_min

        auto make_sphere = [&](double r) -> TopoDS_Shape {
            if (full)
                return BRepPrimAPI_MakeSphere(ax, r).Shape();
            double a3 = full_phi ? kTwoPi : dphi;
            return BRepPrimAPI_MakeSphere(ax, r, a1, a2, a3).Shape();
        };

        // full sphere/sector with no inner cavity
        if (rmin <= 0.0)
            return make_sphere(rmax);

        // hollow: subtract the inner sector (same angular range) to leave the
        // spherical shell segment with its theta-cone and phi-plane faces.
        TopoDS_Shape outer = make_sphere(rmax);
        TopoDS_Shape inner = make_sphere(rmin);
        return BRepAlgoAPI_Cut(outer, inner);
    }

    // --------------------------------------------------------
    // generic polycone
    // --------------------------------------------------------

    if (auto* poly =
        dynamic_cast<G4GenericPolycone*>(solid)) {

        return ConvertGenericPolycone(poly);
    }

    // --------------------------------------------------------
    // union
    // --------------------------------------------------------

    if (auto* uni =
        dynamic_cast<G4UnionSolid*>(solid)) {

        return ConvertUnion(uni);
    }

    // --------------------------------------------------------
    // subtraction
    // --------------------------------------------------------

    if (auto* sub =
        dynamic_cast<G4SubtractionSolid*>(solid)) {

        return ConvertSubtraction(sub);
    }

    // --------------------------------------------------------
    // displaced solid
    // --------------------------------------------------------

    if (auto* disp =
        dynamic_cast<G4DisplacedSolid*>(solid)) {

        return ConvertDisplaced(disp);
    }

    // --------------------------------------------------------
    // unsupported
    // --------------------------------------------------------

    std::cout
        << "Unsupported solid type: "
        << solid->GetEntityType()
        << std::endl;

    return TopoDS_Shape();
}

// ============================================================
// generic polycone
// ============================================================

TopoDS_Shape SolidConverter::ConvertGenericPolycone(
    G4GenericPolycone* solid
) {

    int n =
        solid->GetNumRZCorner();

    if (n < 2)
        return TopoDS_Shape();

    BRepBuilderAPI_MakeWire wireBuilder;

    for (int i = 0; i < n - 1; ++i) {

        auto p1g =
            solid->GetCorner(i);

        auto p2g =
            solid->GetCorner(i + 1);

        gp_Pnt p1(
            p1g.r / mm,
            0.0,
            p1g.z / mm
        );

        gp_Pnt p2(
            p2g.r / mm,
            0.0,
            p2g.z / mm
        );

        auto edge =
            BRepBuilderAPI_MakeEdge(
                p1,
                p2
            );

        wireBuilder.Add(edge);
    }

    TopoDS_Wire wire =
        wireBuilder.Wire();

    gp_Ax1 axis(
        gp_Pnt(0,0,0),
        gp_Dir(0,0,1)
    );

    return BRepPrimAPI_MakeRevol(
        wire,
        axis,
        2.0 * M_PI
    );
}

// ============================================================
// union
// ============================================================

TopoDS_Shape SolidConverter::ConvertUnion(
    G4UnionSolid* solid
) {

    auto left =
        Convert(
            solid->GetConstituentSolid(0)
        );

    auto right =
        Convert(
            solid->GetConstituentSolid(1)
        );

    if (left.IsNull() || right.IsNull())
        return TopoDS_Shape();

    return BRepAlgoAPI_Fuse(
        left,
        right
    );
}

// ============================================================
// subtraction
// ============================================================

TopoDS_Shape SolidConverter::ConvertSubtraction(
    G4SubtractionSolid* solid
) {

    auto left =
        Convert(
            solid->GetConstituentSolid(0)
        );

    auto right =
        Convert(
            solid->GetConstituentSolid(1)
        );

    if (left.IsNull() || right.IsNull())
        return TopoDS_Shape();

    return BRepAlgoAPI_Cut(
        left,
        right
    );
}

// ============================================================
// displaced solid
// ============================================================

TopoDS_Shape SolidConverter::ConvertDisplaced(
    G4DisplacedSolid* solid
) {

    auto* child =
        solid->GetConstituentMovedSolid();

    auto shape =
        Convert(child);

    if (shape.IsNull())
        return TopoDS_Shape();

    gp_Trsf trsf =
        BuildTransform(

            solid->GetObjectRotation(),

            solid->GetObjectTranslation()
        );

    return shape.Moved(
        TopLoc_Location(trsf)
    );
}
