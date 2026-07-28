#include "SolidConverter.hh"

#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4DisplacedSolid.hh>
#include <G4Ellipsoid.hh>
#include <G4GenericPolycone.hh>
#include <G4IntersectionSolid.hh>
#include <G4Orb.hh>
#include <G4Polycone.hh>
#include <G4Polyhedron.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4TessellatedSolid.hh>
#include <G4ThreeVector.hh>
#include <G4Torus.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4VFacet.hh>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>

#include <BRep_Builder.hxx>

#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Solid.hxx>

#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
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
#include <vector>

// ============================================================
// helper: Geant4 transform -> OCC transform
// ============================================================

static gp_Trsf BuildTransform(const G4RotationMatrix &rot,
                              const G4ThreeVector &trans) {
  gp_Trsf trsf;
  double tx = trans.x() / mm;
  double ty = trans.y() / mm;
  double tz = trans.z() / mm;
  trsf.SetValues(rot.xx(), rot.xy(), rot.xz(), tx, rot.yx(), rot.yy(), rot.yz(),
                 ty, rot.zx(), rot.zy(), rot.zz(), tz);
  return trsf;
}

// ============================================================
// RepairAndSolidify — orient-fix a shell, then promote it to a
// genuine TopoDS_Solid.
//
// A source tessellated solid can have inconsistently-wound facets
// (G4's own "negative cubic volume, please check orientation of
// facets!" warning) — some facets face outward, some inward. Sewing
// alone does not fix this: BRepAlgoAPI_Fuse/Cut/Common run against
// an inconsistently-oriented shell/solid lose track of "inside vs
// outside" and can produce chaotic self-intersecting garbage (this
// is what turned DynodeSupportStructure into a triangle soup).
// ShapeFix_Shell re-orients a shell's faces to be mutually
// consistent; ShapeFix_Solid::SolidFromShell then builds a properly
// oriented solid from that (and can bridge small gaps). Falls back
// to the orientation-fixed shell, or the original shape, if either
// step doesn't produce a solid — never worse than before.
// ============================================================

static TopoDS_Shape RepairAndSolidify(const TopoDS_Shape &s) {
  if (s.ShapeType() != TopAbs_SHELL)
    return s;

  TopoDS_Shell shell = TopoDS::Shell(s);

  Handle(ShapeFix_Shell) shellFix = new ShapeFix_Shell();
  shellFix->Init(shell);
  shellFix->Perform();
  TopoDS_Shape fixed = shellFix->Shell();

  if (fixed.ShapeType() != TopAbs_SHELL)
    return fixed;

  Handle(ShapeFix_Solid) solidFix = new ShapeFix_Solid();
  TopoDS_Shape solid = solidFix->SolidFromShell(TopoDS::Shell(fixed));
  if (!solid.IsNull() && solid.ShapeType() == TopAbs_SOLID)
    return solid;

  return fixed;
}

// ============================================================
// FinishSewn — promote a sewn shell to a genuine TopoDS_Solid.
//
// BRepBuilderAPI_Sewing::SewedShape() only ever returns a bare
// SHELL (or a COMPOUND of shells, if the facets didn't all join
// into one connected piece) — never a SOLID. Downstream union/
// subtraction/intersection run BRepAlgoAPI_Fuse/Cut/Common against
// this result; feeding those a shell instead of a proper solid is
// a much more fragile regime and a plausible source of the small
// unsewn/mismatched edges that showed up as open/non-manifold
// edges in meshed interfaces (e.g. G4Ellipsoid-derived solids
// unioned with tubes). Wrapping each closed shell in a properly
// oriented solid (see RepairAndSolidify) gives later booleans a
// real, consistently-oriented solid to work with; a shell that
// can't be repaired falls back to itself unchanged (never worse
// than before).
// ============================================================

static TopoDS_Shape FinishSewn(BRepBuilderAPI_Sewing &sew) {
  sew.Perform();
  TopoDS_Shape sewed = sew.SewedShape();

  auto to_solid = [](const TopoDS_Shape &s) -> TopoDS_Shape {
    return RepairAndSolidify(s);
  };

  if (sewed.ShapeType() == TopAbs_SHELL)
    return to_solid(sewed);

  if (sewed.ShapeType() == TopAbs_COMPOUND) {
    BRep_Builder builder;
    TopoDS_Compound out;
    builder.MakeCompound(out);
    bool any = false;
    for (TopExp_Explorer e(sewed, TopAbs_SHELL); e.More(); e.Next()) {
      builder.Add(out, to_solid(e.Current()));
      any = true;
    }
    if (any)
      return out;
  }

  return sewed;
}

// ============================================================
// dispatcher
// ============================================================

TopoDS_Shape SolidConverter::Convert(G4VSolid *solid) {
  // ── box ─────────────────────────────────────────────────
  if (auto *box = dynamic_cast<G4Box *>(solid)) {
    double dx = 2.0 * box->GetXHalfLength();
    double dy = 2.0 * box->GetYHalfLength();
    double dz = 2.0 * box->GetZHalfLength();
    gp_Pnt corner(-dx / 2.0, -dy / 2.0, -dz / 2.0);
    return BRepPrimAPI_MakeBox(corner, dx, dy, dz);
  }

  // ── tubs ─────────────────────────────────────────────────
  if (auto *tubs = dynamic_cast<G4Tubs *>(solid)) {
    double rmin = tubs->GetInnerRadius();
    double rmax = tubs->GetOuterRadius();
    double h = 2.0 * tubs->GetZHalfLength();
    double sphi = tubs->GetStartPhiAngle();
    double dphi = tubs->GetDeltaPhiAngle();
    gp_Ax2 ax(gp_Pnt(0, 0, -h / 2.0), gp_Dir(0, 0, 1),
              gp_Dir(std::cos(sphi), std::sin(sphi), 0.0));
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    bool full = (dphi <= 0.0) || (dphi >= kTwoPi - 1e-9);
    auto make_cyl = [&](double r) -> TopoDS_Shape {
      return full ? BRepPrimAPI_MakeCylinder(ax, r, h).Shape()
                  : BRepPrimAPI_MakeCylinder(ax, r, h, dphi).Shape();
    };
    TopoDS_Shape outer = make_cyl(rmax);
    if (rmin <= 0.0)
      return outer;
    return BRepAlgoAPI_Cut(outer, make_cyl(rmin));
  }

  // ── cons ─────────────────────────────────────────────────
  if (auto *cons = dynamic_cast<G4Cons *>(solid)) {
    double rmin1 = cons->GetInnerRadiusMinusZ();
    double rmax1 = cons->GetOuterRadiusMinusZ();
    double rmin2 = cons->GetInnerRadiusPlusZ();
    double rmax2 = cons->GetOuterRadiusPlusZ();
    double h = 2.0 * cons->GetZHalfLength();
    double sphi = cons->GetStartPhiAngle();
    double dphi = cons->GetDeltaPhiAngle();
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    bool full_phi = (dphi <= 0.0) || (dphi >= kTwoPi - 1e-9);
    gp_Ax2 ax(gp_Pnt(0, 0, -h / 2.0), gp_Dir(0, 0, 1),
              gp_Dir(std::cos(sphi), std::sin(sphi), 0.0));
    auto make_frustum = [&](double rb, double rt) -> TopoDS_Shape {
      if (std::abs(rb - rt) < 1e-9)
        return full_phi ? BRepPrimAPI_MakeCylinder(ax, rb, h).Shape()
                        : BRepPrimAPI_MakeCylinder(ax, rb, h, dphi).Shape();
      return full_phi ? BRepPrimAPI_MakeCone(ax, rb, rt, h).Shape()
                      : BRepPrimAPI_MakeCone(ax, rb, rt, h, dphi).Shape();
    };
    TopoDS_Shape outer = make_frustum(rmax1, rmax2);
    if (rmin1 <= 0.0 && rmin2 <= 0.0)
      return outer;
    return BRepAlgoAPI_Cut(outer, make_frustum(rmin1, rmin2));
  }

  // ── orb ──────────────────────────────────────────────────
  if (auto *orb = dynamic_cast<G4Orb *>(solid))
    return BRepPrimAPI_MakeSphere(orb->GetRadius() / mm).Shape();

  // ── sphere ───────────────────────────────────────────────
  if (auto *sph = dynamic_cast<G4Sphere *>(solid)) {
    double rmin = sph->GetInnerRadius() / mm;
    double rmax = sph->GetOuterRadius() / mm;
    double sphi = sph->GetStartPhiAngle();
    double dphi = sph->GetDeltaPhiAngle();
    double stheta = sph->GetStartThetaAngle();
    double dtheta = sph->GetDeltaThetaAngle();
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kTwoPi = 2.0 * kPi;
    bool full_phi = (dphi <= 0.0) || (dphi >= kTwoPi - 1e-9);
    bool full_theta = (dtheta <= 0.0) || (dtheta >= kPi - 1e-9);
    bool full = full_phi && full_theta;
    gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1),
              gp_Dir(std::cos(sphi), std::sin(sphi), 0.0));
    double a1 = kPi / 2.0 - (stheta + dtheta);
    double a2 = kPi / 2.0 - stheta;
    auto make_sphere = [&](double r) -> TopoDS_Shape {
      if (full)
        return BRepPrimAPI_MakeSphere(ax, r).Shape();
      double a3 = full_phi ? kTwoPi : dphi;
      return BRepPrimAPI_MakeSphere(ax, r, a1, a2, a3).Shape();
    };
    if (rmin <= 0.0)
      return make_sphere(rmax);
    return BRepAlgoAPI_Cut(make_sphere(rmax), make_sphere(rmin));
  }

  // ── polycone ─────────────────────────────────────────────
  if (auto *poly = dynamic_cast<G4Polycone *>(solid))
    return ConvertPolycone(poly);

  // ── generic polycone ─────────────────────────────────────
  if (auto *poly = dynamic_cast<G4GenericPolycone *>(solid))
    return ConvertGenericPolycone(poly);

  // ── ellipsoid → G4 polyhedron fallback ───────────────────
  if (auto *ell = dynamic_cast<G4Ellipsoid *>(solid))
    return ConvertViaPolyhedron(ell);

  // ── torus ────────────────────────────────────────────────
  if (auto *tor = dynamic_cast<G4Torus *>(solid))
    return ConvertTorus(tor);

  // ── union ────────────────────────────────────────────────
  if (auto *uni = dynamic_cast<G4UnionSolid *>(solid))
    return ConvertUnion(uni);

  // ── subtraction ──────────────────────────────────────────
  if (auto *sub = dynamic_cast<G4SubtractionSolid *>(solid))
    return ConvertSubtraction(sub);

  // ── intersection ─────────────────────────────────────────
  if (auto *isec = dynamic_cast<G4IntersectionSolid *>(solid))
    return ConvertIntersection(isec);

  // ── displaced solid ──────────────────────────────────────
  if (auto *disp = dynamic_cast<G4DisplacedSolid *>(solid))
    return ConvertDisplaced(disp);

  // ── universal fallback: G4 polyhedron tessellation ───────
  // Handles G4TessellatedSolid, any future solid types, and
  // anything not explicitly listed above.
  TopoDS_Shape fallback = ConvertViaPolyhedron(solid);
  if (!fallback.IsNull())
    return fallback;

  std::cout << "Unsupported solid type: " << solid->GetEntityType()
            << std::endl;
  return TopoDS_Shape();
}

// ============================================================
// ConvertViaPolyhedron
//
// Universal fallback: uses G4VSolid::GetPolyhedron() to
// tessellate any solid Geant4 can render — G4Ellipsoid,
// G4Polycone in boolean chains, G4TessellatedSolid, etc.
// Builds an OCC sewn shell from the triangular facets.
// ============================================================

TopoDS_Shape SolidConverter::ConvertViaPolyhedron(G4VSolid *solid) {
  if (!solid)
    return TopoDS_Shape();

  // ── G4TessellatedSolid: exact CAD facets, no approximation ──
  if (auto *tess = dynamic_cast<G4TessellatedSolid *>(solid)) {
    BRepBuilderAPI_Sewing sew;
    sew.SetTolerance(1e-3);
    bool any = false;

    for (int i = 0; i < tess->GetNumberOfFacets(); ++i) {
      const G4VFacet *f = tess->GetFacet(i);
      if (!f || f->GetNumberOfVertices() < 3)
        continue;

      auto v = [&](int k) -> gp_Pnt {
        G4ThreeVector p = f->GetVertex(k);
        return gp_Pnt(p.x() / mm, p.y() / mm, p.z() / mm);
      };

      // Fan triangulation from vertex 0
      for (int k = 1; k + 1 < f->GetNumberOfVertices(); ++k) {
        gp_Pnt p0 = v(0), p1 = v(k), p2 = v(k + 1);
        if (p0.Distance(p1) < 1e-9 || p1.Distance(p2) < 1e-9 ||
            p0.Distance(p2) < 1e-9)
          continue;
        try {
          BRepBuilderAPI_MakePolygon tri;
          tri.Add(p0);
          tri.Add(p1);
          tri.Add(p2);
          tri.Close();
          if (!tri.IsDone())
            continue;
          BRepBuilderAPI_MakeFace face(tri.Wire(), true);
          if (!face.IsDone())
            continue;
          sew.Add(face.Face());
          any = true;
        } catch (...) {
          continue;
        }
      }
    }

    if (!any)
      return TopoDS_Shape();
    return FinishSewn(sew);
  }

  // ── General case: G4Polyhedron tessellation ──────────────
  G4Polyhedron::SetNumberOfRotationSteps(48);
  G4Polyhedron *poly = solid->GetPolyhedron();
  G4Polyhedron::SetNumberOfRotationSteps(24);
  if (!poly) {
    std::cout << "ConvertViaPolyhedron: no G4Polyhedron for "
              << solid->GetEntityType() << std::endl;
    return TopoDS_Shape();
  }

  BRepBuilderAPI_Sewing sew;
  sew.SetTolerance(1e-3);
  bool any = false;

  // GetNextFacet fills nodes[0..n-1] with face vertex positions.
  // Returns false for the last face — so we loop until it returns false,
  // processing that last face too.
  G4int nVerts;
  G4Point3D nodes[4];
  G4bool notLast = true;
  while (true) {
    notLast = poly->GetNextFacet(nVerts, nodes);

    // Fan triangulation from vertex 0
    for (int k = 1; k + 1 < nVerts; ++k) {
      gp_Pnt p0(nodes[0].x() / mm, nodes[0].y() / mm, nodes[0].z() / mm);
      gp_Pnt p1(nodes[k].x() / mm, nodes[k].y() / mm, nodes[k].z() / mm);
      gp_Pnt p2(nodes[k + 1].x() / mm, nodes[k + 1].y() / mm,
                nodes[k + 1].z() / mm);
      if (p0.Distance(p1) < 1e-9 || p1.Distance(p2) < 1e-9 ||
          p0.Distance(p2) < 1e-9)
        continue;
      try {
        BRepBuilderAPI_MakePolygon tri;
        tri.Add(p0);
        tri.Add(p1);
        tri.Add(p2);
        tri.Close();
        if (!tri.IsDone())
          continue;
        BRepBuilderAPI_MakeFace face(tri.Wire(), true);
        if (!face.IsDone())
          continue;
        sew.Add(face.Face());
        any = true;
      } catch (...) {
        continue;
      }
    }

    if (!notLast)
      break;
  }

  if (!any)
    return TopoDS_Shape();
  return FinishSewn(sew);
}

// ============================================================
// G4Polycone — build closed rz profile and revolve
// ============================================================

TopoDS_Shape SolidConverter::ConvertPolycone(G4Polycone *solid) {
  const G4PolyconeHistorical *params = solid->GetOriginalParameters();
  int n = params->Num_z_planes;
  if (n < 2)
    return TopoDS_Shape();

  std::vector<double> zv(n), rminv(n), rmaxv(n);
  for (int i = 0; i < n; ++i) {
    zv[i] = params->Z_values[i] / mm;
    rminv[i] = params->Rmin[i] / mm;
    rmaxv[i] = params->Rmax[i] / mm;
  }

  BRepBuilderAPI_MakeWire wb;
  auto addSeg = [&](double r1, double z1, double r2, double z2) {
    if (std::abs(r1 - r2) < 1e-12 && std::abs(z1 - z2) < 1e-12)
      return;
    wb.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(r1, 0, z1), gp_Pnt(r2, 0, z2)));
  };

  for (int i = 0; i < n - 1; ++i)
    addSeg(rmaxv[i], zv[i], rmaxv[i + 1], zv[i + 1]);
  addSeg(rmaxv[n - 1], zv[n - 1], rminv[n - 1], zv[n - 1]);
  for (int i = n - 1; i > 0; --i)
    addSeg(rminv[i], zv[i], rminv[i - 1], zv[i - 1]);
  addSeg(rminv[0], zv[0], rmaxv[0], zv[0]);

  if (!wb.IsDone())
    return TopoDS_Shape();

  double sphi = solid->GetStartPhi();
  double dphi = solid->GetEndPhi() - sphi;
  constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
  bool full = (dphi >= kTwoPi - 1e-9);

  gp_Ax1 axis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
  return BRepPrimAPI_MakeRevol(wb.Wire(), axis, full ? kTwoPi : dphi).Shape();
}

// ============================================================
// generic polycone
// ============================================================

TopoDS_Shape SolidConverter::ConvertGenericPolycone(G4GenericPolycone *solid) {
  int n = solid->GetNumRZCorner();
  if (n < 2)
    return TopoDS_Shape();

  BRepBuilderAPI_MakeWire wb;
  for (int i = 0; i < n - 1; ++i) {
    auto p1g = solid->GetCorner(i);
    auto p2g = solid->GetCorner(i + 1);
    wb.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(p1g.r / mm, 0.0, p1g.z / mm),
                                   gp_Pnt(p2g.r / mm, 0.0, p2g.z / mm)));
  }

  gp_Ax1 axis(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
  return BRepPrimAPI_MakeRevol(wb.Wire(), axis, 2.0 * M_PI);
}

// ============================================================
// G4Ellipsoid — routed through G4 polyhedron tessellation
// ============================================================

TopoDS_Shape SolidConverter::ConvertEllipsoid(G4Ellipsoid *solid) {
  return ConvertViaPolyhedron(solid);
}

// ============================================================
// G4Torus
// ============================================================

TopoDS_Shape SolidConverter::ConvertTorus(G4Torus *solid) {
  double rmin = solid->GetRmin() / mm;
  double rmax = solid->GetRmax() / mm;
  double rtor = solid->GetRtor() / mm;
  double sphi = solid->GetSPhi();
  double dphi = solid->GetDPhi();
  constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
  bool full = (dphi >= kTwoPi - 1e-9);
  gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1),
            gp_Dir(std::cos(sphi), std::sin(sphi), 0.0));
  auto make_torus = [&](double r) -> TopoDS_Shape {
    return full ? BRepPrimAPI_MakeTorus(ax, rtor, r).Shape()
                : BRepPrimAPI_MakeTorus(ax, rtor, r, dphi).Shape();
  };
  TopoDS_Shape outer = make_torus(rmax);
  if (rmin <= 0.0)
    return outer;
  return BRepAlgoAPI_Cut(outer, make_torus(rmin));
}

// ============================================================
// union
// ============================================================

TopoDS_Shape SolidConverter::ConvertUnion(G4UnionSolid *solid) {
  auto left = Convert(solid->GetConstituentSolid(0));
  auto right = Convert(solid->GetConstituentSolid(1));

  // If either part failed, return what we have rather than null
  if (left.IsNull() && right.IsNull())
    return TopoDS_Shape();
  if (left.IsNull())
    return right;
  if (right.IsNull())
    return left;

  BRepAlgoAPI_Fuse fuse(left, right);
  fuse.SetFuzzyValue(1e-4);
  fuse.Build();
  if (!fuse.IsDone() || fuse.Shape().IsNull()) {
    // Return compound of both parts — not a true boolean
    // but preserves all geometry for meshing purposes
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    builder.Add(compound, left);
    builder.Add(compound, right);
    return compound;
  }
  return fuse.Shape();
}

// ============================================================
// subtraction
// ============================================================

TopoDS_Shape SolidConverter::ConvertSubtraction(G4SubtractionSolid *solid) {
  auto left = Convert(solid->GetConstituentSolid(0));
  auto right = Convert(solid->GetConstituentSolid(1));

  // If subtrahend failed, return base shape uncut
  if (left.IsNull())
    return TopoDS_Shape();
  if (right.IsNull())
    return left;

  BRepAlgoAPI_Cut cut(left, right);
  cut.SetFuzzyValue(1e-4);
  cut.Build();
  if (!cut.IsDone() || cut.Shape().IsNull())
    return left; // return uncut base rather than null
  return cut.Shape();
}

// ============================================================
// intersection
// ============================================================

TopoDS_Shape SolidConverter::ConvertIntersection(G4IntersectionSolid *solid) {
  auto left = Convert(solid->GetConstituentSolid(0));
  auto right = Convert(solid->GetConstituentSolid(1));
  if (left.IsNull() || right.IsNull())
    return TopoDS_Shape();
  BRepAlgoAPI_Common common(left, right);
  common.SetFuzzyValue(1e-4);
  common.Build();
  if (!common.IsDone() || common.Shape().IsNull())
    return TopoDS_Shape();
  return common.Shape();
}

// ============================================================
// displaced solid
// ============================================================

TopoDS_Shape SolidConverter::ConvertDisplaced(G4DisplacedSolid *solid) {
  auto *child = solid->GetConstituentMovedSolid();
  auto shape = Convert(child);
  if (shape.IsNull())
    return TopoDS_Shape();
  gp_Trsf trsf =
      BuildTransform(solid->GetObjectRotation(), solid->GetObjectTranslation());
  return shape.Moved(TopLoc_Location(trsf));
}
