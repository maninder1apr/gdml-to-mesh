#include "InterfaceExtractor.hh"

#include <OpticalInterface.hh>
#include <VolumeInstance.hh>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>

#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>

#include <BRepClass3d.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>

#include <TopAbs_Orientation.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>

#include <BRepAdaptor_Surface.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>

#include <BRepClass3d_SolidClassifier.hxx>

#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using json = nlohmann::json;

// ============================================================
// HasRealSurface — reject null / empty-COMPOUND / degenerate
// (edge-or-vertex-only) intersection results.
//
// An empty TopoDS_Compound is NOT null, so the old
// `if (result.IsNull())` guard let zero-area touching results
// through. Require at least one face and a total area above a
// small floor.
// ============================================================

static bool HasRealSurface(const TopoDS_Shape &s, double area_floor_mm2) {
  if (s.IsNull())
    return false;

  TopExp_Explorer exp(s, TopAbs_FACE);
  if (!exp.More())
    return false; // empty COMPOUND → 0 faces

  GProp_GProps props;
  BRepGProp::SurfaceProperties(s, props);
  return props.Mass() > area_floor_mm2;
}

// ============================================================
// NormalPointsInto — sample a point + (orientation-corrected)
// normal on the first face of `boundary`, step 1 µm along it,
// and classify against `classifyAgainst`.
//
// Returns true if the outward normal points INTO that solid.
// `ok` is false when no face / no defined normal (caller should
// fall back).
// ============================================================

static bool NormalPointsInto(const TopoDS_Shape &boundary,
                             const TopoDS_Shape &classifyAgainst, bool &ok) {
  ok = false;

  TopExp_Explorer exp(boundary, TopAbs_FACE);
  if (!exp.More())
    return false;

  TopoDS_Face face = TopoDS::Face(exp.Current());
  BRepAdaptor_Surface surf(face);

  // BRepTools::UVBounds gives the trimmed-face bounds, tighter
  // than BRepAdaptor_Surface::First/LastUParameter which returns
  // the bounds of the underlying (untrimmed) surface.
  double u1, u2, v1, v2;
  BRepTools::UVBounds(face, u1, u2, v1, v2);
  double u = 0.5 * (u1 + u2);
  double v = 0.5 * (v1 + v2);

  BRepLProp_SLProps props(surf, u, v, 1, 1e-7);
  if (!props.IsNormalDefined())
    return false;

  gp_Pnt pt;
  surf.D0(u, v, pt);
  gp_Dir n = props.Normal();

  // BRepLProp returns the parametric surface normal and ignores
  // the TopoDS orientation flag — apply it manually.
  if (face.Orientation() == TopAbs_REVERSED)
    n.Reverse();

  constexpr double eps = 1e-4; // 0.1 µm in mm
  gp_Pnt test(pt.X() + eps * n.X(), pt.Y() + eps * n.Y(), pt.Z() + eps * n.Z());

  BRepClass3d_SolidClassifier clf;
  clf.Load(classifyAgainst);
  clf.Perform(test, 1e-7);

  ok = true;
  return clf.State() == TopAbs_IN;
}

// ============================================================
// Orientation results: normals point from inside toward outside.
// ============================================================

struct OrientResult {
  TopoDS_Shape boundary;
  const VolumeInstance *inside;
  const VolumeInstance *outside;
};

// ============================================================
// OrientInterfaceClassifier — sibling↔sibling case.
//
// Labels follow the normal direction (the volume the normal
// points toward is lv_outside). For detector interfaces the
// detector is forced to lv_inside; if it would end up outside
// the mesh winding is reversed. This reproduces the original
// OrientInterface behaviour exactly.
// ============================================================

static OrientResult OrientInterfaceClassifier(const TopoDS_Shape &boundary,
                                              const VolumeInstance &A,
                                              const VolumeInstance &B,
                                              bool A_is_detector,
                                              bool B_is_detector) {
  bool ok = false;
  bool test_in_B = NormalPointsInto(boundary, B.shape, ok);
  if (!ok)
    return {boundary, &A, &B}; // fallback

  const VolumeInstance *vol_in = test_in_B ? &A : &B;
  const VolumeInstance *vol_out = test_in_B ? &B : &A;

  if (!A_is_detector && !B_is_detector)
    return {boundary, vol_in, vol_out};

  bool detector_is_inside =
      (A_is_detector && vol_in == &A) || (B_is_detector && vol_in == &B);

  if (detector_is_inside)
    return {boundary, vol_in, vol_out};

  // Detector ended up as lv_outside → flip mesh, swap labels.
  TopoDS_Shape flipped = boundary;
  flipped.Reverse();
  return {flipped, vol_out, vol_in};
}

// ============================================================
// OrientInterfaceForced — mother↔daughter (containment) case.
//
// The mother is always lv_outside, the daughter lv_inside
// (which is also what optical detectors require). The mesh
// winding is made deterministic: the normal must point from the
// daughter into the mother; if not, the boundary is reversed.
// ============================================================

static OrientResult OrientInterfaceForced(const TopoDS_Shape &boundary,
                                          const VolumeInstance &daughter,
                                          const VolumeInstance &mother) {
  bool ok = false;
  bool points_into_mother = NormalPointsInto(boundary, mother.shape, ok);

  TopoDS_Shape oriented = boundary;
  if (ok && !points_into_mother)
    oriented.Reverse();

  return {oriented, &daughter, &mother};
}

// ============================================================
// EnlargedBox — bbox of a shape inflated by the fuzzy gap so
// that exactly-touching / sub-µm-separated shapes are not culled.
// ============================================================

static Bnd_Box EnlargedBox(const TopoDS_Shape &s, double gap) {
  Bnd_Box b;
  BRepBndLib::Add(s, b);
  b.Enlarge(gap);
  return b;
}

// ============================================================
// VolCache — per-volume bbox + exploded faces + per-face bbox,
// all computed ONCE.
//
// Previously the O(N^2) pair loops rebuilt this data on every
// candidate pair: the sibling cull called EnlargedBox() on both
// shapes (a full BRepBndLib::Add each), and FindSharedFaces
// re-exploded BOTH solids into faces + face boxes on every call.
// With one flat mother (LAr) holding ~469 daughters that is
// ~110k pairs, so each volume's box + face list was rebuilt
// hundreds of times — the actual cause of the multi-hour run.
// Caching turns the broad phase into O(1) precomputed box tests.
// ============================================================

struct VolCache {
  Bnd_Box box; // enlarged by fuzzy_mm
  std::vector<TopoDS_Face> faces;
  std::vector<Bnd_Box> faceBox; // each enlarged by fuzzy_mm
};

// ============================================================
// FindSharedFaces — detect coincident faces between two solids.
//
// Iterates face pairs (bbox pre-filtered) and runs a 2D
// face-on-face BRepAlgoAPI_Common with a fuzzy value. Coincident
// faces yield a real-area patch; these are merged into one
// compound. Returns a null shape if no genuine shared face.
//
// Retained for callers that pass DERIVED shapes not present in
// the per-volume cache (steal_for passes flush footprints and
// existing interface boundaries).
// ============================================================

static TopoDS_Shape FindSharedFaces(const TopoDS_Shape &solidA,
                                    const TopoDS_Shape &solidB, double fuzzy_mm,
                                    double area_floor_mm2) {
  std::vector<TopoDS_Face> facesA, facesB;
  std::vector<Bnd_Box> boxA, boxB;

  for (TopExp_Explorer e(solidA, TopAbs_FACE); e.More(); e.Next()) {
    TopoDS_Face f = TopoDS::Face(e.Current());
    facesA.push_back(f);
    boxA.push_back(EnlargedBox(f, fuzzy_mm));
  }
  for (TopExp_Explorer e(solidB, TopAbs_FACE); e.More(); e.Next()) {
    TopoDS_Face f = TopoDS::Face(e.Current());
    facesB.push_back(f);
    boxB.push_back(EnlargedBox(f, fuzzy_mm));
  }

  BRep_Builder builder;
  TopoDS_Compound out;
  builder.MakeCompound(out);
  bool any = false;

  for (size_t a = 0; a < facesA.size(); ++a) {
    for (size_t b = 0; b < facesB.size(); ++b) {

      if (boxA[a].IsOut(boxB[b]))
        continue; // face-level cull

      BRepAlgoAPI_Common common(facesA[a], facesB[b]);
      common.SetFuzzyValue(fuzzy_mm);
      common.Build();
      if (!common.IsDone())
        continue;

      TopoDS_Shape r = common.Shape();
      if (!HasRealSurface(r, area_floor_mm2))
        continue;

      builder.Add(out, r);
      any = true;
    }
  }

  return any ? TopoDS_Shape(out) : TopoDS_Shape();
}

// ============================================================
// FindSharedFacesCached — identical logic to FindSharedFaces but
// takes precomputed faces + per-face boxes for both solids, so
// nothing is re-exploded or re-boxed per pair. Used on the two
// hot loops where both operands are cached volumes.
// ============================================================

static TopoDS_Shape FindSharedFacesCached(
    const std::vector<TopoDS_Face> &facesA, const std::vector<Bnd_Box> &boxA,
    const std::vector<TopoDS_Face> &facesB, const std::vector<Bnd_Box> &boxB,
    double fuzzy_mm, double area_floor_mm2) {
  BRep_Builder builder;
  TopoDS_Compound out;
  builder.MakeCompound(out);
  bool any = false;

  for (size_t a = 0; a < facesA.size(); ++a) {
    for (size_t b = 0; b < facesB.size(); ++b) {

      if (boxA[a].IsOut(boxB[b]))
        continue; // face-level cull

      BRepAlgoAPI_Common common(facesA[a], facesB[b]);
      common.SetFuzzyValue(fuzzy_mm);
      common.Build();
      if (!common.IsDone())
        continue;

      TopoDS_Shape r = common.Shape();
      if (!HasRealSurface(r, area_floor_mm2))
        continue;

      builder.Add(out, r);
      any = true;
    }
  }

  return any ? TopoDS_Shape(out) : TopoDS_Shape();
}

// ============================================================
// ToFaceCompound — flatten any shape into a TopoDS_Compound of
// its faces.
//
// mother.shape / daughter.shape for G4 primitives (box, tubs,
// cons, sphere, ...) are genuine TopoDS_SOLIDs, not shells. When a
// daughter has no sibling/flush contact to cut, the mother↔daughter
// boundary was previously left as that raw solid — SurfaceMesher's
// "skip volumetric intersections" guard then silently dropped the
// WHOLE interface (ShapeType() == TopAbs_SOLID), even though the
// daughter's own outer surface is exactly the interface we want.
// ============================================================

static TopoDS_Shape ToFaceCompound(const TopoDS_Shape &s) {
  if (s.ShapeType() == TopAbs_COMPOUND)
    return s;

  BRep_Builder builder;
  TopoDS_Compound out;
  builder.MakeCompound(out);
  for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
    builder.Add(out, e.Current());
  return out;
}

// ============================================================
// SubtractPatches — remove shared sibling patches from a
// mother↔daughter boundary so each surface region belongs to
// exactly one interface (no double coverage).
//
// Works face-by-face (each cut is a coplanar 2D operation,
// the robust regime). On a failed/ineffective cut the full face
// is kept and a warning is logged — never silently wrong.
// ============================================================

static TopoDS_Shape SubtractPatches(const TopoDS_Shape &boundary,
                                    const std::vector<TopoDS_Shape> &patches,
                                    double fuzzy_mm, double area_floor_mm2,
                                    int iface_id, bool &warned) {
  if (patches.empty())
    return ToFaceCompound(boundary); // fast path: no sibling contact, but
                                      // still normalize a raw solid down to
                                      // its face compound (see above)

  std::vector<Bnd_Box> pbox(patches.size());
  double sum_patch_area = 0.0;
  for (size_t k = 0; k < patches.size(); ++k) {
    pbox[k] = EnlargedBox(patches[k], fuzzy_mm);
    GProp_GProps pp;
    BRepGProp::SurfaceProperties(patches[k], pp);
    sum_patch_area += pp.Mass();
  }

  BRep_Builder builder;
  TopoDS_Compound out;
  builder.MakeCompound(out);

  double removed_area = 0.0; // how much surface the cuts actually removed

  for (TopExp_Explorer e(boundary, TopAbs_FACE); e.More(); e.Next()) {
    TopoDS_Face f = TopoDS::Face(e.Current());
    Bnd_Box fb = EnlargedBox(f, fuzzy_mm);

    // collect patches whose bbox overlaps this face's bbox. Note a
    // bbox overlap does NOT imply coplanarity (a perpendicular
    // neighbour face shares an edge); the Cut simply removes nothing
    // for non-coplanar patches, which is correct.
    TopoDS_Compound local;
    builder.MakeCompound(local);
    bool hasLocal = false;
    for (size_t k = 0; k < patches.size(); ++k) {
      if (fb.IsOut(pbox[k]))
        continue;
      builder.Add(local, patches[k]);
      hasLocal = true;
    }

    if (!hasLocal) {
      builder.Add(out, f); // face fully borders the mother
      continue;
    }

    GProp_GProps fprops;
    BRepGProp::SurfaceProperties(f, fprops);
    double orig_area = fprops.Mass();

    BRepAlgoAPI_Cut cut(f, local);
    cut.SetFuzzyValue(fuzzy_mm);
    cut.Build();

    if (!cut.IsDone()) {
      // fragile cut failed → keep whole face, warn (possible overlap)
      builder.Add(out, f);
      if (!warned) {
        std::cout << "  [WARN] interface " << iface_id
                  << ": coplanar Cut failed; keeping full "
                  << "mother-daughter face (possible double coverage)\n";
        warned = true;
      }
      continue;
    }

    TopoDS_Shape rem = cut.Shape();
    if (rem.IsNull()) { // face fully shared with sibling → drop
      removed_area += orig_area;
      continue;
    }

    GProp_GProps rprops;
    BRepGProp::SurfaceProperties(rem, rprops);
    double rem_area = rprops.Mass();

    removed_area += (orig_area - rem_area);

    if (rem_area <= area_floor_mm2)
      continue; // fully shared → drop

    builder.Add(out, rem);
  }

  // aggregate sanity check: every shared patch should have been removed
  // from exactly one daughter face. If much less area was removed than the
  // patches cover, some coplanar overlap was missed → possible double
  // coverage. (A per-face check would false-positive on perpendicular
  // neighbour faces that merely share an edge with a patch.)
  //
  // The tolerance is RELATIVE (1 % of the patch area): a coincident curved
  // face and its cut counterpart are tessellated slightly differently, so a
  // full removal lands a fraction of a mm² short of the patch area. An
  // absolute floor would false-alarm on every large curved interface.
  double agg_tol = std::max(area_floor_mm2, 0.01 * sum_patch_area);
  if (removed_area < sum_patch_area - agg_tol && !warned) {
    std::cout << "  [WARN] interface " << iface_id << ": subtracted "
              << removed_area << " mm² but siblings share " << sum_patch_area
              << " mm² (possible double coverage)\n";
    warned = true;
  }

  return out;
}

// ============================================================
// Extract
// ============================================================

void InterfaceExtractor::Extract(
    DetectorAssembly &assembly,
    const std::map<std::string, int> &optical_detectors, double fuzzy_mm) {
  auto &volumes = assembly.volumes;

  constexpr double kAreaFloor = 1e-6; // mm² — reject degenerate patches

  std::cout << "\nChecking interfaces (fuzzy = " << fuzzy_mm << " mm)...\n"
            << std::endl;

  // --------------------------------------------------------
  // index maps. volumes[k].id is NOT guaranteed to equal k if
  // a node was skipped, so resolve ids through an explicit map.
  // children_of is ordered (std::map) for deterministic output.
  // --------------------------------------------------------

  std::unordered_map<uint64_t, size_t> id_to_index;
  for (size_t i = 0; i < volumes.size(); ++i)
    id_to_index[volumes[i].id] = i;

  std::map<uint64_t, std::vector<size_t>> children_of;
  for (size_t i = 0; i < volumes.size(); ++i)
    children_of[volumes[i].mother_id].push_back(i);

  // --------------------------------------------------------
  // per-volume cache: bbox + faces + per-face bbox, computed
  // ONCE (indexed by position in `volumes`, the same index
  // stored in children_of / id_to_index). This replaces the
  // repeated EnlargedBox()/face-explosion the pair loops used
  // to do on every candidate pair.
  // --------------------------------------------------------

  std::vector<VolCache> cache(volumes.size());
  for (size_t i = 0; i < volumes.size(); ++i) {
    BRepBndLib::Add(volumes[i].shape, cache[i].box);
    cache[i].box.Enlarge(fuzzy_mm);
    for (TopExp_Explorer e(volumes[i].shape, TopAbs_FACE); e.More(); e.Next()) {
      TopoDS_Face f = TopoDS::Face(e.Current());
      Bnd_Box b;
      BRepBndLib::Add(f, b);
      b.Enlarge(fuzzy_mm);
      cache[i].faces.push_back(f);
      cache[i].faceBox.push_back(std::move(b));
    }
  }

  // shared sibling patches accumulated per daughter VOLUME id, to be
  // cut out of the corresponding mother↔daughter boundaries
  std::unordered_map<uint64_t, std::vector<TopoDS_Shape>> shared_per_daughter;

  // emitted interfaces grouped by participating VOLUME id. A flush
  // daughter looks up its mother's outward interfaces here to
  // re-attribute ("steal") the surface region it covers.
  std::unordered_map<uint64_t, std::vector<int>> interfaces_of_volume;

  auto is_det = [&](const VolumeInstance &v) {
    return optical_detectors.count(v.name) > 0;
  };

  auto emit = [&](const VolumeInstance &A, const VolumeInstance &B,
                  const OrientResult &orient, bool A_det, bool B_det) {
    OpticalInterface iface;
    iface.id = (int)assembly.interfaces.size();
    iface.volumeA = (int)A.id;
    iface.volumeB = (int)B.id;
    iface.nameA = A.name;
    iface.nameB = B.name;
    iface.materialA = A.material;
    iface.materialB = B.material;
    iface.pv_inside = orient.inside->name;
    iface.pv_outside = orient.outside->name;
    iface.lv_inside = orient.inside->lv_name;
    iface.lv_outside = orient.outside->lv_name;
    iface.mat_inside = orient.inside->material;
    iface.mat_outside = orient.outside->material;
    iface.is_detector = A_det || B_det;
    if (A_det)
      iface.detector_channel = optical_detectors.at(A.name);
    else if (B_det)
      iface.detector_channel = optical_detectors.at(B.name);
    iface.boundary = orient.boundary;
    assembly.interfaces.push_back(iface);

    interfaces_of_volume[A.id].push_back(iface.id);
    interfaces_of_volume[B.id].push_back(iface.id);

    std::cout << "INTERFACE " << iface.id << ":\n  " << iface.pv_inside << " ("
              << iface.mat_inside << ")\n"
              << "  --> " << iface.pv_outside << " (" << iface.mat_outside
              << ")\n"
              << (iface.is_detector
                      ? "  detector channel: " +
                            std::to_string(iface.detector_channel) + "\n"
                      : "")
              << std::endl;
  };

  // --------------------------------------------------------
  // steal_for — re-attribute the part of a mother's OUTWARD
  // interfaces that a flush daughter D actually occupies.
  //
  // For every already-emitted interface M↔X where X lies
  // OUTSIDE M (sibling of M, or M's own mother — not a child
  // of M), the region covered by D's flush footprint is cut
  // out of M↔X and re-emitted as D↔X. Relies on top-down
  // (BFS) traversal so all M↔X already exist, and recurses
  // naturally for flush grand-daughters.
  // --------------------------------------------------------

  auto steal_for = [&](const VolumeInstance &D, const VolumeInstance &M,
                       const TopoDS_Shape &flush) {
    auto it = interfaces_of_volume.find(M.id);
    if (it == interfaces_of_volume.end())
      return;

    std::vector<int> snapshot = it->second; // copy: emit() mutates the map

    for (int idx : snapshot) {
      int vA = assembly.interfaces[idx].volumeA;
      int vB = assembly.interfaces[idx].volumeB;
      int x_id = (vA == (int)M.id) ? vB : vA;

      if (x_id == (int)D.id)
        continue; // skip the M↔D interface itself
      auto xit = id_to_index.find((uint64_t)x_id);
      if (xit == id_to_index.end())
        continue;
      const VolumeInstance &X = volumes[xit->second];

      // only outward-facing interfaces: X must not be a child of M
      if (X.mother_id == M.id)
        continue;

      TopoDS_Shape I_boundary = assembly.interfaces[idx].boundary;
      TopoDS_Shape stolen =
          FindSharedFaces(flush, I_boundary, fuzzy_mm, kAreaFloor);
      if (!HasRealSurface(stolen, kAreaFloor))
        continue;

      // (Fix 1) remove the stolen region from M↔X
      bool warned = false;
      TopoDS_Shape shrunk = SubtractPatches(I_boundary, {stolen}, fuzzy_mm,
                                            kAreaFloor, idx, warned);
      assembly.interfaces[idx].boundary = shrunk; // idx valid; no realloc yet

      if (!HasRealSurface(shrunk, kAreaFloor))
        std::cout << "  [INFO] interface " << idx
                  << " fully re-attributed to daughter " << D.name
                  << " (mother no longer touches " << X.name << " here)\n";

      // (Fix 3) emit D↔X for the stolen patch (sibling-style orientation)
      bool D_det = is_det(D);
      bool X_det = is_det(X);
      OrientResult o = OrientInterfaceClassifier(stolen, D, X, D_det, X_det);
      emit(D, X, o, D_det, X_det); // may realloc assembly.interfaces
    }
  };

  // --------------------------------------------------------
  // process mother groups TOP-DOWN (BFS from the root group,
  // mother_id == kNoMother). Required so a flush daughter can
  // steal from interfaces created one level up.
  // --------------------------------------------------------

  std::vector<uint64_t> ordered_mothers;
  {
    std::vector<uint64_t> q;
    if (children_of.count(VolumeInstance::kNoMother))
      q.push_back(VolumeInstance::kNoMother);

    for (size_t h = 0; h < q.size(); ++h) {
      uint64_t m = q[h];
      ordered_mothers.push_back(m);
      auto cit = children_of.find(m);
      if (cit == children_of.end())
        continue;
      for (size_t di : cit->second) {
        uint64_t cid = volumes[di].id;
        if (children_of.count(cid))
          q.push_back(cid);
      }
    }
    // defensive: append any group not reached from the root (tree should
    // reach all; this only guards against a detached hierarchy)
    for (const auto &g : children_of) {
      bool seen = false;
      for (uint64_t m : ordered_mothers)
        if (m == g.first) {
          seen = true;
          break;
        }
      if (!seen)
        ordered_mothers.push_back(g.first);
    }
  }

  for (uint64_t mother_id : ordered_mothers) {
    const std::vector<size_t> &daughters = children_of.find(mother_id)->second;

    // ====================================================
    // (1) sibling ↔ sibling: shared coincident faces
    // ====================================================

    for (size_t a = 0; a < daughters.size(); ++a) {
      for (size_t b = a + 1; b < daughters.size(); ++b) {
        const VolumeInstance &A = volumes[daughters[a]];
        const VolumeInstance &B = volumes[daughters[b]];

        // Same-material siblings produce no OPTICAL interface, but
        // their contact patch still exists physically and must
        // still be subtracted from each daughter's mother boundary.
        // Skipping the geometry check here (as before) silently
        // starves shared_per_daughter, letting the FULL un-subtracted
        // daughter surface leak into the mother interface instead —
        // this is what inflated DynodeSupportStructure/Dynodes/etc.
        bool same_material = (A.material == B.material);

        // gap-inflated bbox cull so touching pairs survive
        // (both cached boxes are already enlarged by fuzzy_mm)
        if (cache[daughters[a]].box.IsOut(cache[daughters[b]].box))
          continue;

        TopoDS_Shape shared = FindSharedFacesCached(
            cache[daughters[a]].faces, cache[daughters[a]].faceBox,
            cache[daughters[b]].faces, cache[daughters[b]].faceBox, fuzzy_mm,
            kAreaFloor);

        if (!HasRealSurface(shared, kAreaFloor))
          continue;

        // always record the contact patch so it gets subtracted
        // from both daughters' mother boundaries, regardless of
        // whether an optical interface is emitted for it
        shared_per_daughter[A.id].push_back(shared);
        shared_per_daughter[B.id].push_back(shared);

        if (same_material)
          continue; // no optical interface between identical media

        bool A_det = is_det(A);
        bool B_det = is_det(B);

        OrientResult orient =
            OrientInterfaceClassifier(shared, A, B, A_det, B_det);

        emit(A, B, orient, A_det, B_det);
      }
    }

    // ====================================================
    // (2) mother ↔ daughter: containment
    //     (top-level volumes under the excluded World have
    //      mother_id == kNoMother and get no mother interface)
    // ====================================================

    if (mother_id == VolumeInstance::kNoMother)
      continue;

    auto mit = id_to_index.find(mother_id);
    if (mit == id_to_index.end())
      continue;
    const VolumeInstance &mother = volumes[mit->second];
    bool mother_det = is_det(mother);

    for (size_t di : daughters) {
      const VolumeInstance &daughter = volumes[di];

      // Skip identical materials
      if (mother.material == daughter.material)
        continue;

      if (cache[mit->second].box.IsOut(cache[di].box))
        continue;

      // faces of the daughter that lie flush against the mother's
      // OUTER wall (the daughter touches M's boundary from inside).
      // There M has zero thickness: the surface there is really
      // daughter↔(whatever borders M), not mother↔daughter.
      TopoDS_Shape flush = FindSharedFacesCached(
          cache[mit->second].faces, cache[mit->second].faceBox, cache[di].faces,
          cache[di].faceBox, fuzzy_mm, kAreaFloor);
      bool has_flush = HasRealSurface(flush, kAreaFloor);

      // re-attribute M's outward surface under the flush
      // footprint to the daughter. Done BEFORE emitting M↔D so the
      // freshly-created M↔D is not itself a steal target.
      if (has_flush)
        steal_for(daughter, mother, flush);

      // For a daughter fully contained in the mother, the interface
      // is the daughter's outer surface — no boolean needed.
      // This works correctly for tessellated (polyhedron) shapes
      // where BRepAlgoAPI_Common often fails.
      TopoDS_Shape result;

      // Try fast path: use daughter shape directly if bbox is contained
      {
        Bnd_Box motherBox, daughterBox;
        BRepBndLib::Add(mother.shape, motherBox);
        BRepBndLib::Add(daughter.shape, daughterBox);
        motherBox.Enlarge(fuzzy_mm);
        if (!motherBox.IsOut(daughterBox)) {
          // Daughter fully inside mother — daughter surface IS the interface
          result = daughter.shape;
        } else {
          // Partial overlap — try boolean
          BRepAlgoAPI_Common common(mother.shape, daughter.shape);
          common.SetFuzzyValue(fuzzy_mm);
          common.Build();
          if (!common.IsDone())
            continue;
          result = common.Shape();
        }
      }

      if (result.IsNull())
        continue;

      if (!HasRealSurface(result, kAreaFloor))
        continue;

      // patches to cut out of M↔D (no double coverage):
      //   - regions shared with touching siblings
      //   - faces flush with the mother's outer wall, which
      //     carry no mother material and were re-attributed above
      std::vector<TopoDS_Shape> patches;
      auto pit = shared_per_daughter.find(daughter.id);
      if (pit != shared_per_daughter.end())
        patches = pit->second;
      if (has_flush)
        patches.push_back(flush);

      bool warned = false;
      TopoDS_Shape boundary =
          SubtractPatches(result, patches, fuzzy_mm, kAreaFloor,
                          (int)assembly.interfaces.size(), warned);

      if (!HasRealSurface(boundary, kAreaFloor))
        continue; // whole surface shared away with siblings / flush

      if (mother_det) {
        std::cout << "  [WARN] mother volume " << mother.name
                  << " is an optical detector but forced to lv_outside\n";
      }

      OrientResult orient = OrientInterfaceForced(boundary, daughter, mother);

      // daughter passed as A so the detector channel resolves to it
      emit(daughter, mother, orient, is_det(daughter), mother_det);
    }
  }

  // --------------------------------------------------------
  // compaction: a mother↔X interface whose region was FULLY
  // re-attributed to a flush daughter is left with an empty
  // boundary. Drop such interfaces and renumber ids 0..N-1 in
  // emission order (SurfaceMesher/WriteInterfacesJSON key the
  // STL filename off iface.id). When nothing was fully stolen
  // this is a no-op and ids are unchanged.
  // --------------------------------------------------------

  {
    std::vector<OpticalInterface> kept;
    kept.reserve(assembly.interfaces.size());
    for (auto &iface : assembly.interfaces) {
      if (!HasRealSurface(iface.boundary, kAreaFloor)) {
        std::cout << "  [INFO] dropping interface " << iface.id << " ("
                  << iface.pv_inside << " -> " << iface.pv_outside
                  << "): fully re-attributed to a flush daughter\n";
        continue;
      }
      iface.id = (int)kept.size();
      kept.push_back(std::move(iface));
    }
    assembly.interfaces = std::move(kept);
  }

  std::cout << "\nTotal interfaces found: " << assembly.interfaces.size()
            << std::endl;
}

// ============================================================
// WriteInterfacesJSON
// called after SurfaceMesher has populated n_triangles/area_mm2
// ============================================================

void InterfaceExtractor::WriteInterfacesJSON(const DetectorAssembly &assembly,
                                             const std::string &outDir) {
  json arr = json::array();

  for (const auto &iface : assembly.interfaces) {

    json det_id = nullptr;
    if (iface.is_detector)
      det_id = iface.detector_channel;

    json entry;
    entry["id"] = iface.id;
    entry["stl"] =
        "cad/interfaces/interface_" + std::to_string(iface.id) + ".stl";
    entry["pv_inside"] = iface.pv_inside;
    entry["pv_outside"] = iface.pv_outside;
    entry["lv_inside"] = iface.lv_inside;
    entry["lv_outside"] = iface.lv_outside;
    entry["material_inside"] = iface.mat_inside;
    entry["material_outside"] = iface.mat_outside;
    entry["surface"] = iface.is_detector ? "detector" : "default";
    entry["detector_id"] = det_id;
    entry["n_triangles"] = iface.n_triangles;
    entry["area_mm2"] = iface.area_mm2;
    entry["mesh_area_mm2"] = iface.mesh_area_mm2;
    entry["occ_area_mm2"] = iface.occ_area_mm2;
    entry["n_open_edges"] = iface.n_open_edges;
    entry["n_nonmanifold_edges"] = iface.n_nonmanifold_edges;
    entry["n_degenerate_tris"] = iface.n_degenerate_tris;
    entry["mesh_empty"] = iface.mesh_empty;

    arr.push_back(entry);
  }

  std::string path = outDir + "/metadata/interfaces.json";
  std::ofstream f(path);

  if (!f) {
    std::cerr << "ERROR: could not write " << path << std::endl;
    return;
  }

  f << std::setw(2) << arr << std::endl;
  std::cout << "Wrote " << path << std::endl;
}