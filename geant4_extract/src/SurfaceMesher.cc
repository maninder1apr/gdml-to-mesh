#include "SurfaceMesher.hh"

#include "DetectorAssembly.hh"
#include "InterfaceExtractor.hh"
#include "OpticalInterface.hh"

#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <StlAPI_Writer.hxx>

#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>
#include <cmath>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>

#include <TopExp_Explorer.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <TopLoc_Location.hxx>

#include <TopAbs_ShapeEnum.hxx>

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

// Set to true to also dump a .brep per interface (debug only). These
// are not needed for plotting or for Theia and roughly double the file
// count, so the default run skips them.
static constexpr bool kWriteBrep = false;

// ============================================================
// mesh-quality audit helpers
//
// A boundary compound's faces each carry their own independent
// Poly_Triangulation (own local node indices), so a shared edge
// between two faces can't be found by comparing node indices —
// it has to be identified by position. Edges are keyed by their
// two endpoints quantized to kEdgeQuantum and canonically ordered
// so a given physical edge hashes the same regardless of which
// triangle/winding visited it first. Count == 1 → open (boundary)
// edge; count > 2 → non-manifold (overlapping/duplicated) edge;
// count == 2 is the expected interior case for a closed surface.
// ============================================================

namespace {

constexpr double kEdgeQuantum = 1e-4; // mm — matches the 0.1 µm eps used
                                       // for normal-direction sampling
constexpr double kDegenerateAreaFloor = 1e-9; // mm^2

struct EdgeKey {
  long long ax, ay, az, bx, by, bz;
  bool operator==(const EdgeKey &o) const {
    return ax == o.ax && ay == o.ay && az == o.az && bx == o.bx &&
           by == o.by && bz == o.bz;
  }
};

struct EdgeKeyHash {
  size_t operator()(const EdgeKey &k) const {
    size_t h = 1469598103934665603ULL;
    auto mix = [&](long long v) {
      h ^= static_cast<size_t>(v);
      h *= 1099511628211ULL;
    };
    mix(k.ax);
    mix(k.ay);
    mix(k.az);
    mix(k.bx);
    mix(k.by);
    mix(k.bz);
    return h;
  }
};

long long QuantizeCoord(double v) {
  return static_cast<long long>(std::llround(v / kEdgeQuantum));
}

EdgeKey MakeEdgeKey(const gp_Pnt &p, const gp_Pnt &q) {
  long long ax = QuantizeCoord(p.X()), ay = QuantizeCoord(p.Y()),
            az = QuantizeCoord(p.Z());
  long long bx = QuantizeCoord(q.X()), by = QuantizeCoord(q.Y()),
            bz = QuantizeCoord(q.Z());
  if (std::tie(ax, ay, az) > std::tie(bx, by, bz)) {
    std::swap(ax, bx);
    std::swap(ay, by);
    std::swap(az, bz);
  }
  return {ax, ay, az, bx, by, bz};
}

// ============================================================
// fiber-layer interface bundling
//
// A fiber system can produce thousands of individually-tiny interfaces
// (every fiber's own core<->cl1, cl1<->cl2, ... contact), each of which
// would otherwise get its own STL file -- LEGEND-1000 alone has ~12000
// fiber placements. Every interface between the same pair of fiber-
// layer materials shares ONE combined STL instead (accumulated across
// threads into a per-pair compound, written once after the main loop);
// each interface's own area/mesh-quality/pv names in interfaces.json
// are completely unaffected, only iface.stl_override changes.
// ============================================================

bool IsFiberLayerMaterial(const std::string &mat) {
  return mat == "tpb_on_fibers" || mat == "pmma_cl2" || mat == "pmma" ||
         mat == "ps_fibers";
}

// single shared filename for every fiber-layer interface, regardless of
// which specific material pair it is (core<->cl1, cl1<->cl2, ...)
const std::string kFiberBundleStlPath = "cad/interfaces/interface_bundle_fiber.stl";

// deterministic, order-independent grouping key for a material pair --
// e.g. "pmma|pmma_cl2" regardless of which side was materialA/materialB,
// so every instance of the same fiber LAYER PAIR (across however many
// individual fibers) accumulates into one shared stats row.
std::string MaterialPairKey(const std::string &matA, const std::string &matB) {
  std::string a = matA, b = matB;
  if (b < a)
    std::swap(a, b);
  return a + "|" + b;
}

// per-fiber-layer-pair timing + mesh-quality accumulator, populated in the
// main meshing loop below and summarized once it's done -- lets a run
// answer "which fiber layer pair is actually slow / actually broken?"
// directly instead of guessing from thousands of individual log lines.
struct FiberLayerStats {
  int count = 0;
  double total_ms = 0.0;
  long open_edges = 0;
  long nonmanifold_edges = 0;
  long degenerate_tris = 0;
  int mesh_empty = 0;
  double total_area_mm2 = 0.0;
};

} // namespace

// ============================================================
// mesh all optical interfaces
// ============================================================

void SurfaceMesher::MeshInterfaces(

    DetectorAssembly &assembly) {

  std::cout << "\nMeshing optical interfaces...\n" << std::endl;

  // --------------------------------------------------------
  // create output directories
  // --------------------------------------------------------

  fs::create_directories("cad/interfaces");
  fs::create_directories("metadata");

  // --------------------------------------------------------
  // loop over interfaces (parallel).
  //
  // Each iteration meshes a DISTINCT iface.boundary and writes
  // its own file — no cross-iteration shared state except the
  // report stream, which is guarded. iface.n_triangles/area_mm2
  // are written per-element, so those are race-free too.
  //
  // schedule(dynamic): interface sizes are very uneven (the
  // LAr-facing boundaries dwarf the small copper patches), so
  // dynamic hand-out keeps all threads busy instead of leaving
  // one thread stuck on the heavy ones under static chunking.
  // --------------------------------------------------------

  const int n = static_cast<int>(assembly.interfaces.size());

  // Written incrementally after each interface completes (see below) so
  // metadata/interfaces.json can be polled mid-run to see audit stats for
  // interfaces finished so far, instead of only appearing once the whole
  // (potentially very long) meshing pass is done.
  InterfaceExtractor extractor;

  // Accumulates every fiber-layer interface's boundary into ONE shared
  // compound (see IsFiberLayerMaterial/kFiberBundleStlPath above), guarded
  // by the same critical section used to build it below. Written out
  // once, after the main loop, instead of per-interface.
  TopoDS_Compound fiber_bundle;
  BRep_Builder fiber_bundle_builder;
  fiber_bundle_builder.MakeCompound(fiber_bundle);
  bool has_fiber_bundle = false;

  // per-fiber-layer-pair timing + mesh-quality breakdown (see
  // FiberLayerStats/MaterialPairKey above) -- printed as a summary table
  // once the main loop below is done.
  std::map<std::string, FiberLayerStats> fiber_layer_stats;

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; ++i) {

    auto &iface = assembly.interfaces[i];
    auto iface_t0 = std::chrono::steady_clock::now();

    // ----------------------------------------------------
    // verify boundary topology
    // ----------------------------------------------------

    TopAbs_ShapeEnum shape_type = iface.boundary.ShapeType();

    // ----------------------------------------------------
    // reject volumetric intersections
    // ----------------------------------------------------

    if (shape_type == TopAbs_SOLID || shape_type == TopAbs_COMPSOLID) {

#pragma omp critical(mesher_cout)
      std::cout << "Skipping volumetric intersection:\n"
                << "  interface " << iface.id << std::endl;

      continue;
    }

    // ----------------------------------------------------
    // tessellate OCC boundary surface
    //
    // Standard_False for the parallel flag: the outer OpenMP loop
    // owns the cores. Leaving BRepMesh's own in-parallel flag on
    // would nest parallelism and oversubscribe.
    // ----------------------------------------------------

    BRepMesh_IncrementalMesh mesher(
        iface.boundary,
        0.1,           // linear deflection (relative fraction when isRelative)
        Standard_True, // isRelative — scale tolerance to each face's size
        0.5,           // angular deflection in radians (~28°)
        Standard_False // NOT in parallel — outer OpenMP loop owns the cores
    );

    mesher.Perform();

    // ----------------------------------------------------
    // export raw OCC BREP (debug only)
    // ----------------------------------------------------

    std::string brep_name =
        "cad/interfaces/interface_" + std::to_string(iface.id) + ".brep";

    if (kWriteBrep) {
      BRepTools::Write(iface.boundary, brep_name.c_str());
    }

    // ----------------------------------------------------
    // export triangulated STL (binary)
    // ----------------------------------------------------

    std::string stl_name =
        "cad/interfaces/interface_" + std::to_string(iface.id) + ".stl";

    bool bundled = IsFiberLayerMaterial(iface.materialA) ||
                   IsFiberLayerMaterial(iface.materialB);

    if (bundled) {
      // Accumulate into the ONE shared fiber compound instead of writing
      // this interface's own file -- area/mesh-quality below are computed
      // from iface.boundary's own triangulation regardless, so nothing
      // about this interface's own metadata is affected.
      stl_name = kFiberBundleStlPath;
      iface.stl_override = stl_name;
#pragma omp critical(fiber_bundle)
      {
        fiber_bundle_builder.Add(fiber_bundle, iface.boundary);
        has_fiber_bundle = true;
      }
    } else {
      StlAPI_Writer stl_writer;
      stl_writer.ASCIIMode() = Standard_False; // write binary STL
      stl_writer.Write(iface.boundary, stl_name.c_str());
    }

    size_t total_vertices = 0;
    size_t total_triangles = 0;
    size_t total_faces = 0;
    double mesh_area_mm2 = 0.0;
    int degenerate_tris = 0;
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edge_count;

    // ----------------------------------------------------
    // iterate surface faces
    // ----------------------------------------------------

    for (TopExp_Explorer explorer(iface.boundary, TopAbs_FACE); explorer.More();
         explorer.Next()) {

      total_faces++;

      TopoDS_Face face = TopoDS::Face(explorer.Current());

      TopLoc_Location location;

      Handle(Poly_Triangulation) triangulation =
          BRep_Tool::Triangulation(face, location);

      if (triangulation.IsNull())
        continue;

      total_vertices += triangulation->NbNodes();
      total_triangles += triangulation->NbTriangles();

      const gp_Trsf &trsf = location.Transformation();

      // Faces are independently triangulated (their own local node
      // indexing), so edges shared across faces can only be matched
      // by transformed 3D position, not by node index.
      for (int t = 1; t <= triangulation->NbTriangles(); ++t) {
        int i1, i2, i3;
        triangulation->Triangle(t).Get(i1, i2, i3);

        gp_Pnt p1 = triangulation->Node(i1).Transformed(trsf);
        gp_Pnt p2 = triangulation->Node(i2).Transformed(trsf);
        gp_Pnt p3 = triangulation->Node(i3).Transformed(trsf);

        gp_Vec e1(p1, p2), e2(p1, p3);
        double tri_area = 0.5 * e1.Crossed(e2).Magnitude();
        mesh_area_mm2 += tri_area;
        if (tri_area < kDegenerateAreaFloor)
          ++degenerate_tris;

        ++edge_count[MakeEdgeKey(p1, p2)];
        ++edge_count[MakeEdgeKey(p2, p3)];
        ++edge_count[MakeEdgeKey(p3, p1)];
      }
    }

    int open_edges = 0, nonmanifold_edges = 0;
    for (const auto &kv : edge_count) {
      if (kv.second == 1)
        ++open_edges;
      else if (kv.second > 2)
        ++nonmanifold_edges;
    }

    // ----------------------------------------------------
    // compute surface area (mm²)
    // ----------------------------------------------------

    GProp_GProps props;
    BRepGProp::SurfaceProperties(iface.boundary, props);
    double area_mm2 = props.Mass();

    // ----------------------------------------------------
    // store back onto the interface struct (this element only)
    // ----------------------------------------------------

    iface.n_triangles = (int)total_triangles;
    iface.area_mm2 = area_mm2;
    iface.mesh_area_mm2 = mesh_area_mm2;
    iface.occ_area_mm2 = area_mm2;
    iface.n_open_edges = open_edges;
    iface.n_nonmanifold_edges = nonmanifold_edges;
    iface.n_degenerate_tris = degenerate_tris;
    iface.mesh_empty = (total_triangles == 0);

    bool broken = iface.mesh_empty || open_edges > 0 || nonmanifold_edges > 0 ||
                 degenerate_tris > 0;

    if (bundled) {
      double elapsed_ms =
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - iface_t0)
              .count();
      std::string key = MaterialPairKey(iface.materialA, iface.materialB);
#pragma omp critical(fiber_layer_stats)
      {
        FiberLayerStats &s = fiber_layer_stats[key];
        ++s.count;
        s.total_ms += elapsed_ms;
        s.open_edges += open_edges;
        s.nonmanifold_edges += nonmanifold_edges;
        s.degenerate_tris += degenerate_tris;
        s.mesh_empty += iface.mesh_empty ? 1 : 0;
        s.total_area_mm2 += area_mm2;
      }
    }

    // ----------------------------------------------------
    // report
    // ----------------------------------------------------

#pragma omp critical(mesher_cout)
    {
      std::cout << "Interface " << iface.id << " : " << iface.materialA
                << " ↔ " << iface.materialB << "\n"
                << "  shape type = " << shape_type << "\n"
                << "  faces      = " << total_faces << "\n"
                << "  vertices   = " << total_vertices << "\n"
                << "  triangles  = " << total_triangles << "\n"
                << "  area_mm2   = " << area_mm2 << " (mesh: " << mesh_area_mm2
                << ")\n"
                << "  brep file  = " << (kWriteBrep ? brep_name : "(skipped)")
                << "\n"
                << "  stl file   = " << stl_name << "\n";
      if (broken) {
        std::cout << "  [WARN] mesh quality: open_edges=" << open_edges
                  << " nonmanifold_edges=" << nonmanifold_edges
                  << " degenerate_tris=" << degenerate_tris
                  << (iface.mesh_empty ? " (EMPTY MESH)" : "") << "\n";
      }
      std::cout << std::endl;

      // stream progress: rewrite interfaces.json now, not just at the end,
      // so it can be polled mid-run. Serialized by this same critical
      // section (concurrent writers to the same file would corrupt it);
      // reading other threads' not-yet-finished entries mid-update is a
      // harmless, self-correcting race for this progress-only purpose.
      extractor.WriteInterfacesJSON(assembly, ".");
    }
  }

  // Write the single accumulated fiber bundle exactly once, now that every
  // interface has contributed its boundary above.
  if (has_fiber_bundle) {
    StlAPI_Writer stl_writer;
    stl_writer.ASCIIMode() = Standard_False;
    stl_writer.Write(fiber_bundle, kFiberBundleStlPath.c_str());
    std::cout << "Wrote fiber-layer bundle: " << kFiberBundleStlPath << std::endl;
  }

  std::cout << "Finished interface meshing.\n" << std::endl;

  // --------------------------------------------------------
  // per-fiber-layer-pair timing + mesh-quality breakdown
  // --------------------------------------------------------

  if (!fiber_layer_stats.empty()) {
    std::cout << "Fiber-layer meshing breakdown (tessellation + quality "
                 "audit time only, not interface discovery):\n";
    std::cout << "  layer pair                            count   total_ms  "
                 "avg_ms   open_edges  nonmanifold  degenerate  empty  "
                 "total_area_mm2\n";
    double grand_total_ms = 0.0;
    int grand_total_count = 0;
    for (const auto &kv : fiber_layer_stats) {
      const FiberLayerStats &s = kv.second;
      double avg_ms = s.count > 0 ? s.total_ms / s.count : 0.0;
      std::cout << "  " << std::left << std::setw(38) << kv.first
                << std::right << std::setw(6) << s.count << std::setw(11)
                << std::fixed << std::setprecision(1) << s.total_ms
                << std::setw(9) << avg_ms << std::setw(12) << s.open_edges
                << std::setw(13) << s.nonmanifold_edges << std::setw(12)
                << s.degenerate_tris << std::setw(7) << s.mesh_empty
                << std::setw(16) << s.total_area_mm2 << "\n";
      grand_total_ms += s.total_ms;
      grand_total_count += s.count;
    }
    std::cout << "  TOTAL: " << grand_total_count << " fiber-layer interfaces, "
              << grand_total_ms << " ms combined meshing time\n"
              << std::endl;
  }

  // --------------------------------------------------------
  // mesh-quality summary
  // --------------------------------------------------------

  {
    int n_broken = 0;
    for (const auto &iface : assembly.interfaces) {
      if (iface.mesh_empty || iface.n_open_edges > 0 ||
          iface.n_nonmanifold_edges > 0 || iface.n_degenerate_tris > 0)
        ++n_broken;
    }
    if (n_broken > 0) {
      std::cout << "[WARN] " << n_broken << " / " << assembly.interfaces.size()
                << " interfaces have mesh-quality issues (open/non-manifold "
                   "edges, degenerate triangles, or empty mesh) — see "
                   "n_open_edges/n_nonmanifold_edges/n_degenerate_tris/"
                   "mesh_empty in interfaces.json\n"
                << std::endl;
    }
  }

  // final rewrite: every interface's stats are now populated (the
  // incremental writes above already kept interfaces.json current
  // throughout, this just guarantees the very last one reflects all of
  // them together).
  extractor.WriteInterfacesJSON(assembly, ".");
}
