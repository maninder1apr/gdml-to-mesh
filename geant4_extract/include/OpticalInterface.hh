#pragma once

#include <string>

#include <TopoDS_Shape.hxx>

struct OpticalInterface {

    int id;

    int volumeA;
    int volumeB;

    std::string nameA;
    std::string nameB;

    std::string materialA;
    std::string materialB;

    // filled by InterfaceExtractor::Extract
    std::string pv_inside;
    std::string pv_outside;
    std::string lv_inside;
    std::string lv_outside;
    std::string mat_inside;
    std::string mat_outside;
    bool        is_detector      = false;
    int         detector_channel = -1;  // GDML channel id, -1 if not a detector

    // normals point from pv_inside toward pv_outside
    TopoDS_Shape boundary;

    // filled by SurfaceMesher
    int    n_triangles = 0;
    double area_mm2    = 0.0;
};
