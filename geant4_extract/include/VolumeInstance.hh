#pragma once

#include <TopoDS_Shape.hxx>

#include <string>
#include <cstdint>

struct VolumeInstance {

    // sentinel for "no mother" (top-level volume directly under the World)
    static constexpr uint64_t kNoMother = ~0ull;

    uint64_t id;

    // id of the mother VolumeInstance, or kNoMother if the mother is the
    // (excluded) World volume
    uint64_t mother_id = kNoMother;

    std::string name;

    std::string lv_name;

    std::string material;

    TopoDS_Shape shape;
};
