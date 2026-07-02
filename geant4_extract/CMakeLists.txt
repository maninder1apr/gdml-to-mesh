cmake_minimum_required(VERSION 3.16)

project(geant4_extract LANGUAGES CXX)

# ============================================================
# C++ standard
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================
# dependencies
# ============================================================

find_package(Geant4 REQUIRED gdml)
find_package(OpenCASCADE REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

# ============================================================
# include directories
# ============================================================

include_directories(
    include
    ${Geant4_INCLUDE_DIRS}
    ${OpenCASCADE_INCLUDE_DIR}
)

# ============================================================
# cnpy
# ============================================================

add_library(cnpy
    cnpy/cnpy.cpp
)

target_include_directories(cnpy PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/cnpy
)

target_link_libraries(cnpy PUBLIC
    ZLIB::ZLIB
)

# ============================================================
# nlohmann_json
# ============================================================

include(FetchContent)

FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.2
)

FetchContent_MakeAvailable(json)

# ============================================================
# OCC mesher
# ============================================================

add_executable(occ_mesher
    src/main_occ_mesher.cc
    src/occ_mesher.cc
    src/SolidConverter.cc
    src/AssemblyBuilder.cc
    src/MaterialExporter.cc
    src/SurfaceExporter.cc
    src/InterfaceExtractor.cc
    src/SurfaceMesher.cc
    src/VolumeMetadataExporter.cc
    src/VolumeSTLExporter.cc
)

target_link_libraries(occ_mesher
    cnpy
    ${Geant4_LIBRARIES}
    ${OpenCASCADE_LIBRARIES}
    nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX
)

# ============================================================
# extract_mesh
# ============================================================

add_executable(extract_mesh
    extract_mesh.cc
)

target_link_libraries(extract_mesh
    cnpy
    ${Geant4_LIBRARIES}
    ${OpenCASCADE_LIBRARIES}
)

# ============================================================
# extract_placements
# ============================================================

add_executable(extract_placements
    extract_placements.cc
)

target_link_libraries(extract_placements
    cnpy
    ${Geant4_LIBRARIES}
    ${OpenCASCADE_LIBRARIES}
)

# ============================================================
# extract_surfaces
# ============================================================

add_executable(extract_surfaces
    extract_surfaces.cc
)

target_link_libraries(extract_surfaces
    cnpy
    ${Geant4_LIBRARIES}
    ${OpenCASCADE_LIBRARIES}
)