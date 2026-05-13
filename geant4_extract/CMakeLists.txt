cmake_minimum_required(VERSION 3.16)

project(geant4_extract LANGUAGES CXX)

# ============================================================
# C++ standard
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================
# Geant4
# ============================================================

find_package(Geant4 REQUIRED gdml)

# ============================================================
# OpenCASCADE
# ============================================================

find_package(OpenCASCADE REQUIRED)

# ============================================================
# zlib (for cnpy)
# ============================================================

find_package(ZLIB REQUIRED)

# ============================================================
# include directories
# ============================================================

include_directories(
    include
    ${Geant4_INCLUDE_DIRS}
    ${OpenCASCADE_INCLUDE_DIRS}
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
# OCC mesher
# ============================================================

add_executable(occ_mesher
    src/main_occ_mesher.cc
    src/occ_mesher.cc
    src/SolidConverter.cc
    src/AssemblyBuilder.cc
    src/MaterialExporter.cc
    src/InterfaceExtractor.cc
    src/SurfaceMesher.cc
)

target_include_directories(occ_mesher PUBLIC
    include
)

target_link_libraries(occ_mesher
    cnpy
    ${Geant4_LIBRARIES}
    ${OpenCASCADE_LIBRARIES}
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
    ${OpenCASCADE_LIBRARIES}       # optional (if mesh functions use OCC)
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
    ${OpenCASCADE_LIBRARIES}       # optional (only needed if OCC used)
)

# ============================================================
# extract_surfaces (NEEDS OpenCascade)
# ============================================================

add_executable(extract_surfaces
    extract_surfaces.cc
)

target_link_libraries(extract_surfaces
    cnpy
    ${Geant4_LIBRARIES}
    ${OpenCASCADE_LIBRARIES}
)
