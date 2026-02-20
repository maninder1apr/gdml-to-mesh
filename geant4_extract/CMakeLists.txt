cmake_minimum_required(VERSION 3.16)
project(geant4_extract LANGUAGES CXX)

# ------------------------------------------------------------
# C++ standard
# ------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------
# Geant4
# ------------------------------------------------------------
find_package(Geant4 REQUIRED gdml)

# ------------------------------------------------------------
# zlib (for cnpy)
# ------------------------------------------------------------
find_package(ZLIB REQUIRED)

# ------------------------------------------------------------
# cnpy (local)
# ------------------------------------------------------------
add_library(cnpy
    cnpy/cnpy.cpp
)

target_include_directories(cnpy PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/cnpy
)

target_link_libraries(cnpy PUBLIC
    ZLIB::ZLIB
)

# ------------------------------------------------------------
# extract_mesh
# ------------------------------------------------------------
add_executable(extract_mesh
    extract_mesh.cc
)

target_link_libraries(extract_mesh
    cnpy
    ${Geant4_LIBRARIES}
)

# ------------------------------------------------------------
# extract_placements
# ------------------------------------------------------------
add_executable(extract_placements
    extract_placements.cc
)

target_link_libraries(extract_placements
    cnpy
    ${Geant4_LIBRARIES}
)