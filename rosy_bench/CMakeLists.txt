# Headless build of the RoSy (orientation field) part of
# "Robust Hex-Dominant Mesh Generation using Field-Guided Polyhedral Agglomeration".
#
# Builds a single command-line tool, `rosy_field`, that runs
#   load OBJ -> TetGen tet meshing -> multi-resolution hierarchy -> RoSy (quaternion) smoothing
# and exports the frame field + singularity graph. No GUI / OpenGL / nanogui is built.
#
#   cmake -S rosy_bench -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
cmake_minimum_required(VERSION 3.16)

project(rosy_bench CXX C)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(RHDM_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(RHDM_EIGEN_DIR ${RHDM_ROOT}/ext/nanoguig/ext/eigen)

if(NOT EXISTS ${RHDM_EIGEN_DIR}/Eigen/Core OR NOT EXISTS ${RHDM_ROOT}/ext/pcg32/pcg32.h)
  message(FATAL_ERROR "Submodules missing: run `git submodule update --init --recursive` in ${RHDM_ROOT}")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------- TBB
# An installed oneTBB (e.g. `brew install tbb`, `apt install libtbb-dev`).
# The TBB 2017 submodule in ext/tbb cannot be used: its CMakeLists forces policies
# that CMake >= 4 rejects, and it has no arm64-macOS machine port. The RoSy smoother
# only uses parallel_for/blocked_range/spin_mutex, so results do not depend on the TBB version.
find_package(TBB REQUIRED)
set(RHDM_TBB_TARGET TBB::tbb)

# ---------------------------------------------------------------- TetGen
# Compiled here (instead of ext/tetgen/CMakeLists.txt) so that Shewchuk's exact
# predicates are built without FMA contraction / optimisation, as TetGen requires.
add_library(rhdm_tetgen STATIC
  ${RHDM_ROOT}/ext/tetgen/tetgen.cxx
  ${RHDM_ROOT}/ext/tetgen/predicates.cxx)
target_compile_definitions(rhdm_tetgen PUBLIC TETLIBRARY)
target_include_directories(rhdm_tetgen PUBLIC ${RHDM_ROOT}/ext/tetgen)
set_source_files_properties(${RHDM_ROOT}/ext/tetgen/predicates.cxx
  PROPERTIES COMPILE_OPTIONS "-O0;-ffp-contract=off")
target_compile_options(rhdm_tetgen PRIVATE -w)

# ---------------------------------------------------------------- RoSy core
# Only the translation units reachable from load/tet_meshing/build/smoothOrientationsTet.
add_library(rhdm_core STATIC
  ${RHDM_ROOT}/src/hierarchy.cpp
  ${RHDM_ROOT}/src/orientations.cpp
  ${RHDM_ROOT}/src/meshio.cpp
  ${RHDM_ROOT}/src/orient_triangle_mesh.cpp
  ${RHDM_ROOT}/src/dedge.cpp
  ${RHDM_ROOT}/src/subdivide.cpp
  ${RHDM_ROOT}/ext/rply/rply.c)
# SINGLE_PRECISION matches the original build (Float = float).
target_compile_definitions(rhdm_core PUBLIC RHDM_HEADLESS SINGLE_PRECISION TETLIBRARY)
target_include_directories(rhdm_core PUBLIC
  ${RHDM_ROOT}/src
  ${RHDM_EIGEN_DIR}
  ${RHDM_ROOT}/ext/pcg32
  ${RHDM_ROOT}/ext/rply)
target_link_libraries(rhdm_core PUBLIC ${RHDM_TBB_TARGET} rhdm_tetgen)
# Old Eigen 3.2 / original sources are noisy with a modern clang.
target_compile_options(rhdm_core PUBLIC -Wno-deprecated-declarations -Wno-deprecated
  -Wno-unused-parameter -Wno-ignored-attributes -Wno-unknown-warning-option)

# ---------------------------------------------------------------- tool
add_executable(rosy_field ${CMAKE_CURRENT_SOURCE_DIR}/rosy_field.cpp)
target_link_libraries(rosy_field PRIVATE rhdm_core)
