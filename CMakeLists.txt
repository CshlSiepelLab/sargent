cmake_minimum_required(VERSION 3.16)

# Single source of truth for the version number: the top-level `version` file.
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/version" SARGENT_VERSION_FILE)
string(STRIP "${SARGENT_VERSION_FILE}" SARGENT_VERSION_FILE)

project(sargent
  VERSION ${SARGENT_VERSION_FILE}
  LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---- tskit C library (compiled in statically) ---------------------------
# TSKIT_SOURCE_DIR is the tskit C source tree: the c/ directory of the git
# submodule by default, or an unpacked tskit C release tarball
# (tskit-X.Y.Z.tar.xz from github.com/tskit-dev/tskit/releases), which has
# the same layout.
set(TSKIT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/tskit/c"
    CACHE PATH "tskit C source tree (submodule c/ or unpacked C release)")
if(NOT EXISTS "${TSKIT_SOURCE_DIR}/tskit.h")
  message(FATAL_ERROR
    "tskit sources not found in ${TSKIT_SOURCE_DIR}.  Clone with "
    "--recurse-submodules (or run `git submodule update --init`), or set "
    "-DTSKIT_SOURCE_DIR to an unpacked tskit C release.")
endif()

add_library(tskit STATIC
  ${TSKIT_SOURCE_DIR}/tskit/core.c
  ${TSKIT_SOURCE_DIR}/tskit/tables.c
  ${TSKIT_SOURCE_DIR}/tskit/trees.c
  ${TSKIT_SOURCE_DIR}/tskit/genotypes.c
  ${TSKIT_SOURCE_DIR}/tskit/stats.c
  ${TSKIT_SOURCE_DIR}/tskit/convert.c
  ${TSKIT_SOURCE_DIR}/tskit/haplotype_matching.c
  ${TSKIT_SOURCE_DIR}/subprojects/kastore/kastore.c)
target_include_directories(tskit SYSTEM PUBLIC
  ${TSKIT_SOURCE_DIR}
  ${TSKIT_SOURCE_DIR}/subprojects/kastore)
target_link_libraries(tskit PUBLIC m)
add_library(tskit::tskit ALIAS tskit)

file(READ "${TSKIT_SOURCE_DIR}/VERSION.txt" TSKIT_VERSION)
string(STRIP "${TSKIT_VERSION}" TSKIT_VERSION)
message(STATUS "tskit ${TSKIT_VERSION} from ${TSKIT_SOURCE_DIR}")

# ---- Project layout ----------------------------------------------------
# OpenMP detection happens inside src/lib/CMakeLists.txt via the standard
# find_package(OpenMP).  Building with Homebrew llvm clang (see README) gives
# OpenMP for free; Apple Clang users get a non-parallel build.
set(SARGENT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(src/lib)
add_subdirectory(src/progs)
add_subdirectory(test/lib)
