cmake_minimum_required(VERSION 3.16)

project(sargent
  VERSION 0.1.0
  LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---- tskit dependency (built from submodule via meson) -----------------
include(ExternalProject)

set(TSKIT_SOURCE_DIR  "${CMAKE_SOURCE_DIR}/third_party/tskit/c")
set(TSKIT_BINARY_DIR  "${CMAKE_BINARY_DIR}/tskit-build")
set(TSKIT_INSTALL_DIR "${CMAKE_BINARY_DIR}/tskit-install")
set(TSKIT_LIB         "${TSKIT_INSTALL_DIR}/lib/libtskit.dylib")

find_program(MESON_EXECUTABLE meson REQUIRED)
find_program(NINJA_EXECUTABLE ninja REQUIRED)

ExternalProject_Add(tskit_external
  SOURCE_DIR  ${TSKIT_SOURCE_DIR}
  BINARY_DIR  ${TSKIT_BINARY_DIR}
  INSTALL_DIR ${TSKIT_INSTALL_DIR}
  PATCH_COMMAND
    ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/third_party/patches/tskit_meson_options.txt
      ${TSKIT_SOURCE_DIR}/meson_options.txt
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/third_party/patches/tskit_meson.build
      ${TSKIT_SOURCE_DIR}/meson.build
  CONFIGURE_COMMAND
    ${MESON_EXECUTABLE} setup ${TSKIT_BINARY_DIR} ${TSKIT_SOURCE_DIR}
      --prefix=${TSKIT_INSTALL_DIR}
      --wipe
      -Dbuild_tests=false
      -Dbuild_examples=false
  BUILD_COMMAND
    ${NINJA_EXECUTABLE} -C ${TSKIT_BINARY_DIR} libtskit.dylib
  INSTALL_COMMAND
    ${NINJA_EXECUTABLE} -C ${TSKIT_BINARY_DIR} install
    COMMAND ${CMAKE_COMMAND} -E copy
      ${TSKIT_SOURCE_DIR}/subprojects/kastore/kastore.h
      ${TSKIT_INSTALL_DIR}/include/kastore.h
  BUILD_BYPRODUCTS ${TSKIT_LIB}
)

file(MAKE_DIRECTORY "${TSKIT_INSTALL_DIR}/include")

add_library(tskit::tskit IMPORTED SHARED GLOBAL)
add_dependencies(tskit::tskit tskit_external)
set_target_properties(tskit::tskit PROPERTIES
  IMPORTED_LOCATION     ${TSKIT_LIB}
  INTERFACE_INCLUDE_DIRECTORIES ${TSKIT_INSTALL_DIR}/include
  INTERFACE_LINK_LIBRARIES m)

message(STATUS "tskit will be built from ${TSKIT_SOURCE_DIR}")

# ---- Project layout ----------------------------------------------------
# OpenMP detection happens inside src/lib/CMakeLists.txt via the standard
# find_package(OpenMP).  Building with Homebrew llvm clang (see README) gives
# OpenMP for free; Apple Clang users get a non-parallel build.
set(SARGENT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(src/lib)
add_subdirectory(src/progs)
add_subdirectory(test/lib)

# ---- Installation ------------------------------------------------------
install(DIRECTORY "${SARGENT_INCLUDE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
