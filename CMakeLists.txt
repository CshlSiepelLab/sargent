cmake_minimum_required(VERSION 3.16)

project(sargent
  VERSION 0.1.0
  LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# allows use of upper case PHAST_ROOT / TSKIT_ROOT
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

# ---- PHAST dependency --------------------------------------------------
set(PHAST_ROOT "" CACHE PATH "PHAST install prefix (contains lib/cmake/phast)")

if(PHAST_ROOT)
  list(PREPEND CMAKE_PREFIX_PATH "${PHAST_ROOT}")
  find_package(phast CONFIG REQUIRED)
else()
  message(STATUS "PHAST_ROOT not set — phast support disabled")
endif()

# ---- tskit dependency --------------------------------------------------
set(TSKIT_ROOT "/usr/local/Caskroom/miniforge/base/envs/sargent"
    CACHE PATH "tskit install prefix (conda env or other prefix)")

find_path(TSKIT_INCLUDE_DIR tskit.h
  PATHS "${TSKIT_ROOT}/include"
  REQUIRED)

find_library(TSKIT_LIB tskit
  PATHS "${TSKIT_ROOT}/lib"
  REQUIRED)

add_library(tskit::tskit IMPORTED INTERFACE)
set_target_properties(tskit::tskit PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${TSKIT_INCLUDE_DIR}"
  INTERFACE_LINK_LIBRARIES      "${TSKIT_LIB};m")

message(STATUS "tskit include: ${TSKIT_INCLUDE_DIR}")
message(STATUS "tskit lib:     ${TSKIT_LIB}")

# ---- Project layout ----------------------------------------------------
set(SARGENT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(src/lib)
add_subdirectory(src/progs)

# ---- Installation ------------------------------------------------------
install(DIRECTORY "${SARGENT_INCLUDE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
