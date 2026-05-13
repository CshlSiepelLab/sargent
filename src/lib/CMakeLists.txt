# src/lib/CMakeLists.txt

add_library(sargent_core INTERFACE)

target_include_directories(sargent_core INTERFACE
    "${SARGENT_INCLUDE_DIR}"
)

target_link_libraries(sargent_core INTERFACE
    tskit::tskit
)

if(TARGET phast::phast_core)
  target_link_libraries(sargent_core INTERFACE phast::phast_core)
endif()

find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(sargent_core INTERFACE OpenMP::OpenMP_C)
endif()
