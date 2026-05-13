# src/lib/CMakeLists.txt

add_library(sargent_core STATIC
    bed.c
    trees.c
)

target_include_directories(sargent_core PUBLIC
    "${SARGENT_INCLUDE_DIR}"
)

target_link_libraries(sargent_core PUBLIC
    tskit::tskit
    m
)

if(TARGET phast::phast_core)
  target_link_libraries(sargent_core PUBLIC phast::phast_core)
endif()

find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(sargent_core PUBLIC OpenMP::OpenMP_C)
endif()
