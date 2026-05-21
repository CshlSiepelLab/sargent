# src/lib/CMakeLists.txt

add_library(sargent_core STATIC
    bed.c
    trees.c
    variational.c
    conjgrad.c
    hutchinson.c
)

target_include_directories(sargent_core PUBLIC
    "${SARGENT_INCLUDE_DIR}"
)

target_link_libraries(sargent_core PUBLIC
    tskit::tskit
    m
)


find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(sargent_core PUBLIC OpenMP::OpenMP_C)
  target_compile_definitions(sargent_core PUBLIC SARGENT_HAS_OPENMP)
endif()
