# src/lib/CMakeLists.txt

add_library(sargent_core STATIC
    bed.c
    trees.c
    variational.c
    conjgrad.c
    hutchinson.c
    gall.c
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

# LAPACK / BLAS for eigendecomposition (gall.c: dsyev_)
if(APPLE)
  # Accelerate framework ships with every macOS install
  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  target_link_libraries(sargent_core PUBLIC ${ACCELERATE_FRAMEWORK})
  target_compile_definitions(sargent_core PUBLIC SARGENT_HAS_LAPACK)
else()
  find_package(LAPACK)
  if(LAPACK_FOUND)
    target_link_libraries(sargent_core PUBLIC ${LAPACK_LIBRARIES})
    target_compile_definitions(sargent_core PUBLIC SARGENT_HAS_LAPACK)
    message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
  else()
    message(WARNING "LAPACK not found — gall precompute and EMMA mode disabled")
  endif()
endif()
