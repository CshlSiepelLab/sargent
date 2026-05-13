# src/progs/CMakeLists.txt

# ---- sargent.help: regenerate from sargent.help_src when it changes ----
add_custom_command(
    OUTPUT  "${CMAKE_CURRENT_SOURCE_DIR}/sargent.help"
    COMMAND bash "${CMAKE_SOURCE_DIR}/src/munge-help.sh"
                 "${CMAKE_CURRENT_SOURCE_DIR}/sargent.help_src"
            > "${CMAKE_CURRENT_SOURCE_DIR}/sargent.help"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/sargent.help_src"
    COMMENT "Generating sargent.help from sargent.help_src"
    VERBATIM
)
add_custom_target(sargent_help
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/sargent.help"
)

# ---- executables ----
add_executable(sargent sargent.c)

add_dependencies(sargent sargent_help)

# Read version file and inject as compile definition
file(READ "${CMAKE_SOURCE_DIR}/version" SARGENT_VERSION_RAW)
string(STRIP "${SARGENT_VERSION_RAW}" SARGENT_VERSION)
message(STATUS "SARGENT_VERSION='${SARGENT_VERSION}'")

target_compile_definitions(sargent PRIVATE SARGENT_VERSION="${SARGENT_VERSION}")

target_compile_options(sargent PRIVATE -Wall -Wextra)
target_link_libraries(sargent PRIVATE sargent_core phast::phast_core)

install(TARGETS sargent RUNTIME DESTINATION bin)
