# src/progs/CMakeLists.txt

# ---- sargent.help.h: generate from sargent_help.txt at configure time ----
# Re-run cmake when the text file changes to regenerate the header.
set(_help_txt "${CMAKE_CURRENT_SOURCE_DIR}/sargent_help.txt")
set(_help_hdr "${CMAKE_CURRENT_SOURCE_DIR}/sargent.help.h")

file(READ "${_help_txt}" _help_raw)
string(REPLACE "\\" "\\\\" _help_esc "${_help_raw}")
string(REPLACE "\"" "\\\"" _help_esc "${_help_esc}")
string(REPLACE "\n" "\\n\"\n\"" _help_esc "${_help_esc}")
file(WRITE "${_help_hdr}"
    "/* Auto-generated from sargent_help.txt — do not edit directly. */\n"
    "static const char HELP[] =\n"
    "\"${_help_esc}\";\n")

# ---- executables ----
add_executable(sargent sargent.c)

# Read version file and inject as compile definition
file(READ "${CMAKE_SOURCE_DIR}/version" SARGENT_VERSION_RAW)
string(STRIP "${SARGENT_VERSION_RAW}" SARGENT_VERSION)
message(STATUS "SARGENT_VERSION='${SARGENT_VERSION}'")

target_compile_definitions(sargent PRIVATE SARGENT_VERSION="${SARGENT_VERSION}")

target_compile_options(sargent PRIVATE -Wall -Wextra)
target_link_libraries(sargent PRIVATE sargent_core)

install(TARGETS sargent RUNTIME DESTINATION bin)
