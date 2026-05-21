# test/lib/CMakeLists.txt — Phase 5 paranoia harnesses
#
# Built alongside the main sargent target; not installed.

add_executable(check_tree_matvec check_tree_matvec.c)
target_compile_options(check_tree_matvec PRIVATE -Wall -Wextra)
target_link_libraries(check_tree_matvec PRIVATE sargent_core)
