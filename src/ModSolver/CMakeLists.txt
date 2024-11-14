# src/ModSolver/CMakeLists.txt
add_library(ModSolver ModSolver.cpp)

# Link dependencies (e.g., LibMath)
target_link_libraries(ModSolver PRIVATE LibMath)
