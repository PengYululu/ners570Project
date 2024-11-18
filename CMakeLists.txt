cmake_minimum_required(VERSION 3.15)
project(SimulationCode LANGUAGES Fortran)

# Define the output directories for binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Add subdirectories
add_subdirectory(srcMath)
add_subdirectory(src)