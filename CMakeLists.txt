# cmake_minimum_required(VERSION 3.15)
# project(SimulationCode LANGUAGES Fortran)

# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

# # Define the output directories for binaries and libraries
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# # Add subdirectories
# add_subdirectory(srcMath)
# add_subdirectory(src)
# add_subdirectory(tests)

# enable_testing()


# Top of CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(Simulation LANGUAGES Fortran)
set(CMAKE_Fortran_COMPILER mpif90)

#update MPI setup
find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS} -O3")

# Define the output directories for binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Add subdirectories
add_subdirectory(srcMath)
add_subdirectory(src)
add_subdirectory(tests)

enable_testing()