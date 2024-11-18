# Glob all Fortran source files in the directory
file(GLOB ALL_SOURCES *.f90)

# Exclude ModCommunicator.f90 (requires MPI)
list(FILTER ALL_SOURCES EXCLUDE REGEX ".*ModCommunicator.f90")

# Find MPI
find_package(MPI REQUIRED)

# MPI source
set(MPI_SOURCES ModCommunicator.f90)

# Non-MPI sources (remaining files)
set(NON_MPI_SOURCES ${ALL_SOURCES})

# Create a library for non-MPI sources
add_library(NonMPI_Library ${NON_MPI_SOURCES})

# Create a library for the MPI source
add_library(MPI_Library ${MPI_SOURCES})
target_compile_options(MPI_Library PRIVATE ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(MPI_Library PRIVATE ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(MPI_Library PRIVATE ${MPI_Fortran_LIBRARIES})
set_target_properties(MPI_Library PROPERTIES Fortran_COMPILER mpif90)

# Combine both libraries into the final executable
#add_executable(Simulation main.f90)
#target_link_libraries(Simulation PRIVATE MPI_Library NonMPI_Library)
