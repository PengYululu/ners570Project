# Glob all Fortran source files in the directory
file(GLOB ALL_SOURCES *.f90)

# Exclude ModCommunicator.f90 and ModAdvance.f90 (requires MPI)
list(FILTER ALL_SOURCES EXCLUDE REGEX ".*(main|ModCommunicator|ModAdvance)\\.f90")

# Find MPI
find_package(MPI REQUIRED)

# Define MPI sources
set(MPI_SOURCES 
    ModCommunicator.f90
    ModAdvance.f90
)

# Non-MPI sources (remaining files)
set(NON_MPI_SOURCES ${ALL_SOURCES})

# Create a library for non-MPI sources
add_library(NonMPI_Library ${NON_MPI_SOURCES})
target_include_directories(NonMPI_Library PRIVATE ${CMAKE_BINARY_DIR}/srcMath)
target_link_libraries(NonMPI_Library PRIVATE Math)

# Create a library for MPI sources
add_library(MPI_Library ${MPI_SOURCES})
target_compile_options(MPI_Library PRIVATE ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(MPI_Library PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_BINARY_DIR}/srcMath)
target_link_libraries(MPI_Library PRIVATE Math ${MPI_Fortran_LIBRARIES})
set_target_properties(MPI_Library PROPERTIES Fortran_COMPILER mpif90)

# Add the main executable
add_executable(Simulation main.f90)
target_include_directories(Simulation PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${CMAKE_BINARY_DIR}/srcMath)
target_link_libraries(Simulation PRIVATE MPI_Library NonMPI_Library ${MPI_Fortran_LIBRARIES})
target_compile_options(Simulation PRIVATE ${MPI_Fortran_COMPILE_FLAGS})
set_target_properties(Simulation PROPERTIES Fortran_COMPILER mpif90)



# Define the run directory
set(RUN_DIR ${CMAKE_BINARY_DIR}/run)

# Create the run directory and copy input files
file(MAKE_DIRECTORY ${RUN_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/input/PARAM.in.test DESTINATION ${RUN_DIR})

# Copy the executable to the run directory after it's built
add_custom_command(
    TARGET Simulation POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Simulation> ${RUN_DIR}
)