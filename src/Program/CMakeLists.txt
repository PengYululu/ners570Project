# Root CMakeLists.txt (in /project directory)
cmake_minimum_required(VERSION 3.12)
project(NumericalModel Fortran)

# Set Fortran compiler flags
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -fbacktrace")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -traceback")
endif()

# Add source directory
add_subdirectory(src)

# src/CMakeLists.txt
# Create library for LibMath
add_library(LibMath STATIC
    LibMath.f90
)

# Create libraries for each module
add_library(ModCore STATIC
    ModParameters.f90
    ModParamReader.f90
    ModBlock.f90
    ModGhostCell.f90
    ModJobAssigner.f90
    ModInitiator.f90
    ModSolver.f90
    ModCommunicator.f90
    ModAdvance.f90
    ModWritePlot.f90
)

# Set dependencies between modules based on the diagram
set_source_files_properties(ModParamReader.f90 PROPERTIES DEPENDS ModParameters.f90)
set_source_files_properties(ModBlock.f90 PROPERTIES DEPENDS "ModParameters.f90;ModParamReader.f90")
set_source_files_properties(ModGhostCell.f90 PROPERTIES DEPENDS ModBlock.f90)
set_source_files_properties(ModSolver.f90 PROPERTIES DEPENDS "ModBlock.f90;LibMath.f90")
set_source_files_properties(ModCommunicator.f90 PROPERTIES DEPENDS ModParameters.f90)
set_source_files_properties(ModAdvance.f90 PROPERTIES DEPENDS ModCommunicator.f90)
set_source_files_properties(ModWritePlot.f90 PROPERTIES DEPENDS "ModInitiator.f90;ModAdvance.f90")

# Link dependencies
target_link_libraries(ModCore PUBLIC LibMath)

# Add program subdirectory
add_subdirectory(program)

# src/program/CMakeLists.txt
# Create the executable
add_executable(numerical_model
    main.f90
)

# Link with the libraries
target_link_libraries(numerical_model
    PRIVATE
    ModCore
    LibMath
)