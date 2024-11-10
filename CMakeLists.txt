cmake_minimum_required(VERSION 3.10)
project(NumericalModelSoftware)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Add all modules
add_subdirectory(src/ModParameters)
add_subdirectory(src/ModParamReader)
add_subdirectory(src/ModJobAssigner)
add_subdirectory(src/ModBlock)
add_subdirectory(src/ModGhostCell)
add_subdirectory(src/ModSolver)
add_subdirectory(src/LibMath)
add_subdirectory(src/ModCommunicator)
add_subdirectory(src/ModAdvance)
add_subdirectory(src/ModWritePlot)
add_subdirectory(src/ModInitiator)
add_subdirectory(src/Program)
