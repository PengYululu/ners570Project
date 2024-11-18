# Add all Fortran source files in this directory
file(GLOB SRC_FILES *.f90)

# Define the executable
add_executable(Simulation ${SRC_FILES})

# Link the Math library
target_link_libraries(Simulation PRIVATE Math)

# Optionally specify include directories if needed
# target_include_directories(Simulation PUBLIC ../include/)
