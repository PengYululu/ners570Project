# Add all Fortran source files in this directory to the library
file(GLOB SRC_MATH *.f90)

# Create a static library called "Math"
add_library(Math STATIC ${SRC_MATH})

# Specify the language for the Math library
set_target_properties(Math PROPERTIES LINKER_LANGUAGE Fortran)