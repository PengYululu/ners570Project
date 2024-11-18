# Enable CTest for test automation
enable_testing()

# Find MPI for tests that require it
find_package(MPI REQUIRED)

# Include the src directory for module files
include_directories(${CMAKE_BINARY_DIR}/src)

# Define the test source files
set(TEST_SOURCES
    test_JobAssigner.f90
)

# Loop through the test files
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    if(${TEST_NAME} STREQUAL "test_JobAssigner")
        # Create the executable for the MPI test
        add_executable(${TEST_NAME} ${TEST_SOURCE})

        # Add MPI compile flags, include paths, and libraries
        target_compile_options(${TEST_NAME} PRIVATE ${MPI_Fortran_COMPILE_FLAGS})
        target_include_directories(${TEST_NAME} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
        target_link_libraries(${TEST_NAME} PRIVATE MPI_Library NonMPI_Library ${MPI_Fortran_LIBRARIES})

        # Add the test to CTest
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    endif()
endforeach()