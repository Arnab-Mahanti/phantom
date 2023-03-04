# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.12)

project ("phantom")

# Include sub-projects.
add_subdirectory(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/matplotplusplus)
add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/highFive)