# CMakeList.txt : CMake project for phantom, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (phantom "phantom.cpp" "phantom.h")

set_property(TARGET phantom PROPERTY CXX_STANDARD 20)

# TODO: Add tests and install targets if needed.
target_include_directories(phantom PUBLIC ${CMAKE_SOURCE_DIR}/dependencies)

target_link_libraries(phantom PUBLIC matplot HighFive)


