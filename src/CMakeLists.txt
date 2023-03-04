# CMakeList.txt : CMake project for phantom, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (phantom "phantom.cpp" "phantom.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET phantom PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
find_package(PythonLibs REQUIRED)

target_include_directories(phantom PUBLIC ${CMAKE_SOURCE_DIR}/dependencies)
target_include_directories(phantom PUBLIC ${Python3_INCLUDE_DIRS})

target_link_libraries(phantom PUBLIC matplot)

if (WIN32 AND MSVC)
	target_link_libraries(phantom ${Python3_LIBRARIES_RELEASE})
endif()
