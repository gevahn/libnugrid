cmake_minimum_required (VERSION 3.5)
project (libNuGrid)
# The version number.
set (Tutorial_VERSION_MAJOR 0)
set (Tutorial_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "-std=c++11")

 
# configure a header file to pass some of the CMake settings
# to the source code
# configure_file (
#   )
 
include_directories("${PROJECT_BINARY_DIR}")


add_library(nugrid nugrid.cpp lebedev_grids/lebedevData.cpp) 

install (TARGETS nugrid DESTINATION ${PROJECT_SOURCE_DIR}/tests/lib)
install (FILES nugrid.h DESTINATION ${PROJECT_SOURCE_DIR}/tests/include)
