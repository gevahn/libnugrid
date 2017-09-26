cmake_minimum_required (VERSION 3.5)
project (libNuGrid)

set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories("${PROJECT_BINARY_DIR}")

add_library(nugrid SHARED nugrid.cpp lebedev_grids/lebedevData.cpp) 

install (TARGETS nugrid DESTINATION ${PROJECT_SOURCE_DIR}/tests/lib)
install (FILES nugrid.h DESTINATION ${PROJECT_SOURCE_DIR}/tests/include)


# compile a pybind11 module of libnugrid

if("$ENV{CONDA_PREFIX}" STREQUAL "" )

    message(FATAL_ERROR "This must be run within the MolSSI SSS conda prefix. 
                       Delete the build directory and try again")

endif()
 

set(PYTHON_EXECUTABLE $ENV{CONDA_PREFIX}/bin/python3)

# Find the pybind11 in the conda path

set(PYBIND11_CPP_STANDARD -std=c++11)

find_package(pybind11 CONFIG REQUIRED

    PATHS $ENV{CONDA_PREFIX}

             NO_DEFAULT_PATH)

         message(STATUS "Found pybind11: ${pybind11_CONFIG}" )

# Creates a python module named module_name

pybind11_add_module(pynugrid MODULE pynugrid.cpp)
target_link_libraries(pynugrid nugrid)
