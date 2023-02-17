cmake_minimum_required(VERSION 2.9 FATAL_ERROR)
project(triangle_group)

enable_testing()

#---------------------------------------------------------------------#
# Compiler setup
#---------------------------------------------------------------------#

enable_language(CXX)
enable_language(Fortran)
set (CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


#---------------------------------------------------------------------#
# Copy input files to both ./build and ./build/bin 
#---------------------------------------------------------------------#

file(COPY ${CMAKE_SOURCE_DIR}/wls/ring_reduction.inp
     DESTINATION ${CMAKE_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/group_specs.inp
     DESTINATION ${CMAKE_BINARY_DIR}/bin)

file(COPY ${CMAKE_SOURCE_DIR}/wls/generators
DESTINATION ${CMAKE_BINARY_DIR})

#---------------------------------------------------------------------#
# Directory setup
#---------------------------------------------------------------------#

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(INCLUDE_PATH ${CMAKE_BINARY_DIR}/include) 


#---------------------------------------------------------------------#
# Find the C++ boost library
#---------------------------------------------------------------------#

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

#---------------------------------------------------------------------#
# Source directories
#---------------------------------------------------------------------#

add_subdirectory(src)
add_subdirectory(tests)

