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

if(APPLE)
     set(CMAKE_CXX_FLAGS "-Wall -ld_classic")
     set(CMAKE_CXX_FLAGS_DEBUG "-g")
     set(CMAKE_CXX_FLAGS_RELEASE "-O3")
else()
     set(CMAKE_CXX_FLAGS "-Wall")
     set(CMAKE_CXX_FLAGS_DEBUG "-g")
     set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()
#---------------------------------------------------------------------#
# Copy input files to both ./build and ./build/bin 
#---------------------------------------------------------------------#

file(COPY ${CMAKE_SOURCE_DIR}/wls/ring_reduction.inp
     DESTINATION ${CMAKE_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/wls/generators.tar.gz
     DESTINATION ${CMAKE_BINARY_DIR})

add_custom_target( UntarGenerators ALL)
add_custom_command(TARGET UntarGenerators PRE_BUILD
     COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_BINARY_DIR}/generators.tar.gz
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     DEPENDS ${CMAKE_BINARY_DIR}/generators.tar.gz
     COMMENT "Extracting generators.tar.gz"
     VERBATIM
)

#---------------------------------------------------------------------#
# Directory setup
#---------------------------------------------------------------------#

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(INCLUDE_PATH ${CMAKE_BINARY_DIR}/include) 


#---------------------------------------------------------------------#
# Find the C++ boost library
#---------------------------------------------------------------------#

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenBLAS64 QUIET)

message((${Boost_INCLUDE_DIR}))

#---------------------------------------------------------------------#
# Source directories
#---------------------------------------------------------------------#

add_subdirectory(src)
add_subdirectory(tests)