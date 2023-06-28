# include_directories(../include ${Boost_INCLUDE_DIR})
# set(SOURCES main.cpp)
# add_executable(triangle_group.x ${SOURCES})
# target_include_directories(group_algebra.x PUBLIC ../include)

# -- fuchsian generate group
set(SOURCES fuchsian_generate_group.cpp ../src/fuchsian_group.cpp ../src/fuchsian_ring.cpp ../include/fuchsian_ring.h ../include/fuchsian_group.h )
add_executable(fuchsian_generate_group ${SOURCES})
target_include_directories(fuchsian_generate_group PRIVATE ../include)

#-- fuchsian_generate_reg.cpp

include_directories(../include ${Boost_INCLUDE_DIR} ${Eigen3_INCLUDE_DIR})
set(SOURCES fuchsian_generate_reg.cpp ../src/fuchsian_group.cpp ../src/fuchsian_ring.cpp)
add_executable(fuchsian_generate_reg ${SOURCES})
target_link_libraries(fuchsian_generate_reg OpenMP::OpenMP_CXX)
target_compile_options(fuchsian_generate_reg PRIVATE -ftrapv)

#-- generate_group.cpp

include_directories(../include ${Boost_INCLUDE_DIR} ${Eigen3_INCLUDE_DIR})
set(SOURCES generate_group.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp)
add_executable(generate_group ${SOURCES})

#-- generate_reg.cpp

include_directories(../include ${Boost_INCLUDE_DIR} ${Eigen3_INCLUDE_DIR})
set(SOURCES generate_reg.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp)
add_executable(generate_reg ${SOURCES})
target_link_libraries(generate_reg OpenMP::OpenMP_CXX)
target_compile_options(generate_reg PRIVATE -ftrapv)
#target_include_directories(generate_laplacian PRIVATE ${OpenMP_INCLUDE_DIRS})

# include_directories(../include ${Boost_INCLUDE_DIR})
# set(SOURCES generate_group_open.cpp)
# add_executable(generate_group_open.x ${SOURCES})

add_library (m_npy npy.f90)

#set(OpenBLAS64_INCLUDE /home/fox/Git/OpenBLAS/include)
#set(OpenBLAS64_LIB /home/fox/Git/OpenBLAS/lib/libopenblas_haswellp-r0.3.21.dev.a)
#${OpenBLAS64_LIB} ${OpenBLAS64_LIBRARIES}
if(OpenBLAS64_FOUND)
  set(SOURCES dsyev_hamiltonian.f90)
  add_executable(dsyev_hamiltonian ${SOURCES})
  target_link_libraries(dsyev_hamiltonian  OpenBLAS64::OpenBLAS64  Threads::Threads m_npy)
  target_include_directories(dsyev_hamiltonian PRIVATE ${Threads_INCLUDE_DIRS} ${OpenBLAS64_INCLUDE_DIRS})
  target_compile_options(dsyev_hamiltonian PRIVATE -finteger-4-integer-8)
endif()
