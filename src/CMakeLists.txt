# include_directories(../include ${Boost_INCLUDE_DIR})
# set(SOURCES main.cpp)
# add_executable(triangle_group.x ${SOURCES})
# target_include_directories(group_algebra.x PUBLIC ../include)

#-- generate_group.cpp

include_directories(../include ${Boost_INCLUDE_DIR} ${Eigen3_INCLUDE_DIR})
set(SOURCES generate_group.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp)
add_executable(generate_group ${SOURCES})

#-- generate_laplacian.cpp

include_directories(../include ${Boost_INCLUDE_DIR} ${Eigen3_INCLUDE_DIR})
set(SOURCES generate_laplacian.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp)
add_executable(generate_laplacian ${SOURCES})
target_link_libraries(generate_laplacian OpenMP::OpenMP_CXX)
#target_include_directories(generate_laplacian PRIVATE ${OpenMP_INCLUDE_DIRS})

# include_directories(../include ${Boost_INCLUDE_DIR})
# set(SOURCES generate_group_open.cpp)
# add_executable(generate_group_open.x ${SOURCES})

add_library (m_npy npy.f90)

#set(OpenBLAS64_INCLUDE /home/fox/Git/OpenBLAS/include)
#set(OpenBLAS64_LIB /home/fox/Git/OpenBLAS/lib/libopenblas_haswellp-r0.3.21.dev.a)
#${OpenBLAS64_LIB} ${OpenBLAS64_LIBRARIES}
set(SOURCES dsyev_hamiltonian.f90)
add_executable(dsyev_hamiltonian ${SOURCES})
target_link_libraries(dsyev_hamiltonian  OpenBLAS64::OpenBLAS64  Threads::Threads m_npy)
target_include_directories(dsyev_hamiltonian PRIVATE ${Threads_INCLUDE_DIRS} ${OpenBLAS64_INCLUDE_DIRS})
target_compile_options(dsyev_hamiltonian PRIVATE -finteger-4-integer-8)