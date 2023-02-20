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


# include_directories(../include ${Boost_INCLUDE_DIR})
# set(SOURCES generate_group_open.cpp)
# add_executable(generate_group_open.x ${SOURCES})

add_library (m_npy npy.f90)

set(SOURCES dsyev_hamiltonian.f90)
add_executable(dsyev_hamiltonian ${SOURCES})
target_link_libraries(dsyev_hamiltonian LAPACK::LAPACK m_npy)
target_include_directories(dsyev_hamiltonian PRIVATE ${LAPACK_INCLUDE_DIRS})
