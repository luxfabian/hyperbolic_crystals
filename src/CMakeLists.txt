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