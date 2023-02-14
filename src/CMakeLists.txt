# include_directories(../include ${Boost_INCLUDE_DIR})
# set(SOURCES main.cpp)
# add_executable(triangle_group.x ${SOURCES})
# target_include_directories(group_algebra.x PUBLIC ../include)

# -- generate_group.cpp

include_directories(../include ${Boost_INCLUDE_DIR})
set(SOURCES generate_group.cpp ../src/generate_group.cpp ../src/triangle_group.cpp)
add_executable(generate_group.x ${SOURCES})