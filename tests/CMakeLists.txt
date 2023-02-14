# -- test_numberfield -------------------------------------------------

include_directories(../include ${Boost_INCLUDE_DIR})
set(SOURCES test_boost_string.cpp)
add_executable(test_boost_string ${SOURCES})
add_test(NAME test_boost_string COMMAND test_boost_string)

include_directories(../include ${Boost_INCLUDE_DIR})
set(SOURCES test_triangle_group.cpp ../src/triangle_group.cpp)
add_executable(test_triangle_group ${SOURCES})
add_test(NAME test_triangle_group COMMAND test_triangle_group)

include_directories(../include ${Boost_INCLUDE_DIR})
set(SOURCES test_ring_extension.cpp)
add_executable(test_ring_extension ${SOURCES})
add_test(NAME test_ring_extension COMMAND test_ring_extension)
