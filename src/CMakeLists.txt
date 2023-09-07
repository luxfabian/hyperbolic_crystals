# -- fuchsian generate group

set(SOURCES fuchsian_generate_group.cpp ../src/fuchsian_group.cpp ../src/fuchsian_ring.cpp ../include/fuchsian_ring.h ../include/fuchsian_group.h )
add_executable(fuchsian_generate_group ${SOURCES})
target_include_directories(fuchsian_generate_group PRIVATE ../include)

#-- fuchsian_generate_reg.cpp

include_directories(../include)
set(SOURCES fuchsian_generate_reg.cpp ../src/fuchsian_group.cpp ../src/fuchsian_ring.cpp)
add_executable(fuchsian_generate_reg ${SOURCES})
target_link_libraries(fuchsian_generate_reg OpenMP::OpenMP_CXX)
target_compile_options(fuchsian_generate_reg PRIVATE -ftrapv)

#-- generate_group.cpp

include_directories(../include)
set(SOURCES generate_group.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp)
add_executable(generate_group ${SOURCES})

#-- generate_reg.cpp

include_directories(../include)
set(SOURCES generate_reg.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp)
add_executable(generate_reg ${SOURCES})
target_link_libraries(generate_reg OpenMP::OpenMP_CXX)
target_compile_options(generate_reg PRIVATE -ftrapv)

#-- tffg_generate_group.cpp

include_directories(../include)
set(SOURCES tffg_generate_group.cpp ../src/tffg_group.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp )
add_executable(tffg_generate_group ${SOURCES})


#-- tffg_generate_reg.cpp

include_directories(../include)
set(SOURCES tffg_generate_reg.cpp ../src/tffg_group.cpp ../src/triangle_group.cpp ../src/ring_extension.cpp )
add_executable(tffg_generate_reg ${SOURCES})
target_link_libraries(tffg_generate_reg OpenMP::OpenMP_CXX)
target_compile_options(tffg_generate_reg PRIVATE -ftrapv)

#-- dsyev_hamiltonian.cpp

if(OpenBLAS64_FOUND)
  add_library (m_npy npy.f90)
  set(SOURCES dsyev_hamiltonian.f90)
  add_executable(dsyev_hamiltonian ${SOURCES})
  target_link_libraries(dsyev_hamiltonian  OpenBLAS64::OpenBLAS64  m_npy)
  target_include_directories(dsyev_hamiltonian PRIVATE  ${OpenBLAS64_INCLUDE_DIRS})
  target_compile_options(dsyev_hamiltonian PRIVATE -finteger-4-integer-8)
endif()
