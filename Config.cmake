# Set the name and version of the package
set(PACKAGE_NAME "OpenBLAS64")
set(PACKAGE_VERSION "0.32.1")

# Set the include directory
set(PACKAGE_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")

# Set the library directory and name
set(PACKAGE_LIBRARY_DIRS "${CMAKE_CURRENT_LIST_DIR}/lib")
set(PACKAGE_LIBRARIES "mypackage")

# Create an imported target for the package
add_library(${PACKAGE_LIBRARIES} INTERFACE IMPORTED)
set_target_properties(${PACKAGE_LIBRARIES} PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${PACKAGE_INCLUDE_DIRS}"
  INTERFACE_LINK_DIRECTORIES "${PACKAGE_LIBRARY_DIRS}"
)

# Export the target
export(TARGETS ${PACKAGE_LIBRARIES}
  NAMESPACE ${PACKAGE_NAME}::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Targets.cmake"
)

# Export the package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}ConfigVersion.cmake"
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion
)
configure_package_config_file(
  "${CMAKE_CURRENT_LIST_DIR}/Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake"
  INSTALL_DESTINATION "lib/cmake/${PACKAGE_NAME}"
)
