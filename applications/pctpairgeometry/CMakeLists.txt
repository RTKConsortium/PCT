wrap_ggo(pctpairgeometry_GGO_C pctpairgeometry.ggo)
add_executable(
  pctpairgeometry
  pctpairgeometry.cxx
  ${pctpairgeometry_GGO_C}
)
target_link_libraries(pctpairgeometry PCT)

# Ensure binary is placed in the central PCT binary dir
set_target_properties(
  pctpairgeometry
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_BINARY_DIR}/bin"
)

# Installation code
install(
  TARGETS
    pctpairgeometry
  RUNTIME
    DESTINATION ${PCT_INSTALL_RUNTIME_DIR}
    COMPONENT Runtime
  LIBRARY
    DESTINATION ${PCT_INSTALL_LIB_DIR}
    COMPONENT RuntimeLibraries
  ARCHIVE
    DESTINATION ${PCT_INSTALL_ARCHIVE_DIR}
    COMPONENT Development
)
