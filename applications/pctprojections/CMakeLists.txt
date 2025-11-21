wrap_ggo(pctprojections_GGO_C pctprojections.ggo)
add_executable(
  pctprojections
  pctprojections.cxx
  ${pctprojections_GGO_C}
)
target_link_libraries(pctprojections PCT)

# Ensure binary is placed in the central PCT binary dir
set_target_properties(
  pctprojections
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_BINARY_DIR}/bin"
)

# Installation code
install(
  TARGETS
    pctprojections
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
