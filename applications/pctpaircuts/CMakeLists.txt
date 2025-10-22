wrap_ggo(pctpaircuts_GGO_C pctpaircuts.ggo)
add_executable(
  pctpaircuts
  pctpaircuts.cxx
  ${pctpaircuts_GGO_C}
)
target_link_libraries(pctpaircuts PCT)

# Ensure binary is placed in the central PCT binary dir
set_target_properties(
  pctpaircuts
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${PCT_BINARY_DIR}/bin"
)

# Installation code
install(
  TARGETS
    pctpaircuts
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
