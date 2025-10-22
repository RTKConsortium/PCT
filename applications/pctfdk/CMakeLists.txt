wrap_ggo(pctfdk_GGO_C pctfdk.ggo)
add_executable(
  pctfdk
  pctfdk.cxx
  ${pctfdk_GGO_C}
)
target_link_libraries(pctfdk PCT)

# Ensure binary is placed in the central PCT binary dir
set_target_properties(
  pctfdk
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${PCT_BINARY_DIR}/bin"
)

# Installation code
install(
  TARGETS
    pctfdk
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
