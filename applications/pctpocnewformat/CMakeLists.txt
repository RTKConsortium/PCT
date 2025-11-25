wrap_ggo(pctpocnewformat_GGO_C pctpocnewformat.ggo)
add_executable(
  pctpocnewformat
  pctpocnewformat.cxx
  ${pctpocnewformat_GGO_C}
)
target_link_libraries(pctpocnewformat PCT)

# Ensure binary is placed in the central PCT binary dir
set_target_properties(
  pctpocnewformat
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${PCT_BINARY_DIR}/bin"
)

# Installation code
install(
  TARGETS
    pctpocnewformat
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
