wrap_ggo(pctzengbackprojections_GGO_C pctzengbackprojections.ggo)
add_executable(
  pctzengbackprojections
  pctzengbackprojections.cxx
  ${pctzengbackprojections_GGO_C}
)
target_link_libraries(pctzengbackprojections PCT)

# Ensure binary is placed in the central PCT binary dir
set_target_properties(
  pctzengbackprojections
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${PCT_BINARY_DIR}/bin"
)

# Installation code
install(
  TARGETS
    pctzengbackprojections
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
