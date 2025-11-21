wrap_ggo(pctbackprojectionbinning_GGO_C pctbackprojectionbinning.ggo)
add_executable(
  pctbackprojectionbinning
  pctbackprojectionbinning.cxx
  ${pctbackprojectionbinning_GGO_C}
)
target_link_libraries(pctbackprojectionbinning PCT)

# Ensure binary is placed in the central PCT binary dir
set_target_properties(
  pctbackprojectionbinning
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
      "${CMAKE_BINARY_DIR}/bin"
)

# Installation code
install(
  TARGETS
    pctbackprojectionbinning
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
