WRAP_GGO(pctbackprojectionbinning_GGO_C pctbackprojectionbinning.ggo)
add_executable(pctbackprojectionbinning pctbackprojectionbinning.cxx ${pctbackprojectionbinning_GGO_C})
target_link_libraries(pctbackprojectionbinning PCT)

# Installation code
install(TARGETS pctbackprojectionbinning
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
