WRAP_GGO(pctbinning_GGO_C pctbinning.ggo)
add_executable(pctbinning pctbinning.cxx ${pctbinning_GGO_C})
target_link_libraries(pctbinning PCT)

# Installation code
install(TARGETS pctbinning
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
