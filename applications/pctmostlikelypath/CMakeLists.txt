WRAP_GGO(pctmostlikelypath_GGO_C pctmostlikelypath.ggo)
add_executable(pctmostlikelypath pctmostlikelypath.cxx ${pctmostlikelypath_GGO_C})
target_link_libraries(pctmostlikelypath PCT)

# Installation code
install(TARGETS pctmostlikelypath
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
