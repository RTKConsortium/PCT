WRAP_GGO(pctparkershortscanweighting_GGO_C pctparkershortscanweighting.ggo)
add_executable(pctparkershortscanweighting pctparkershortscanweighting.cxx ${pctparkershortscanweighting_GGO_C})
target_link_libraries(pctparkershortscanweighting PCT)

# Installation code
install(TARGETS pctparkershortscanweighting
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
