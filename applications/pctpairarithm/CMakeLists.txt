WRAP_GGO(pctpairarithm_GGO_C pctpairarithm.ggo)
add_executable(pctpairarithm pctpairarithm.cxx ${pctpairarithm_GGO_C})
target_link_libraries(pctpairarithm PCT)

# Installation code
install(TARGETS pctpairarithm
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
