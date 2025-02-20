WRAP_GGO(pctprojections_GGO_C pctprojections.ggo)
add_executable(pctprojections pctprojections.cxx ${pctprojections_GGO_C})
target_link_libraries(pctprojections PCT)

# Installation code
install(TARGETS pctprojections
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
