WRAP_GGO(pctbackprojections_GGO_C pctbackprojections.ggo)
add_executable(pctbackprojections pctbackprojections.cxx ${pctbackprojections_GGO_C})
target_link_libraries(pctbackprojections PCT)

# Installation code
install(TARGETS pctbackprojections
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
