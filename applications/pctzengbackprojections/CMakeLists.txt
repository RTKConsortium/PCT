WRAP_GGO(pctzengbackprojections_GGO_C pctzengbackprojections.ggo)
add_executable(pctzengbackprojections pctzengbackprojections.cxx ${pctzengbackprojections_GGO_C})
target_link_libraries(pctzengbackprojections PCT)

# Installation code
install(TARGETS pctzengbackprojections
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
