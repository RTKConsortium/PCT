WRAP_GGO(pctfdk_GGO_C pctfdk.ggo)
add_executable(pctfdk pctfdk.cxx ${pctfdk_GGO_C})
target_link_libraries(pctfdk PCT)

# Installation code
install(TARGETS pctfdk
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
