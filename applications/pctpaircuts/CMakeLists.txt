WRAP_GGO(pctpaircuts_GGO_C pctpaircuts.ggo)
add_executable(pctpaircuts pctpaircuts.cxx ${pctpaircuts_GGO_C})
target_link_libraries(pctpaircuts PCT)

# Installation code
install(TARGETS pctpaircuts
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
