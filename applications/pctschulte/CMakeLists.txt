WRAP_GGO(pctschulte_GGO_C pctschulte.ggo)
add_executable(pctschulte pctschulte.cxx ${pctschulte_GGO_C})
target_link_libraries(pctschulte PCT)

# Installation code
install(TARGETS pctschulte
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
