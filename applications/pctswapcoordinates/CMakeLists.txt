WRAP_GGO(pctswapcoordinates_GGO_C pctswapcoordinates.ggo)
add_executable(pctswapcoordinates pctswapcoordinates.cxx ${pctswapcoordinates_GGO_C})
target_link_libraries(pctswapcoordinates PCT)

# Installation code
install(TARGETS pctswapcoordinates
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
