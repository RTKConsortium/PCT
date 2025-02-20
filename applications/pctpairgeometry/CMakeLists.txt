WRAP_GGO(pctpairgeometry_GGO_C pctpairgeometry.ggo)
add_executable(pctpairgeometry pctpairgeometry.cxx ${pctpairgeometry_GGO_C})
target_link_libraries(pctpairgeometry PCT)

# Installation code
install(TARGETS pctpairgeometry
  RUNTIME DESTINATION ${PCT_INSTALL_RUNTIME_DIR} COMPONENT Runtime
  LIBRARY DESTINATION ${PCT_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${PCT_INSTALL_ARCHIVE_DIR} COMPONENT Development)
