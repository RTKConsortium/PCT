if(PCT_BUILD_SPHINX)
  # Build sphinx documentation in the build directory
  set(PCT_DOC_OUTPUT_DIR ${CMAKE_BINARY_DIR}/sphinx)

  # Copy the source tree documentation to the build tree before calling sphinx.
  # The images are downloaded when conf.py is run by sphinx.
  add_custom_target(copy_sources ALL
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PCT_DOC_OUTPUT_DIR}/applications"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PCT_DOC_OUTPUT_DIR}/documentation"
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${PCT_SOURCE_DIR}/applications" "${PCT_DOC_OUTPUT_DIR}/applications"
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${PCT_SOURCE_DIR}/documentation" "${PCT_DOC_OUTPUT_DIR}/documentation"
    COMMAND ${CMAKE_COMMAND} -E copy "${PCT_SOURCE_DIR}/conf.py" "${PCT_DOC_OUTPUT_DIR}/conf.py"
    COMMAND ${CMAKE_COMMAND} -E copy "${PCT_SOURCE_DIR}/index.md" "${PCT_DOC_OUTPUT_DIR}/index.md"
    COMMENT "Copying documentation sources"
    )

  find_program(SPHINX_EXECUTABLE
    NAMES sphinx-build sphinx-build.exe
    DOC "Sphinx documentation generator"
  )

  add_custom_target( sphinx_doc ALL
    COMMAND ${SPHINX_EXECUTABLE}
    -c ${PCT_DOC_OUTPUT_DIR}                   # conf.py directory
    ${PCT_DOC_OUTPUT_DIR}                      # Source directory containing root index.md
    "${PCT_DOC_OUTPUT_DIR}/_build"             # Output directory
    COMMENT "Generating sphinx documentation"
    DEPENDS copy_sources
    WORKING_DIRECTORY ${PCT_DOC_OUTPUT_DIR}
    )

  set_property(
    DIRECTORY APPEND PROPERTY
    ADDITIONAL_MAKE_CLEAN_FILES
    "${PCT_DOC_OUTPUT_DIR}/_build"
    )
endif()
