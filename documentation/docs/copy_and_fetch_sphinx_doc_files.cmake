# List of files to be copied to PCT_DOC_OUTPUT_DIR
list(APPEND DOCSRCFILES
      "conf.py"
      "index.md"
    )

# List of subdirectories to be copied to PCT_DOC_OUTPUT_DIR
set(DOCSUBDIRS
      "examples"
      "applications"
      "documentation"
   )

# Append the files in the subdirectories
foreach(subdir ${DOCSUBDIRS})
  file(GLOB_RECURSE SUBSUBDIR LIST_DIRECTORIES false RELATIVE "${PCT_SOURCE_DIR}" "${PCT_SOURCE_DIR}/${subdir}/*")
  list(APPEND DOCSRCFILES ${SUBSUBDIR})
endforeach()

# Make the directory and copy the file if they are different
foreach(file ${DOCSRCFILES})
  get_filename_component(dir "${PCT_DOC_OUTPUT_DIR}/${file}" DIRECTORY)
  file(MAKE_DIRECTORY ${dir})
  file(COPY_FILE "${PCT_SOURCE_DIR}/${file}" "${PCT_DOC_OUTPUT_DIR}/${file}" ONLY_IF_DIFFERENT)
endforeach()

# Fetch external data from https://data.kitware.com
# Used by Sphinx conf.py during setup.
#
if(POLICY CMP0009)
    cmake_policy(SET CMP0009 NEW)
endif()
set(link_content sha512)
file(GLOB_RECURSE content_files
  "${PCT_DOC_OUTPUT_DIR}/*${link_content}")
foreach(content_file ${content_files})
  file(RELATIVE_PATH content_file_rel_path ${PCT_DOC_OUTPUT_DIR} ${content_file})
  get_filename_component(last_ext "${content_file_rel_path}" LAST_EXT)
  string(REPLACE "${last_ext}" "" content_file_rel_path "${content_file_rel_path}")

  if(NOT EXISTS "${PCT_DOC_OUTPUT_DIR}/${content_file_rel_path}" OR
     NOT "${PCT_DOC_OUTPUT_DIR}/${content_file_rel_path}" IS_NEWER_THAN "${content_file}")
    file(READ "${content_file}" hash)
    set(URL "https://data.kitware.com:443/api/v1/file/hashsum/${link_content}/${hash}/download")
    string(REGEX REPLACE "\n" "" URL "${URL}")
    file(DOWNLOAD "${URL}" "${PCT_DOC_OUTPUT_DIR}/${content_file_rel_path}")
  endif()
  list(APPEND DOCSRCFILES "${content_file_rel_path}")
endforeach()

# Clean-up: remove all files which are not in the input files or in the _build and documentation/docs/ExternalData/ subfolders
file(GLOB_RECURSE filesInOutput LIST_DIRECTORIES false RELATIVE "${PCT_DOC_OUTPUT_DIR}" "${PCT_DOC_OUTPUT_DIR}/*")
list(REMOVE_ITEM filesInOutput ${DOCSRCFILES})
file(GLOB_RECURSE filesInOutputBuild LIST_DIRECTORIES false RELATIVE "${PCT_DOC_OUTPUT_DIR}" "${PCT_DOC_OUTPUT_DIR}/_build/*")
list(REMOVE_ITEM filesInOutput ${filesInOutputBuild})
foreach(file ${filesInOutput})
  file(REMOVE "${PCT_DOC_OUTPUT_DIR}/${file}")
endforeach()
