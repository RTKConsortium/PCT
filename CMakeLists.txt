cmake_minimum_required(VERSION 3.16.3)
project(PCT)

set(PCT_LIBRARIES PCT)

if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if(NOT PCT_INSTALL_RUNTIME_DIR)
  set(PCT_INSTALL_RUNTIME_DIR bin)
endif()
if(NOT PCT_INSTALL_LIB_DIR)
  set(PCT_INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR})
endif()
if(NOT PCT_INSTALL_ARCHIVE_DIR)
  set(PCT_INSTALL_ARCHIVE_DIR ${CMAKE_INSTALL_LIBDIR})
endif()

# --------------------------------------------------------
# Build applications
option(PCT_BUILD_APPLICATIONS "Build PCT applications" ${RTK_BUILD_APPLICATIONS})
if(PCT_BUILD_APPLICATIONS)
  add_subdirectory(applications)
endif()

#=========================================================
# Build doc
option(PCT_BUILD_SPHINX "Build Sphinx Documentation" OFF)
if(PCT_BUILD_SPHINX)
  add_subdirectory(documentation/docs)
endif()

#---------------------------------------------------------
# Git hooks
if(EXISTS "${PCT_SOURCE_DIR}/.git/config" AND
    NOT EXISTS "${PCT_SOURCE_DIR}/.git/hooks/pre-commit")
  # Silently ignore the error if the hooks directory is read-only.
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy ${PCT_SOURCE_DIR}/cmake/Hooks/pre-commit
                                     ${PCT_SOURCE_DIR}/.git/hooks/pre-commit
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE  _output
    RESULT_VARIABLE _result
    )
  if(_result AND NOT "${_output}" MATCHES "Error copying file")
    message("${_output}")
  endif()
endif()