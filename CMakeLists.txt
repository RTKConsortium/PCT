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

itk_module_examples()

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
