cmake_minimum_required(VERSION 3.16.3)
project(PCTExamples)

set(ExampleSpecificComponents
  PCT
  )

if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED COMPONENTS ITKImageIO ITKTransformIO ${ExampleSpecificComponents})
else()
  # When being built as part of ITK, ITKImageIO and ITKTransformIO
  # lists of modules are not yet ready, causing a configure error
  find_package(ITK REQUIRED COMPONENTS ${ExampleSpecificComponents})
endif()
include(${ITK_USE_FILE})


