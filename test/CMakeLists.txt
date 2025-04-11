itk_module_test()

#-----------------------------------------------------------------------------
# CXX tests

set(PCTTests
  pctProtonPairsToDistanceDrivenProjectionTest.cxx
  )

CreateTestDriver(PCT "${PCT-Test_LIBRARIES}" "${PCTTests}")

itk_add_test(NAME pctProtonPairsToDistanceDrivenProjectionTest
  COMMAND PCTTestDriver pctProtonPairsToDistanceDrivenProjectionTest
  )

#-----------------------------------------------------------------------------
# Python tests
if(ITK_WRAP_PYTHON)
  itk_python_add_test(NAME pctPythonWrappingInstantiationTest COMMAND pctPythonWrappingInstantiationTest.py)
  itk_python_add_test(NAME pctOutputArgumentWrappingTest COMMAND pctOutputArgumentWrapping.py)
  itk_python_add_test(NAME pctPairProtonsTest COMMAND pctPairProtonsTest.py)
endif()

#-----------------------------------------------------------------------------
# Applications testing
# This must come before the call to ExternalData_add_target(ITKData)
# in ITKModuleExternal.cmake for DATA{} to be referenced, so we put it here
# rather than in the applications subdirectory.
if(PCT_BUILD_APPLICATIONS)

  itk_add_test(NAME pctappbinningtest
    COMMAND itkTestDriver
      --compare
      DATA{Baseline/projections0000.mhd,projections0000.raw}
      ${ITK_TEST_OUTPUT_DIR}/projections.mhd
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pctbinning
      -i DATA{Input/pairs0000.mhd,pairs0000.raw} -o ${ITK_TEST_OUTPUT_DIR}/projections.mhd --source 1000 --dimension=400,20,220
  )

endif()
