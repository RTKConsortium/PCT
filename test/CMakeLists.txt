itk_module_test()

#-----------------------------------------------------------------------------
# CXX tests

set(
  PCTTests
  pctProtonPairsToDistanceDrivenProjectionTest.cxx
  pctHoleFillingImageFilterTest.cxx
  pctPolynomialMLPFunctionTest.cxx
)

createtestdriver(PCT "${PCT-Test_LIBRARIES}" "${PCTTests}")

itk_add_test(NAME pctProtonPairsToDistanceDrivenProjectionTest
  COMMAND PCTTestDriver pctProtonPairsToDistanceDrivenProjectionTest
)
itk_add_test(NAME pctPolynomialMLPFunctionTest
  COMMAND PCTTestDriver pctPolynomialMLPFunctionTest
)

itk_add_test(NAME pctHoleFillingImageFilterTest COMMAND PCTTestDriver pctHoleFillingImageFilterTest)

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
      -i DATA{Input/pairs0000.mhd,pairs0000.raw} -o ${ITK_TEST_OUTPUT_DIR}/projections.mhd --source 1000 --size=400,20,220
  )
  itk_add_test(NAME pctapppaircutstest
    COMMAND itkTestDriver
      --compare
      DATA{Baseline/baseline_paircuts.mhd,baseline_paircuts.raw}
      ${ITK_TEST_OUTPUT_DIR}/paircuts.mhd
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pctpaircuts
      -i DATA{Input/input_paircuts.mhd,input_paircuts.raw} -o ${ITK_TEST_OUTPUT_DIR}/paircuts.mhd --roiradius 200 --minwepl 70 --maxwepl 180 --fluence .5 --seed 1234
  )
endif()
