find_package(Gengetopt REQUIRED)

#-----------------------------------------------------------------------------
# Find ITK.
# Required to include ITK_USE_FILE in order to Register IO factories
# Force requested modules to be PCT dependencies only, otherwise all
# available factories will try to register themselves.
if (NOT ITK_DIR)
  set(ITK_DIR ${ITK_BINARY_DIR}/CMakeTmp)
endif()
find_package(ITK REQUIRED COMPONENTS ${ITK_MODULE_PCT_DEPENDS} ${ITK_MODULE_PCT-Test_DEPENDS})
include(${ITK_USE_FILE})

add_subdirectory(pctpaircuts)
add_subdirectory(pctpairgeometry)
add_subdirectory(pctmostlikelypath)
add_subdirectory(pctbinning)
add_subdirectory(pctfillholl)
add_subdirectory(pctbackprojectionbinning)
add_subdirectory(pctzengbackprojections)
add_subdirectory(pctparkershortscanweighting)
add_subdirectory(pctprojections)
add_subdirectory(pctbackprojections)
add_subdirectory(pctfdktwodweights)
add_subdirectory(pctfdk)
add_subdirectory(pctschulte)
add_subdirectory(pctswapcoordinates)
add_subdirectory(pctpairarithm)
