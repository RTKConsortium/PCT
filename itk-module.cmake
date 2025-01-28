set(DOCUMENTATION "")

#=========================================================
# Module PCT
#=========================================================
itk_module(PCT
  ENABLE_SHARED
  EXCLUDE_FROM_DEFAULT
  DEPENDS
    ITKIOMeta
    RTK
  TEST_DEPENDS
    ITKTestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
)
