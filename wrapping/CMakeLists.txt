itk_wrap_module(PCT)
itk_auto_load_submodules()
itk_end_wrap_module()

set(PCT_VERSION_SCRIPT "${PCT_BINARY_DIR}/applications/pctversion.py")
configure_file(
  ${PCT_SOURCE_DIR}/applications/pctversion.py.in
  ${PCT_VERSION_SCRIPT}
)

file(
  GLOB PCT_PYTHON_APP
  CONFIGURE_DEPENDS
  "${PCT_SOURCE_DIR}/applications/pct*/pct*.py"
)

wrap_itk_python_bindings_install(/itk "PCT"
  __init_pct__.py
  ${PCT_PYTHON_APP}
  ${PCT_VERSION_SCRIPT}
  ${PCT_SOURCE_DIR}/wrapping/pctExtras.py
  ${PCT_SOURCE_DIR}/applications/pctargumentparser.py
)

# Copy python applications to the ITK wrapping directory to ensure they can be imported in tests.
# This directory is added to the PYTHONPATH by itk_python_add_test.
if(ITK_DIR)
  set(itk_wrap_python_binary_dir "${ITK_DIR}/Wrapping/Generators/Python")
else()
  set(itk_wrap_python_binary_dir "${ITK_BINARY_DIR}/Wrapping/Generators/Python")
endif()
file(COPY ${PCT_PYTHON_APP} DESTINATION "${itk_wrap_python_binary_dir}")
