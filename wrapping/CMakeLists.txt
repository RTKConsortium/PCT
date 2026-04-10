itk_wrap_module(PCT)
itk_auto_load_submodules()
itk_end_wrap_module()

set(PCT_VERSION_SCRIPT "${PCT_BINARY_DIR}/applications/pctversion.py")
configure_file(
  ${PCT_SOURCE_DIR}/applications/pctversion.py.in
  ${PCT_VERSION_SCRIPT}
)

file(
  GLOB PCT_GROUP_SCRIPTS
  CONFIGURE_DEPENDS
  "${PCT_SOURCE_DIR}/applications/pct*_group.py"
)
file(
  GLOB PCT_APP_SCRIPTS
  CONFIGURE_DEPENDS
  "${PCT_SOURCE_DIR}/applications/pct*/pct*.py"
)

wrap_itk_python_bindings_install(/itk "PCT"
  __init_pct__.py
  ${PCT_GROUP_SCRIPTS}
  ${PCT_APP_SCRIPTS}
  ${PCT_VERSION_SCRIPT}
  ${PCT_SOURCE_DIR}/wrapping/pctExtras.py
  ${PCT_SOURCE_DIR}/applications/pctargumentparser.py
)
