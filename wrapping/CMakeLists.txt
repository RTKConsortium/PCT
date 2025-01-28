itk_wrap_module(PCT)
itk_auto_load_submodules()
itk_end_wrap_module()

wrap_itk_python_bindings_install(/itk "RTK"
    ${PCT_SOURCE_DIR}/applications/pctfdktwodweights/pctfdktwodweights.py
)
