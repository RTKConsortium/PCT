import sys
import importlib

itk_module = sys.modules["itk"]
pct_module = getattr(itk_module, "PCT")

# Load the CMake-generated version and assign it to `itk.PCT.__version__`.
pct_version = importlib.import_module("itk.pctConfig").PCT_GLOBAL_VERSION_STRING
setattr(pct_module, "__version__", pct_version)

# Import PCT submodules
pct_submodules = [
    "itk.pctargumentparser",
    "itk.pctExtras",
]
for mod_name in pct_submodules:
    mod = importlib.import_module(mod_name)
    for a in dir(mod):
        if a[0] != "_":
            setattr(pct_module, a, getattr(mod, a))

# Application modules
_app_modules = [
    "pctbinning",
    "pctfdk",
    "pctpairprotons",
    "pctweplfit",
]

# Dynamically access make_application_func from pctExtras
pct_extras = importlib.import_module("itk.pctExtras")
make_application_func = getattr(pct_extras, "make_application_func")

# Dynamically register applications
for app_name in _app_modules:
    setattr(pct_module, app_name, make_application_func(app_name))
