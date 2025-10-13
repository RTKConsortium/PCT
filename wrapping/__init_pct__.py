import sys
import importlib

itk_module = sys.modules["itk"]
pct_module = getattr(itk_module, "PCT")

# Import PCT submodules
pct_submodules = [
    "itk.pctversion",
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
    "pctfdktwodweights",
    "pctpairprotons",
    "pctweplfit",
]

# Dynamically access make_application_func from pctExtras
pct_extras = importlib.import_module("itk.pctExtras")
make_application_func = getattr(pct_extras, "make_application_func")

# Dynamically register applications
for app_name in _app_modules:
    setattr(pct_module, app_name, make_application_func(app_name))
