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
