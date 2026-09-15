#!/usr/bin/env python

import os
import sys
import warnings
import glob
from itk import PCT as pct
from itk import RTK as rtk
from opengate.contrib.protonct.protonct import protonct

if len(sys.argv) < 2:
    print("Usage: python GradientDescentReconstruction.py <outputfolder>")
    sys.exit(1)

output_folder = sys.argv[1]

number_of_projections = 360

# Generate some data
gate_folder = os.path.join(output_folder, "gate")
protonct(
    gate_folder,
    projections=number_of_projections,
    protons_per_projection=1000,
    verbose=False,
)

# Add noise to generated data
for file in ["PhaseSpaceIn", "PhaseSpaceOut"]:
    pct.pctaddnoise(
        input=os.path.join(gate_folder, f"{file}.root"),
        output=os.path.join(gate_folder, f"{file}_noisy.root"),
        tree=file,
        material_budget=0.01,
        tracker_distance=10.0,
        noise_position=1.0,
        noise_energy=1.0,
        seed=1234,
    )

# Convert GATE data to PCT list-mode
pairs_folder = os.path.join(output_folder, "pairs")
os.makedirs(pairs_folder, exist_ok=True)
pct.pctpairprotons(
    input_in=os.path.join(gate_folder, "PhaseSpaceIn.root"),
    input_out=os.path.join(gate_folder, "PhaseSpaceOut.root"),
    output=os.path.join(pairs_folder, "pairs.mhd"),
    psin="PhaseSpaceIn",
    psout="PhaseSpaceOut",
    plane_in=-110.0,
    plane_out=110.0,
    verbose=True,
    store_time=True,
)

# TODO cut the pairs (pctpaircuts is not converted to Python yet)

# Build geometry
geometry = os.path.join(output_folder, "geometry.xml")
rtk.rtksimulatedgeometry(
    nproj=number_of_projections, output=geometry, sdd=1000.0 + 110.0, sid=1000.0
)

# Generate stopping power fit
sp_fit = os.path.join(output_folder, "sp_fit.txt")
pct.pctstoppingpower(output=sp_fit)

# Reconstruct
pct.pctgradientdescent(
    path=pairs_folder,
    regexp="pairs.*\\.mhd",
    geometry=geometry,
    sp_fit=sp_fit,
    output_dir=os.path.join(output_folder, "gradient_descent"),
    size=[110, 3, 110],
    spacing=[2.0, 1.0, 2.0],
    optimizer="Adagrad",
    verbose=True,
    number_of_iterations=10,
)
