#!/usr/bin/env python

import os
import itk
from itk import PCT as pct
from itk import RTK as rtk
from opengate.contrib.protonct.protonct import protonct

output_folder = "/tmp/output"
number_of_projections = 720

# # Generate some data
gate_folder = os.path.join(output_folder, "gate")
protonct(gate_folder, projections=number_of_projections, verbose=False)

# TODO example on how to make data noisy

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
)

# TODO cut the pairs (pctpaircuts is not converted to Python yet)

# Bin the pairs into projections
projections_folder = os.path.join(output_folder, "projections")
os.makedirs(projections_folder, exist_ok=True)
for p in range(number_of_projections):
    pct.pctbinning(
        input=os.path.join(pairs_folder, f"pairs{p:04d}.mhd"),
        output=os.path.join(projections_folder, f"projections{p}.mhd"),
        source=-1000.0,
        size=[200, 1, 200],
        spacing=[2.0, 1.0, 1.0],
        verbose=True,
    )

# Build geometry
geometry = os.path.join(output_folder, "geometry.xml")
rtk.rtksimulatedgeometry(
    nproj=number_of_projections, output=geometry, sdd=1000.0 + 110.0, sid=1000.0
)

# Reconstruct
pct.pctfdk(
    geometry=geometry,
    path=projections_folder,
    regexp=r"projections.*\.mhd",
    output=os.path.join(output_folder, "recon.mhd"),
    size=[210, 1, 210],
    verbose=True,
)
