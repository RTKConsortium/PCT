import torch
from itk import PCT as pct
from itk.pctgradientdescent import Velocity, EnergyLoss, RSPInterp, SPWater


def test_gradientdescent_velocity():
    e = torch.rand(10, requires_grad=True, dtype=torch.float64) * 1000
    torch.autograd.gradcheck(Velocity.apply, (e))


def test_gradientdescent_energyloss():
    number_of_protons = 100
    es = torch.rand(number_of_protons, dtype=torch.float64) * 200.0
    rsps = torch.rand(number_of_protons, requires_grad=True, dtype=torch.float64)
    sp_h2o = torch.rand(number_of_protons, dtype=torch.float64) * 10.0
    d = torch.rand(number_of_protons, dtype=torch.float64) * 10.0
    torch.autograd.gradcheck(EnergyLoss.apply, (es, rsps, sp_h2o, d))


def test_gradientdescent_rspinterp():
    number_of_protons = 100
    rsp = torch.rand((10, 10, 10), requires_grad=True, dtype=torch.float64) * 2.0
    rsp_spacing = (10.0, 10.0, 10.0)
    coords = torch.rand((3, number_of_protons), dtype=torch.float64) * 100.0 - 50.0
    torch.autograd.gradcheck(
        RSPInterp.apply, (rsp, rsp_spacing, coords[0], coords[1], coords[2])
    )


def test_gradientdescent_spwater():
    energies = torch.linspace(0.0, 200.0, 200, dtype=torch.float64)
    sps = torch.rand(energies.shape, dtype=torch.float64)
    es = torch.rand((100,), requires_grad=True, dtype=torch.float64) * 200.0
    torch.autograd.gradcheck(SPWater.apply, (energies, sps, es))
