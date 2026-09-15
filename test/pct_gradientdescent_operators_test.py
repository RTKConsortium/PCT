import torch
from itk import PCT as pct
from itk.pctgradientdescent import Velocity, EnergyLoss, RSPInterp, SPWater


def test_gradientdescent_velocity():
    torch.set_default_dtype(torch.float64)

    e = torch.rand(10, requires_grad=True) * 1000
    torch.autograd.gradcheck(Velocity.apply, (e))


def test_gradientdescent_energyloss():
    torch.set_default_dtype(torch.float64)

    number_of_protons = 100
    es = torch.rand(number_of_protons) * 200
    rsps = torch.rand(number_of_protons, requires_grad=True)
    sp_h2o = torch.rand(number_of_protons) * 10
    d = torch.rand(number_of_protons) * 10
    torch.autograd.gradcheck(EnergyLoss.apply, (es, rsps, sp_h2o, d))


def test_gradientdescent_rspinterp():
    torch.set_default_dtype(torch.float64)

    number_of_protons = 100
    rsp = torch.rand((10, 10, 10), requires_grad=True) * 2
    rsp_spacing = (10, 10, 10)
    ps = torch.rand((number_of_protons, 3)) + -5 * 10
    torch.autograd.gradcheck(RSPInterp.apply, (rsp, rsp_spacing, ps[0], ps[1], ps[2]))


def test_gradientdescent_spwater():
    torch.set_default_dtype(torch.float64)

    energies = torch.linspace(0.0, 200.0, 200)
    sps = torch.rand(energies.shape)
    es = torch.rand((100,), requires_grad=True) * 200.0
    torch.autograd.gradcheck(SPWater.apply, (energies, sps, es))
