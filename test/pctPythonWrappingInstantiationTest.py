import itk
from itk import PCT as pct

imageType = itk.Image[itk.F, 3]
pct.FDKDDWeightProjectionFilter[imageType].New()
pct.FDKDDConeBeamReconstructionFilter[imageType].New()
pct.FDKDDBackProjectionImageFilter[imageType, imageType].New()
pct.ProtonPairsToDistanceDrivenProjection[imageType, imageType].New()
pct.SmallHoleFiller[imageType]()

imageType = itk.Image[itk.F, 4]
pct.DDParkerShortScanImageFilter[imageType, imageType].New()
pct.ProtonPairsToBackProjection[imageType, imageType].New()
pct.ZengWeightedBackProjectionImageFilter[imageType].New()
pct.SmallHoleFiller[imageType]()

pct.ThirdOrderPolynomialMLPFunction[itk.D].New()
pct.ThirdOrderPolynomialMLPFunction[itk.F].New()

pct.SchulteMLPFunction.New()
pct.PolynomialMLPFunction.New()

for t1 in [itk.F, itk.D]:
    for t2 in [itk.F, itk.D]:
        pct.EnergyStragglingFunctor[t1, t2]()
        pct.BetheBlochProtonStoppingPower[t1, t2]()
        pct.IntegratedBetheBlochProtonStoppingPowerInverse[t1, t2](1, 1)
