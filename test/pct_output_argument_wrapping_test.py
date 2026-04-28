# Test for MostLikelyPathFunction::Evaluate

from itk import PCT as pct


def test_output_argument_wrapping():
    pos_in = [25.0, -25.0, -110.0]
    pos_out = [-25.0, 25.0, 110.0]
    dir_in = [0.0, 0.0, 1.0]
    dir_out = [0.0, 0.0, 1.0]

    mlp = pct.PolynomialMLPFunction.New()
    mlp.SetPolynomialDegree(2)
    mlp.Init(pos_in, pos_out, dir_in, dir_out)

    result = mlp.Evaluate([-50.0, 0.0, 50.0])
    print(result)
