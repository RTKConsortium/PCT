import pytest
import itk
from itk import PCT as pct
from itk import RTK as rtk


@pytest.fixture()
def parser():
    p = pct.PCTArgumentParser()
    p.add_argument("--string-single", type=str, required=True)
    p.add_argument("--number-single", type=int, required=True)
    p.add_argument("--string-many", type=str, nargs="+", required=True)
    p.add_argument("--number-many", type=int, nargs="+", required=True)
    return p


def assert_args(args):
    assert args.string_single == "a"
    assert args.number_single == 1
    assert args.string_many == ["a", "b", "c"]
    assert args.number_many == [1, 2, 3]


def test_parser_inherits_from_rtk():
    assert issubclass(pct.PCTArgumentParser, rtk.RTKArgumentParser)


def test_version_flag_uses_pct_version(parser):
    assert parser._version == pct.__version__


def test_verbose_flag_inherited(parser):
    args = parser.parse_args(
        [
            "--string-single",
            "a",
            "--number-single",
            "1",
            "--string-many",
            "a,b,c",
            "--number-many",
            "1,2,3",
            "-v",
        ]
    )
    assert args.verbose is True
    assert_args(args)


def test_pctargumentparser_cli(parser):
    args = parser.parse_args(
        [
            "--string-single",
            "a",
            "--number-single",
            "1",
            "--string-many",
            "a,b,c",
            "--number-many",
            "1,2,3",
        ]
    )
    assert_args(args)


def test_pctargumentparser_cli_commas(parser):
    args = parser.parse_args(
        [
            "--string-single",
            "a",
            "--number-single",
            "1",
            "--string-many",
            "a",
            "b",
            "c",
            "--number-many",
            "1",
            "2",
            "3",
        ]
    )
    assert_args(args)


def test_pctargumentparser_python(parser):
    args = parser.parse_kwargs(
        string_single="a",
        number_single=1,
        string_many=["a", "b", "c"],
        number_many=[1, 2, 3],
    )
    assert_args(args)


def test_pctargumentparser_python_commas(parser):
    args = parser.parse_kwargs(
        string_single="a", number_single=1, string_many="a,b,c", number_many="1,2,3"
    )
    assert_args(args)
