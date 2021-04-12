import pytest
from pathlib import Path
import pandas as pd

from isoanalyst.input_spec import InputSpec

TESTDIR = Path(__file__).parent
TESTFILE = TESTDIR / "test_input_spec.csv"


def is_lower_str(s):
    assert s.lower() == s


def test_csv_factory():
    inp = InputSpec.from_csv(TESTFILE)
    assert isinstance(inp, InputSpec)
    assert isinstance(inp.df, pd.DataFrame)
    inp.df.applymap(lambda s: is_lower_str(s) if type(s) == str else s)


def test_get_conditions():
    inp = InputSpec.from_csv(TESTFILE)
    expected = ["cond1"]
    assert expected == inp.get_conditions()


def test_get_feature_filepaths_blanks():
    inp = InputSpec.from_csv(TESTFILE)
    expected = ["blanks_rep1.csv", "blanks_rep2.csv"]
    assert expected == inp.get_feature_filepaths("blank")


def test_get_feature_filepaths_condition():
    inp = InputSpec.from_csv(TESTFILE)
    expected = [
        "feature_list_rep1.csv",
        "feature_list_rep2.csv",
        "feature_list_rep3.csv",
    ]
    assert expected == inp.get_feature_filepaths("cond1")