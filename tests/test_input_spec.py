from pathlib import Path
import pandas as pd

from isoanalyst.input_spec import InputSpec

TESTDIR = Path(__file__).parent
TESTFILE = TESTDIR / "test_input_spec.csv"


def test_csv_factory():
    inp = InputSpec.from_csv(TESTFILE)
    assert isinstance(inp, InputSpec)
    assert isinstance(inp.df, pd.DataFrame)


def test_get_conditions():
    inp = InputSpec.from_csv(TESTFILE)
    expected = ["COND1"]
    assert expected == inp.get_conditions()


def test_get_feature_filepaths_blanks():
    inp = InputSpec.from_csv(TESTFILE)
    expected = list(
        map(
            Path,
            ["tests/test_files/blanks_rep1.csv", "tests/test_files/blanks_rep2.csv"],
        )
    )
    assert expected == inp.get_feature_filepaths("blank")


def test_get_feature_filepaths_condition():
    inp = InputSpec.from_csv(TESTFILE)
    expected = list(
        map(
            Path,
            [
                "tests/test_files/feature_list_rep1.csv",
                "tests/test_files/feature_list_rep2.csv",
                "tests/test_files/feature_list_rep3.csv",
            ],
        )
    )
    assert expected == inp.get_feature_filepaths("cond1")


def test_get_scan_filepaths():
    inp = InputSpec.from_csv(TESTFILE)
    expected = list(
        map(
            Path,
            [
                "tests/test_files/full_scan_nat_rep1.mzml",
                "tests/test_files/full_scan_nat_rep2.mzml",
                "tests/test_files/full_scan_nat_rep3.mzml",
            ],
        )
    )
    assert expected == inp.get_scan_filepaths("12cond1")


def test_validate_success():
    inp = InputSpec.from_csv(TESTFILE)
    print(inp.df.head())
    assert inp.validate()


def test_validate_fails_missing():
    df = pd.DataFrame(
        [
            {
                "filepath": Path("missing"),
                "organism": "org1",
                "type": "f",
                "element": "c",
                "isotope": 13,
                "condition": "cond1",
                "replicate": 1,
            }
        ]
    )
    inp = InputSpec(df)
    assert inp.validate() == False


def test_validate_fails_duplicated():
    df = pd.DataFrame(
        [
            {
                "filepath": Path("tests/test_files/full_scan_nat_rep1.mzml"),
                "organism": "org1",
                "type": "f",
                "element": "c",
                "isotope": 13,
                "condition": "cond1",
                "replicate": 1,
            },
            {
                "filepath": Path("tests/test_files/full_scan_nat_rep1.mzml"),
                "organism": "org1",
                "type": "f",
                "element": "c",
                "isotope": 13,
                "condition": "cond1",
                "replicate": 2,
            },
        ]
    )
    inp = InputSpec(df)
    assert inp.validate() == False


def test_validate_fails_missing():
    df = pd.DataFrame(
        [
            {
                "filepath": Path("tests/test_files/full_scan_nat_rep1.mzml"),
            }
        ]
    )
    inp = InputSpec(df)
    assert inp.validate() == False


def test_validate_fails_all_scan_blank():
    df = pd.DataFrame(
        [
            {
                "filepath": Path("tests/test_files/full_scan_nat_rep1.mzml"),
                "organism": "org1",
                "type": "f",
                "element": "c",
                "isotope": 13,
                "condition": "blank",
                "replicate": 1,
            },
            {
                "filepath": Path("tests/test_files/full_scan_nat_rep2.mzml"),
                "organism": "org1",
                "type": "s",
                "element": "c",
                "isotope": 13,
                "condition": "blank",
                "replicate": 1,
            },
        ]
    )
    inp = InputSpec(df)
    assert inp.validate() == False


def test_get_filepath_info():
    inp = InputSpec.from_csv(TESTFILE)
    data = inp.get_filepath_info(Path("tests/test_files/full_scan_labelled_rep3.mzml"))
    assert data["organism"] == "ORG1"
    assert data["type"] == "s"
    assert data["isotope"] == 13
    assert data["element"] == "C"
    assert data["condition"] == "COND1"


def test_get_filepath_info_missing():
    inp = InputSpec.from_csv(TESTFILE)
    data = inp.get_filepath_info("MISSING")
    assert len(data) == 0
