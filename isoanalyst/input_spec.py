import pandas as pd
from pathlib import Path
from typing import Union, List, Dict


class InputSpec:
    def __init__(self, df: pd.DataFrame):
        self.df = df

    @classmethod
    def from_csv(cls, fname: Union[str, Path]) -> "InputSpec":
        """Factory method to load InputSpec from CSV"""
        df = pd.read_csv(fname)
        df["filepath"] = df.filepath.apply(Path)
        return cls(df)

    def validate(self):
        """Make sure all specified input files exist and the input meets some simple criteria."""
        df = self.df
        seen_cols = set(df.columns)
        expected_cols = set(
            [
                "filepath",
                "organism",
                "type",
                "isotope",
                "element",
                "condition",
                "replicate",
            ]
        )
        if seen_cols != expected_cols:
            print(f"Missing input filespec cols {expected_cols - seen_cols}")
            return False
        if not df.filepath.is_unique:
            print("Filepath duplication")
            return False
        filepaths = set(f for f in self.df.filepath)
        exists = set(x for x in filepaths if x.exists())
        if not filepaths == exists:
            print(f"Filepath missing {filepaths - exists}")
            return False
        if df[df.condition.isin(self.get_conditions())].isnull().sum().sum() > 0:
            print("Missing values in DF")
            return False
        if "blank" in set(x.lower() for x in df.condition):
            if len(self.get_scan_filepaths("blank")) != 0:
                print("Please do not include blanks in the ALL SCAN data!")
                return False
        return True

    def get_conditions(self) -> List[str]:
        """Get list of conditions as list. Automatically removes "blank" as a condition."""
        remove_blank = lambda lst: [
            i for i in lst if i.lower() != "blank" and not pd.isna(i)
        ]
        return remove_blank(self.df["condition"].unique())

    def get_feature_filepaths(self, condition_name: str) -> List[str]:
        """Get list of feature filepaths for a condition name.

        Args:
            condition_name (str): Name of the condition.
        """
        return self.df[
            (self.df["condition"].str.lower() == condition_name.lower())
            & (self.df["type"] == "f")
        ]["filepath"].values.tolist()

    def get_scan_filepaths(self, condition_name: str) -> List[str]:
        """Get list of all-scan filepaths for condition name.

        Args:
            condition_name (str): Name of the condition.
        """
        return self.df[
            (self.df["condition"].str.lower() == condition_name.lower())
            & (self.df["type"] == "s")
        ]["filepath"].values.tolist()

    def get_filepath_info(self, filepath: str) -> Dict:
        try:
            return dict(self.df[self.df.filepath == filepath].iloc[0])
        except IndexError:
            return {}