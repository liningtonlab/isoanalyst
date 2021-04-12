import pandas as pd
from pathlib import Path
from typing import Union, List


class InputSpec:
    def __init__(self, df: pd.DataFrame):
        self.df = df

    @classmethod
    def from_csv(cls, fname: Union[str, Path]) -> "InputSpec":
        """Factory method to load InputSpec from CSV"""
        df = pd.read_csv(fname)
        # Cleanup strings
        # make them all lowercase for consitency
        df = df.applymap(lambda s: s.lower() if type(s) == str else s)
        return cls(df)

    def get_conditions(self) -> List[str]:
        """Get list of conditions as list. Automatically removes "blank" as a condition."""
        remove_blank = lambda lst: [i for i in lst if i != "blank" and not pd.isna(i)]
        return remove_blank(self.df["condition"].unique())

    def get_feature_filepaths(self, condition_name: str) -> List[str]:
        """Get list of feature filepaths from a condition name.

        Args:
            condition_name (str): Name of the condition.
        """
        return self.df[
            (self.df["condition"] == condition_name.lower()) & (self.df["type"] == "f")
        ]["filepath"].values.tolist()
