import pandas as pd
from pathlib import Path
from typing import Union


def feature_list(file_path: Union[str, Path], scan_rate: float) -> pd.DataFrame:
    df = pd.read_csv(
        file_path,
        usecols=[1, 2, 3, 4],
        names=["precmz", "rettime", "RtStart", "RtEnd"],
        skiprows=1,  # Skip MzMine Headers
    )
    print(df.head())


if __name__ == "__main__":
    feature_list(
        r"Z:\Linington\working\isoanalyst_example\generalized\feature_lists\20180409_BLANK_M24.csv"
    )
