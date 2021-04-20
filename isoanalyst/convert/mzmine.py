import pandas as pd
from pathlib import Path
from typing import Union


def feature_list(file_path: Union[str, Path]) -> pd.DataFrame:
    ### ADD A SCAN WINDOW AROUND CENTER SCAN DURING SCRAPE
    df = pd.read_csv(file_path, usecols=["row m/z", "row retention time"])
    df.rename(
        columns={"row m/z": "precmz", "row retention time": "rettime"}, inplace=True
    )
    df["sample"] = file_path.name
    # Have to assume that MzMine features are all Z=1
    # because it does not properly return charge data...
    df["precz"] = 1
    return df
