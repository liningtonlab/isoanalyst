#!/usr/bin/env python3

"""Tool to convert source data in to isoanalyst input format"""

import pandas as pd
import pymzml


def mzml(file_path):
    """Import mzML files derived from applying MSConvert to .raw files."""
    headers = ["ScanIndex", "RetTime", "Mz", "Intensity"]

    # Borrowed and modified from https://github.com/rlinington/ms2analyte/blob/master/ms2analyte/converters/waters.py
    # Waters data includes the lockspray internal calibrant scans as 'MS1' data. These are differentiated from true
    # MS1 data by the 'function' attribute in the spectrum element. Data MS1 scans are function 1. Lockspray scans are
    # assigned the highest possible function number (floating, depends on how many DDA scans were permitted during
    # acquisition setup). Commonly lockspray function=5. This is always 3 for MSe (DIA) data.
    # NOTE: If MS2 functionality is added, this is not an issue, because all MS2 data have ms level = 2, and are
    # therefore all legitimate for inclusion.

    # Parse mzML file and format appropriate scan data as Dataframe
    run = pymzml.run.Reader(file_path)
    input_data = []
    for spec in run:
        # Skip over non-MS1 data
        if spec.ms_level != 1:
            continue
        # Skip lockspray or other functions
        if spec.id_dict.get("function") != 1:
            continue
        scan_number = spec.ID
        retention_time = round(spec.scan_time_in_minutes(), 2)
        for peak in spec.peaks("raw"):
            mz = round(peak[0], 4)
            intensity = int(peak[1])
            input_data.append([scan_number, retention_time, mz, intensity])

        # Print import progress. Useful because importing large mzML files can be slow.
        if spec.index % 100 == 0 and spec.index > 0:
            print("Completed import of scan " + str(spec.index))

    return pd.DataFrame(input_data, columns=headers)


def func001(file_path):
    df = pd.read_csv(
        file_path, usecols=["FunctionScanIndex", "ScanTimeMin", "Mz", "Intensity"]
    )
    return df.rename(
        {"FunctionScanIndex": "ScanIndex", " ScanTimeMin": "RetTime"}, copy=True
    )


def cppis(file_path):
    df = pd.read_csv(
        file_path,
        usecols=[
            "Sample",
            "PrecMz",
            "PrecZ",
            "PrecIntensity",
            "RetTime",
            "ScanRangeLow",
            "ScanRangeHigh",
        ],
    )
    return df.rename({"ScanTimeMin": "RetTime"}, copy=True)