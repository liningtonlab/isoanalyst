#!/usr/bin/env python
# coding: utf-8
import re
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd

import dereplicator


"""NOTES
See:
https://pandas.pydata.org/pandas-docs/stable/getting_started/basics.html#iteration

1) pandas `.ix` indexer is deprecated and will be removed soon
2) Iterating through DataFrames is "slow" -> vectorize with df.COL.apply
3) If need to iterate through df, use df.itertuples(). This returns a namedtuple
   so you can access attributes using row.COL and the index with row.Index
4) I like to use pathlib Path for OS operations

Extra speedup?
https://pandas.pydata.org/pandas-docs/stable/user_guide/enhancingperf.html#cython-writing-c-extensions-for-pandas
"""


def ppm_tolerance(mass, error=10):
    """Determine high,low error range for a given mass
    Range of (mass-10ppm, mass+10ppm)

    Args:
        mass (float): Mass to calculate error range for

    Returns:
        tuple: Low, High error range
    """
    ppm = mass * (error*1E-6)
    high = mass+ppm
    low = mass-ppm

    return (low, high)


def c_isotope(mass, sign=1):
    """Add mass difference from 12C to 13C to mass

    Args:
        mass (float): Mass to calculate against
        sign (int): Optional: Plus or minus 1. Default=1

    Returns:
        float: Mass plus mass difference from 12C to 13C
    """
    return mass + sign*1.00335


def combine_dfs(dfs):
    """Combine any number of pandas dataframe

    Args:
        dfs (list): List of pd.DataFrame objects

    Returns:
        pd.DataFrame: A combined pandas dataframe
    """
    return pd.concat(dfs, sort=False).reset_index(drop=True)


'''All functions for munging cppis.csv files to create an mz masterlist containing
all real features in unlabelled control samples'''

<<<<<<< HEAD
def drop_cppis(df, filter_rt=0.8, **kwargs):
=======
def drop_cppis(df, inplace=True):
>>>>>>> parent of 4ccd875... Add RT filter to pre-processing
    """ Take cppis csv from MSeXpress and drops unnecessary columns and duplicates

    Args:
        df (pd.DataFrame): cppis type dataframe
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
<<<<<<< HEAD
    # Validate filter
    if not filter_rt >= 0.0:
        raise exc.InvalidFilter()

    # Flexiblility for ignore cols
    ignore_cols = kwargs.get("ignore_cols",
        ['PrecMHplus', 'CPPIS', 'PccChain', 'Mode', 'Func', 'Scan', 'Sequence', 'IsVirtual',
        'FragEff', 'IsoA0Ratio','IsoA1Ratio','IsoA2Ratio','IsoA0RatioCV','IsoA1RatioCV',
        'IsoA2RatioCV','IsoA3RatioCV','IsoA3Ratio','ProdMHplus','ProdMz', 'ProdIntensity',
        'Ar1','Ar3','A0ProdMzBindex'])


    data = df[df['RetTime'] > filter_rt].copy()

    data.drop(ignore_cols, errors='ignore', axis=1, inplace=True)

    return data.drop_duplicates()
=======
    if not inplace:
        data = df.copy()
    else:
        data = df
    if not inplace:
        data = data.drop(['PrecMHplus', 'CPPIS', 'PccChain', 'Mode', 'Func', 'Scan', 'Sequence', 'IsVirtual',
                    'FragEff', 'IsoA0Ratio','IsoA1Ratio','IsoA2Ratio','IsoA0RatioCV','IsoA1RatioCV',
                    'IsoA2RatioCV','IsoA3RatioCV','IsoA3Ratio','ProdMHplus','ProdMz', 'ProdIntensity',
                    'Ar1','Ar3','A0ProdMzBindex'], errors='ignore', axis=1, inplace=inplace)
    else:
        data.drop(['PrecMHplus', 'CPPIS', 'PccChain', 'Mode', 'Func', 'Scan', 'Sequence', 'IsVirtual',
                    'FragEff', 'IsoA0Ratio','IsoA1Ratio','IsoA2Ratio','IsoA0RatioCV','IsoA1RatioCV',
                    'IsoA2RatioCV','IsoA3RatioCV','IsoA3Ratio','ProdMHplus','ProdMz', 'ProdIntensity',
                    'Ar1','Ar3','A0ProdMzBindex'], errors='ignore', axis=1, inplace=inplace)
    # Will return None if inplace=True, else return DataFrame
    return data.drop_duplicates(inplace=inplace)
>>>>>>> parent of 4ccd875... Add RT filter to pre-processing


def getfilenames_cppis(cppis_dir, conditions):
    """Collect files and specify conditions
    This ensures that all the replicates get filed and munged together.
    This is used within the munge_cppis() function
    # TODO: Add robustness
    # To reconsider - overhead of DF may not be worth it

    Args:
        cppis_dir (str or Path): Path to input cppis files

    Returns:
        pd.DataFrame: DataFrame containing files and conditions
    """
    condition_regexp = re.compile(f"({'|'.join(conditions)}|BLANK)")
    find_condition = lambda x: condition_regexp.search(x).group()

    # look in cppis folder for all .csv files
    cppis_dir = Path(cppis_dir)
    cppis_csvs = list(cppis_dir.glob("*.csv"))

    files = pd.DataFrame({
        "path": cppis_csvs,
        "filename": [x.name for x in cppis_csvs],
        "condition": [find_condition(x.name) for x in cppis_csvs]
    })
    return files


def munge_cppis(files, cond, out_dir):
    """Pre-process data for specific condition before detecting isotope labelling

    Args:
        files (pd.DataFrame): DataFrame containing list of filenames from getfilenames_cppis function
        cond (str): condition ie ACE, ACED0, BLANK, etc
        out_dir (str or Path): Directory to put outputs in

    Returns:
        pd.DataFrame: DataFrame which has been pre-processed for IsoTracing
    """
    out_dir = Path(out_dir)
    replicates = files[files['condition'] == cond] # get filenames associated with current condition
    dfs = [pd.read_csv(f) for f in replicates['path']]
    _ = [drop_cppis(df) for df in dfs]
    # All reps dataframe
    # Will be editted inplace along the way
    print("Combining DFs")
    df = combine_dfs(dfs)
    df.to_csv(out_dir.joinpath('All_ions_sorted.csv'), index=False)

    # R-tree based replicate comparison
    print("Performing replicate comparison")
    averaged = dereplicator.replicate_compare(df)
    averaged.to_csv(out_dir.joinpath('All_ions_averaged.csv'), index=False)

    return averaged


def blank_subtract(blank_df, df, inplace=True):
    """Given a DF of blank data and another DF, remove blanks from DF

    Args:
        blank_df (pd.DataFrame): Blanks data frame
        df (pd.DataFrame): Other CPPIS like DF
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
    # Make sure indices are sequential integers 0,1,2,etc...
    df.reset_index(inplace=True, drop=True)
    # Keep a set of IDs to drop
    drop_me = dereplicator.find_overlap(df, blank_df)
    # Will return None if inplace=True, else return DataFrame
    return df.drop(drop_me, inplace=inplace)


def prep_func(fname, **kwargs):
    """Drop unnecessary columns and add metadata to func001-like DF

    Args:
        fname (str or Path): Path to func001-like CSV

    Returns:
        pd.DataFrame: Pre-processed func001-like DF
    """
    # Some defaults with flexibility for kwargs
    sname_char = kwargs.get("sname_char", "_")
    sname_index = kwargs.get("sname_index", 1)
    ignore_cols = kwargs.get("ignore_cols", ['drift','DriftFwhm','QuadMass'])

    fname = Path(fname)
    sname = fname.name.split(sname_char)[sname_index]
    df = pd.read_csv(fname).drop(ignore_cols, axis=1)
    df["Organism"], df["Isotope"], df["Condition"] = split_samplename(sname)
    return df


def split_samplename(s):
    """Take sample name and split into organism, isotope, and condition tuple

    Args:
        s (string): Sample name in format "RLUS1234ISODN

    Returns:
        tuple: organism, isotope, condition strings tuple
    """
    return s[:8], s[8:10], s[10:]


def get_func_slice(df, mz, low_scan, high_scan, seen):
    """Return the indices of func DF given mz, scan range, and seen list/et

    Arguments:
        df (pd.DataFrame): func001-like DataFrame to slice
        mz (float): Mass to query
        low_scan (int): LowScan
        high_scan (int): HighScan
        seen (set): Set of seen indices to ignore

    Returns:
        iterable: List-like of indices in func001-like DF
    """
    # print(len(seen))
    mz_tol = ppm_tolerance(mz)
    masks = (
        df['MZ'] >= mz_tol[0],
        df['MZ'] <= mz_tol[1],
        df['FunctionScanIndex'] >= low_scan,
        df['FunctionScanIndex'] <= high_scan,
        [idx not in seen for idx in df.index],
    )
    func_slice = df[reduce(np.logical_and, masks)]
    seen.update(func_slice.index)
    return func_slice.index


def calc_exp(g):
    """Calculate the avgPrecMZ, Lowscan and HighScan for an ExpId

    Args:
        g (pd.DataFrame): Exp_ID DataFrame

    Returns:
        tuple: (avgPrecMz, LowScan, HighScan)
    """
    return round(g.PrecMz.mean(), 4), g['LowScan'].min(), g['LowScan'].max()


def isotope_slicer(df, mz, low_scan, high_scan, exp_id, seen):
    """Iteratively find all isotope data associated with a given mass.
    Continues until a slice has less than five datapoints.

    Args:
        df (pd.DataFrame): Func001 dataframe
        mz (float): Precursor mass to start scanning from
        low_scan (int): Low scan value in CPPIS
        high_scan (int): High scan value in CPPIS

    Labels DataFrame inplace
    """
    counter = 0
    to_mark = {}
    # Find base ion,
    to_mark[counter] = get_func_slice(df, mz, low_scan, high_scan, seen)

    # find isotopes +/- C13
    # Initialize while loop
    # Need to look forwards only because CPPIS only contains M0 peaks
    this_mz = mz
    while True:
        # print(f"Sign = {sign}")
        counter += 1
        mn = c_isotope(this_mz)
        this_mz = mn
        # print(f"{exp_id} - {this_mz} - {counter}")
        indices = get_func_slice(df, mn, low_scan, high_scan, seen)
        to_mark[counter] = indices
        # print(f"Found {len(indices)} features")
        # Stop condition
        if len(indices) < 5:
            break

    for c, idc in to_mark.items():
        df.loc[idc, "Isotopomer"] = f"M{c}"
        df.loc[idc, "Exp_ID"] = exp_id
