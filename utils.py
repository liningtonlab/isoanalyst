#!/usr/bin/env python
# coding: utf-8
import glob
import os
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


def checkdir(my_dir):
    """Create a directory if is doesn't exist
    Useful for creating output directories.
    
    !! May be unnecessary
    Args:
        my_dir (str or Path): Path to directory
    """
    # Make sure directory exists
    my_dir = Path(my_dir)
    if not my_dir.exists():
        my_dir.mkdir()
        

def ppm_tolerance(mass):
    """Determine high,low error range for a given mass
    Range of (mass-10ppm, mass+10ppm)
    
    Args:
        mass (float): Mass to calculate error range for
    
    Returns:
        tuple: Low, High error range
    """
    ppm = (mass)/((1/10)*1000000) # Equiv 10 PPM
    high = mass+ppm
    low = mass-ppm

    return (low, high)

    
def rt_tolerance(rt):
    """Determine high,low error range for a given retention time
    
    Args:
        rt (float): Retention time to calculate error range for
    
    Returns:
        tuple: Low, High error range
    """
    high = rt + 0.03
    low = rt - 0.03
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


def outerjoin_dfs(dfs):
    """DEPRECATE - this solution works but is very CPU and time intensive
    Much better to use built in pd.concat

    Combine any number of pandas dataframe by outer join
    Useful for combining replicates or experimental conditions
    
    Args:
        dfs (list): List of pd.DataFrame objects
    
    Returns:
        pd.DataFrame: A combined pandas dataframe
    """
    return reduce(lambda left,right: pd.merge(left,right,how='outer'), dfs)


'''All functions for munging cppis.csv files to create an mz masterlist containing all real features in 
unlabelled control samples'''

def drop_cppis(df, inplace=True):
    """ Take cppis csv from MSeXpress and drops unnecessary columns and duplicates
    
    Args:
        df (pd.DataFrame): cppis type dataframe
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
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


def averages(df, col="reps"):
    grouped = df.groupby(col)
    new_data =[]
    for idx, grp in grouped:
        new_data.append(
        {
            "Sample": re.split('[_-]', grp.Sample.iloc[0])[1],
            "PrecMz": round(grp.PrecMz.mean(), 4),
            "RetTime": round(grp.RetTime.mean(), 3),
            "PrecZ": grp.PrecZ.iloc[0], 
            "LowScan": grp.ScanLowRange.min(),
            "HighScan": grp.ScanHighRange.max(),
            "reps": idx
        })
    return pd.DataFrame(new_data)

def cat_averages(df, gcol="reps"):
    """
    index must be C_#
    calculates average mz & RT for each C_# and adds this data to new cols
    also grabs low and high scans for the range to use in all replicates 
    """
    groups = df.set_index(gcol)
    for i in groups.index:
        temp = groups.loc[i,:]
        reps_mz = [temp['PrecMz']]
        reps_RT = [temp['RetTime']]
        groups.loc[i,'PrecMz'] = round(np.average(reps_mz), 4)
        groups.loc[i,'RetTime'] = round(np.average(reps_RT), 3)
        groups.loc[i,'LowScan'] = temp['ScanLowRange'].min()
        groups.loc[i,'HighScan'] = temp['ScanHighRange'].max()
    
    final = groups.reset_index()
    return final

# takes output from averages(df) and munges the table to be used for 
# next steps in processing- combine for label condition & combine for all labels. 
def collapse_avg(df):
    cond = df.loc[0,'Sample']
    name = re.split('[_ -]',cond)
    drop_reps = df.drop_duplicates('reps')
    drop_reps['Sample'] = name[1]
    final = drop_reps.reset_index(drop=True)
    return final


def group_cppis(df, c, new):
    """Apply grouping to cppis-like dataframe inplace
    
    This function takes a df in 'cppis-like' format (uses , loops through each row(ions) and 
    identifies ions that are the same based on mz and RT tolerances groups them together in 
    new column 'new' as c_1,c_2,c_3,etc. new is the new column name, c is the group # assignment typically
    based on experimental condition.

    Args:
        df (pd.DataFrame): cppis-like dataframe to process
        c (str): group # assignment typically based on experimental condition.
        new (str): Name for grouping column
    """
    # create column for "compound" c_#
    df[new] = np.nan
    counter = 1
    # loop through data one row at a time
    seen = set()
    for row in df.itertuples():
        # if in seen, continue
        if row.Index in seen:
            continue

        mz_range = ppm_tolerance(row.PrecMz)
        rt_range = rt_tolerance(row.RetTime)

        # Create masks for selecting data
        low_mz = df['PrecMz'] > mz_range[0]
        high_mz = df['PrecMz'] < mz_range[1]
        low_rt = df['RetTime'] > rt_range[0]
        high_rt = df['RetTime'] < rt_range[1]
        matching_z = df['PrecZ'] == row.PrecZ
        # removes indices which have been labelled from slics
        not_in_seen = df.apply(lambda x: x.name not in seen, axis=1)
        
        # find all other rows with compatible mz, rt, and Z, and assign a c_#
        # Assign a value to the pandas slice
        # df.loc[low_rt&high_rt&low_mz&high_mz&matching_z&not_in_seen, new] = f"{c}_{counter}"
        df.loc[low_rt&high_rt&low_mz&high_mz&matching_z&not_in_seen, new] = f"{c}_{counter}"
        # any_seen = any(row.Index in seen for row in df.loc[low_rt&high_rt&low_mz&high_mz&matching_z].itertuples())
        # if any_seen:
            # print("Some have been seen!")
        seen.update(df.loc[low_rt&high_rt&low_mz&high_mz&matching_z&not_in_seen].index)
        counter += 1

    # All the above edits df inplace
    # return df


# removes ions that are in less than 3 samples based on c_# in column t 
# c_# is assigned in the group_cppis() function
def replicate_filter(df, col, min_reps=3, inplace=True):
    """Remove fragments where there are less than min_reps values for each
    identifer=col
    
    Args:
        df (pd.DataFrame): DataFrame to work on
        col (str): Column containing identifier info
        min_reps (int, optional): Min number of replicated to consider. Defaults to 3.
    """
    grouped = df.groupby(col)
    to_drop = []
    for _, g in grouped:
        if g.shape[0] < min_reps:
            to_drop.extend(g.index)
    # Will return None if inplace=True, else return DataFrame
    return df.drop(to_drop, inplace=inplace)


def OLD_munge_cppis(files, cond, out_dir):
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
    # Will be editted inplace along the wa=y 
    df = combine_dfs(dfs)
    df.to_csv(out_dir.joinpath('All_ions_sorted.csv'), index=False)
    
    # Could be replaces with Rtree/Conn-comp like approach
    # which would likely be faster but less transparent
    group_cppis(df, cond, 'reps')
    df.to_csv(out_dir.joinpath('All_ions_grouped.csv'), index=False)

    replicate_filter(df, 'reps')
    df.to_csv(out_dir.joinpath('All_ions_filtered.csv'), index=False)

    # Reduce to a single function which averages and collapses dataframe
    # averages does return a new dataframe
    averaged = averages(df)
    averaged.to_csv(out_dir.joinpath('All_ions_averaged.csv'), index=False)

    return averaged


# def munge_cppis(files, cond, out_dir, **kwargs):
def munge_cppis(files, cond, out_dir):
    """Pre-process data for specific condition before detecting isotope labelling
    
    Args:
        files (pd.DataFrame): DataFrame containing list of filenames from getfilenames_cppis function
        cond (str): condition ie ACE, ACED0, BLANK, etc
        out_dir (str or Path): Directory to put outputs in
    
    REMOVED:
        cppis_dir (str or Path): Input directory with CPPIS files (now built into getfilenames)
    
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


def OLD_blank_subtract(blank_df, df, inplace=True):
    """Given a DF of blank data and another DF, remove blanks from DF
    
    Args:
        blank_df (pd.DataFrame): Blanks data frame
        df (pd.DataFrame): Other CPPIS like DF
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
    # Keep a set of IDs to drop 
    drop_me = set()
    
    for data in blank_df.itertuples():
        mz_range = ppm_tolerance(data.PrecMz)
        rt_range = rt_tolerance(data.RetTime)
        low_mz = df['PrecMz'] > mz_range[0]
        high_mz = df['PrecMz'] < mz_range[1]
        low_rt = df['RetTime'] > rt_range[0]
        high_rt = df['RetTime'] < rt_range[1]
        matching_z = df['PrecZ'] == data.PrecZ
        
        # This will add unique indices to set for dropping
        drop_me.update(df[low_rt&high_rt&low_mz&high_mz&matching_z].index)

    # Will return None if inplace=True, else return DataFrame
    return df.drop(list(drop_me), inplace=inplace)

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
    # Some defaults with flexibility for kwargs
    sname_char = kwargs.get("sname_char", "_")
    sname_index = kwargs.get("sname_index", 1)
    ignore_cols = kwargs.get("ignore_cols", ['drift','DriftFwhm','QuadMass'])

    fname = Path(fname)
    sname = fname.name.split("_")[sname_index]
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


def mark_func_slice(df, mz, low_scan, high_scan, exp_id, seen, counter=0):
    mz_tol = ppm_tolerance(mz)
    masks = (
        df['MZ'] >= mz_tol[0],
        df['MZ'] <= mz_tol[1],
        df['FunctionScanIndex'] >= low_scan,
        df['FunctionScanIndex'] <= high_scan,
        [idx not in seen for idx in df.index],
    )
    func_slice = df[reduce(np.logical_and, masks)]
    df.loc[func_slice.index, "Isotopomer"] = f"M{counter}"
    df.loc[func_slice.index, "Exp_ID"] = exp_id
    seen.update(func_slice.index)
    return len(func_slice.index) >= 5


def isotope_slicer(df, mz, low_scan, high_scan, exp_id, seen):
    """Recursively find all isotope data associated with a given mass.
    Continues until a slice has less than five datapoints.
    
    Args:
        df (pd.DataFrame): Func001 dataframe
        mz (float): Precursor mass to start scanning from
        low_scan (int): Low scan value in CPPIS
        high_scan (int): High scan value in CPPIS
    
    Labels DataFrame inplace
    """
    df['Isotopomer'] = None
    df['Exp_ID'] = None
    counter = 0
    # Find base ion, 
    mark_func_slice(df, mz, low_scan, high_scan, exp_id, seen, counter)
    counter +=1
    def iso_recur(df, mz, low_scan, high_scan, sign=1):
        # find isotopes +/- C13
        mn = c_isotope(mz, sign=sign)
        should_continue = mark_func_slice(df, mz, low_scan, high_scan, exp_id, seen, counter)
        if should_continue:
            iso_recur(df, mn, low_scan, high_scan)
    # Find all isotopes +1 C13
    iso_recur(df, mz, low_scan, high_scan, sign=1)
    # Find all isotopes -1 C13
    iso_recur(df, mz, low_scan, high_scan, sign=-1)