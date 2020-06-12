#!/usr/bin/env python
# coding: utf-8
import re
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress, ttest_ind

import isotracer.dereplicator as dereplicator
import isotracer.exceptions as exc


def ppm_tolerance(mass, error=10):
    """Determine high,low error range for a given mass
    Range of (mass-10ppm, mass+10ppm)

    Args:
        mass (float): Mass to calculate error range for

    Returns:
        tuple: Low, High error range
    """
    ppm = mass * (error * 1e-6)
    high = mass + ppm
    low = mass - ppm

    return (low, high)


def c_isotope(mass, sign=1):
    """Add mass difference from 12C to 13C to mass

    Args:
        mass (float): Mass to calculate against
        sign (int): Optional: Plus or minus 1. Default=1

    Returns:
        float: Mass plus mass difference from 12C to 13C
    """
    return mass + sign * 1.00335


def n_isotope(mass, sign=1):
    """Add mass difference from 14N to 15N to mass
    14N = 14.003074
    15N = 15.000109

    Args:
        mass (float): Mass to calculate against
        sign (int): Optional: Plus or minus 1. Default=1

    Returns:
        float: Mass plus mass difference from 14N to 15N
    """
    return mass + sign * 0.9970349


def combine_dfs(dfs):
    """Combine any number of pandas dataframe

    Args:
        dfs (list): List of pd.DataFrame objects

    Returns:
        pd.DataFrame: A combined pandas dataframe
    """
    return pd.concat(dfs, sort=False).reset_index(drop=True)


"""All functions for munging cppis.csv files to create an mz masterlist containing
all real features in unlabelled control samples"""


def prep_cppis(fname, min_rt=0.8, **kwargs):
    """ Take cppis csv from MSeXpress and drops unnecessary columns and duplicates

    Args:
        fname (str or Path): Path to cppis type dataframe
        min_rt(float, optional): Filter RT less than this value. MUST BE >= 0
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
    # Validate filter
    if not min_rt >= 0.0:
        raise exc.InvalidFilter()

    # Flexiblility for ignore cols
    ignore_cols = kwargs.get(
        "ignore_cols",
        [
            "PrecMHplus",
            "CPPIS",
            "PccChain",
            "Mode",
            "Func",
            "Scan",
            "Sequence",
            "IsVirtual",
            "FragEff",
            "IsoA0Ratio",
            "IsoA1Ratio",
            "IsoA2Ratio",
            "IsoA0RatioCV",
            "IsoA1RatioCV",
            "IsoA2RatioCV",
            "IsoA3RatioCV",
            "IsoA3Ratio",
            "ProdMHplus",
            "ProdMz",
            "ProdIntensity",
            "Ar1",
            "Ar3",
            "A0ProdMzBindex",
        ],
    )
    rt_col = kwargs.get("rt_col", "RetTime")
    path = Path(fname)

    df = pd.read_csv(path).drop(ignore_cols, errors="ignore", axis=1).drop_duplicates()

    return df[df[rt_col] > min_rt].copy()


def prep_func(fname, exp_name, min_rt=0.8, **kwargs):
    """Drop unnecessary columns and add metadata to func001-like DF

    Args:
        fname (str or Path): Path to func001-like CSV
        exp_name (str): Global experiment ID, used for splitting sample
        min_rt(float, optional): Filter RT less than this value. MUST BE >= 0

    Returns:
        pd.DataFrame: Pre-processed func001-like DF
    """
    # Validate filter
    if not min_rt >= 0.0:
        raise exc.InvalidFilter()

    # Some defaults with flexibility for kwargs
    sname_char = kwargs.get("sname_char", "_")
    sname_index = kwargs.get("sname_index", 1)
    ignore_cols = kwargs.get(
        "ignore_cols", ["drift", "DriftFwhm", "QuadMass", "FunctionIndex"]
    )
    rt_col = kwargs.get("rt_col", "RT")

    fname = Path(fname)
    sname = fname.name.split(sname_char)[sname_index]
    df = pd.read_csv(fname).drop(ignore_cols, axis=1)
    df["Organism"], df["Isotope"], df["Condition"] = split_samplename(
        sname, sname_len=len(exp_name)
    )

    return df[df[rt_col] > min_rt].copy()


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

    files = pd.DataFrame(
        {
            "path": cppis_csvs,
            "filename": [x.name for x in cppis_csvs],
            "condition": [find_condition(x.name) for x in cppis_csvs],
        }
    )
    return files


def munge_cppis(files, cond, out_dir, config=None):
    """Pre-process data for specific condition before detecting isotope labelling

    Args:
        files (pd.DataFrame): DataFrame containing list of filenames from getfilenames_cppis function
        cond (str): condition ie ACE, ACED0, BLANK, etc
        out_dir (str or Path): Directory to put outputs in

    Returns:
        pd.DataFrame: DataFrame which has been pre-processed for IsoTracing
    """
    out_dir = Path(out_dir)
    replicates = files[
        files["condition"] == cond
    ]  # get filenames associated with current condition
    dfs = [prep_cppis(f) for f in replicates["path"]]
    # All reps dataframe
    # Will be editted inplace along the way
    print("Combining DFs")
    df = combine_dfs(dfs)
    df.to_csv(out_dir.joinpath("All_ions_sorted.csv"), index=False)

    # R-tree based replicate comparison
    print("Performing replicate comparison")
    averaged = dereplicator.replicate_compare(df, config=config)
    averaged.to_csv(out_dir.joinpath("All_ions_averaged.csv"), index=False)

    return averaged


def blank_subtract(blank_df, df, inplace=True, config=None):
    """Given a DF of blank data and another DF, remove blanks from DF

    Args:
        blank_df (pd.DataFrame): Blanks data frame
        df (pd.DataFrame): Other CPPIS like DF
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
    # Make sure indices are sequential integers 0,1,2,etc...
    df.reset_index(inplace=True, drop=True)
    # Keep a set of IDs to drop
    drop_me = dereplicator.find_overlap(df, blank_df, config=config)
    # Will return None if inplace=True, else return DataFrame
    return df.drop(drop_me, inplace=inplace)


def split_samplename(s, sname_len=8):
    """Take sample name and split into organism, isotope, and condition tuple

    Args:
        s (string): Sample name in format "RLUS1234ISODN

    Returns:
        tuple: organism, isotope, condition strings tuple
    """
    i1 = sname_len
    i2 = sname_len + 2
    return s[:i1], s[i1:i2], s[i2:]


def get_func_slice(df, mz, low_scan, high_scan, mz_tol=10):
    """Return the indices of func DF given mz, scan range

    Arguments:
        df (pd.DataFrame): func001-like DataFrame to slice
        mz (float): Mass to query
        low_scan (int): LowScan
        high_scan (int): HighScan

    Returns:
        iterable: List-like of indices in func001-like DF
    """
    # print(len(seen))
    mz_tol = ppm_tolerance(mz, error=mz_tol)
    masks = (
        df["MZ"] >= mz_tol[0],
        df["MZ"] <= mz_tol[1],
        df["FunctionScanIndex"] >= low_scan,
        df["FunctionScanIndex"] <= high_scan,
        # [idx not in seen for idx in df.index],
    )
    func_slice = df[reduce(np.logical_and, masks)]
    indices = list(func_slice.index)
    # seen.update(indices)
    return indices


def calc_exp(g):
    """Calculate the avgPrecMZ, Lowscan and HighScan for an ExpId

    Args:
        g (pd.DataFrame): Exp_ID DataFrame

    Returns:
        tuple: (avgPrecMz, LowScan, HighScan)
    """
    return round(g.PrecMz.mean(), 4), g["LowScan"].min(), g["HighScan"].max()


def isotope_slicer(df, mz, low_scan, high_scan, min_scans=5, iso="C13", mz_tol=10):
    """Iteratively find all isotope data associated with a given mass.
    Continues until a slice has less than five datapoints.

    Args:
        df (pd.DataFrame): Func001 dataframe
        mz (float): Precursor mass to start scanning from
        low_scan (int): Low scan value in CPPIS
        high_scan (int): High scan value in CPPIS

    Returns:
        list: indices to mark after processing
    """
    results = []
    # Find base ion,
    # results.append(get_func_slice(df, mz, low_scan, high_scan, seen))
    results.append(get_func_slice(df, mz, low_scan, high_scan))

    # find isotopes +/- C13
    # Initialize while loop
    # Need to look forwards only because CPPIS only contains M0 peaks
    if iso == "N15":
        iso_func = n_isotope
    else:
        iso_func = c_isotope

    this_mz = mz
    while True:
        mn = iso_func(this_mz)
        this_mz = mn
        indices = get_func_slice(df, mn, low_scan, high_scan, mz_tol=10)
        results.append(indices)
        # print(f"Found {len(indices)} features")
        # Stop condition
        if len(indices) < min_scans:
            break

    return results


def mark_func(df, results):
    """Marks isotope analysis DF"""
    # data = df.copy()
    data = []
    seen = set()
    # for e_id, marks in results.items():
    # results as a list of tuples
    # each of marks is a list of indices
    for e_id, marks in results:
        hits = filter(lambda x: x, marks)
        for c, idc in enumerate(hits):
            idc_filtered = list(filter(lambda x: x not in seen, idc))
            seen.update(idc_filtered)
            data.extend(
                {"idx": i, "Isotopomer": f"M{c}", "Exp_ID": e_id} for i in idc_filtered
            )
            # data.loc[idc_filtered, "Isotopomer"] = f"M{c}"
            # data.loc[idc_filtered, "Exp_ID"] = e_id

    ddf = pd.DataFrame(data).set_index("idx")
    # return data[-data['Exp_ID'].isna()].copy()
    df1 = df.join(ddf)
    return df1[-df1["Exp_ID"].isna()].copy()


def calc_rep_stats(df, exp_id, iso, cond, min_scans=3):
    data = []

    isos = sorted(df.Isotopomer.unique(), key=lambda x: int(x.strip("M")))
    if df.empty:
        return data
    # Calculate the slope data for each replicate
    for i in range(len(isos) - 1):
        g = df.drop_duplicates(["FunctionScanIndex", "Isotopomer"]).set_index(
            "FunctionScanIndex"
        )
        mi = g[g.Isotopomer == isos[i]]
        mj = g[g.Isotopomer == isos[i + 1]]
        scans = np.intersect1d(mi.index, mj.index)
        # TODO: logging
        if len(scans) < min_scans:
            continue
        slope, intercept, _, _, std_err = linregress(
            mi.loc[scans, "Intensity"].values, mj.loc[scans, "Intensity"].values
        )

        data.append(
            {
                "Exp_ID": exp_id,
                "Condition": cond,
                "Isotope": iso,
                "Isotopomer": f"{isos[i]}v{isos[i+1]}",
                "Slope": slope,
                "Intercept": intercept,
                "std_err": std_err,
                "scan_count": len(scans),
            }
        )
    return data


def aggregate_results(df):
    ic = df.groupby(["Exp_ID", "Isotope", "Isotopomer"])["Slope"].agg(
        [("mean", "mean"), ("stdev", "std"), ("rep_count", "count")]
    )
    return pd.merge(
        df, ic, left_on=["Exp_ID", "Isotope", "Isotopomer"], right_index=True
    )


def get_cppis_masterlist(source_dir):
    source_dir = Path(source_dir)
    all_features_csv = source_dir.glob("*_All_features.csv")
    return pd.read_csv(next(all_features_csv))  # load mz masterlist (1 file)


def conf_test(t, p, alpha=0.05):
    # Simple p-test plus negative t-score (slope larger than M0)
    # Use this do define labelled/enriched
    return p <= alpha and t < 0.0


def get_summary_df(master):
    temp = master[["Exp_ID", "RetTime", "PrecMz"]].copy().set_index("Exp_ID")
    # Sort by ExpID
    temp["indexNumber"] = [int(i.split("_")[-1]) for i in temp.index]
    temp.sort_values("indexNumber", ascending=True, inplace=True)
    temp.drop("indexNumber", axis=1, inplace=True)
    # Average PrecMz, RetTime
    data = [
        {
            "Exp_ID": idx,
            "RetTime": temp.loc[idx, "RetTime"].mean(),
            "PrecMz": temp.loc[idx, "PrecMz"].mean(),
        }
        for idx in temp.index.unique()
    ]
    return pd.DataFrame(data).set_index("Exp_ID")


def condition_stats(data_dir, idxs, cond):
    """Calculate summary for a given condition
    Takes as input the `All_isotope_conditions` data dir
    """
    res_csv = data_dir.glob(f"*{cond}.csv")
    df = pd.read_csv(next(res_csv))

    unlabelled = []
    label_list = []
    label_count = []
    for idx in idxs:
        slc = df[df["Exp_ID"] == idx]
        unlabelled.append(any(x in slc.Isotope.unique() for x in (12, 14)))
        label_list.append(len(slc[slc.labelled == True]) > 0)
        label_count.append(isotope_count(slc))
    return {
        f"{cond}_unlabelled": unlabelled,
        f"{cond}_labelled": label_list,
        f"{cond}_labelled_count": label_count,
    }


def isotope_count(slc):
    # Takes slice of data with sig p-val and neg tstat
    # Sort by isotopomer ratio
    enriched = slc[(-slc.Isotope.isin([12, 14])) & (slc.labelled == True)]
    if len(enriched) < 1:
        return 0
    return max(int(i.split("v")[1].strip("M")) for i in enriched.Isotopomer)


def summarize_labels(data_dir, master, conditions):
    data_dir = Path(data_dir)
    print("Summarizing data")
    # New dataframe for data
    df = get_summary_df(master)
    idxs = list(df.index)
    for cond in conditions:
        df = df.assign(**condition_stats(data_dir, idxs, cond))
    return df


def get_mz_tol(config):
    return config.get("Tolerances", {}).get("PrecMz", ("ppm", 10))[1]


def filter_summary(df, num_conditions):
    # This function takes the summary df and returns
    # a df with only ions labeled in at least x conditions
    print(f"Filtering summary to a minimum of {num_conditions} conditions")
    idxs = df.index.values
    labelled_cols = [x for x in df.columns[df.columns.str.contains("_labelled")]]
    if num_conditions > len(labelled_cols):
        print(
            "\033[93mWARNING: num conditions for filtering is greater than number of possible conditions\033[0m"
        )
    labelled = df.loc[:, labelled_cols]
    labelled_idxs = []
    for i in idxs:
        i_labels = labelled.loc[i, labelled_cols]
        if sum(i_labels) >= num_conditions:
            labelled_idxs.append(i)
    return df.loc[labelled_idxs]


def run_label_analysis(df, cond, out_dir):
    out_dir = Path(out_dir)
    print(f"Analyzing labels in {cond}")
    unique_isotopes = df["Isotope"].unique()  # isotopes (U/L, 12/13, 14/15)
    nat_iso = unique_isotopes.min()  # unlabeled 12 or 14
    label_iso = unique_isotopes.max()  # labeled 13 or 15
    # Initialize results columns in case data is all null for a condition
    df["labelled"] = np.nan
    df["pval"] = np.nan
    df["tstat"] = np.nan
    exps = df.groupby("Exp_ID")
    for _, grp in exps:
        # if only 1 isotope condition (U/L), skip
        # or data is len 1
        if len(grp.Isotope.unique()) < 2 or len(grp) == 1:
            continue
        nat_ratio = grp.loc[
            (grp.Isotope == nat_iso) & (grp.Isotopomer == "M0vM1"), "Slope"
        ]
        # if only two instances of slope calculated, skip
        if len(nat_ratio) < 2:
            continue

        labelled_slc = grp[grp.Isotope == label_iso]
        for idx in labelled_slc.Isotopomer.unique():
            labelled_ratio = labelled_slc.loc[labelled_slc.Isotopomer == idx, "Slope"]
            t, p = ttest_ind(
                nat_ratio.values,
                labelled_ratio.values,
                equal_var=False,
                nan_policy="omit",
            )
            labelled = conf_test(t, p, alpha=0.05)
            indx = labelled_ratio.index.values
            df.loc[indx, "labelled"] = labelled
            df.loc[indx, "pval"] = p
            df.loc[indx, "tstat"] = t
    df.to_csv(out_dir.joinpath(f"iso_analysis_{cond}.csv"), index=False)
