#!/usr/bin/env python
# coding: utf-8
import re
from functools import reduce
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import linregress, ttest_ind

import isoanalyst.dereplicator as dereplicator
import isoanalyst.exceptions as exc
from isoanalyst.input_spec import InputSpec
from isoanalyst.convert import waters, mzmine

############################
##      General
############################
def ppm_tolerance(mass: float, error: float = 10.0):
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


def c_isotope(mass: float, sign: int = 1):
    """Add mass difference from 12C to 13C to mass

    Args:
        mass (float): Mass to calculate against
        sign (int): Optional: Plus or minus 1. Default=1

    Returns:
        float: Mass plus (minus if sign = -1) mass difference from 12C to 13C
    """
    return mass + sign * 1.00335


def n_isotope(mass: float, sign: int = 1):
    """Add mass difference from 14N to 15N to mass
    14N = 14.003074
    15N = 15.000109

    Args:
        mass (float): Mass to calculate against
        sign (int): Optional: Plus or minus 1. Default=1

    Returns:
        float: Mass plus (minus if sign = -1) mass difference from 14N to 15N
    """
    return mass + sign * 0.9970349


def combine_dfs(dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """Combine any number of pandas dataframe

    Args:
        dfs (list): List of pd.DataFrame objects

    Returns:
        pd.DataFrame: A combined pandas dataframe
    """
    return pd.concat(dfs, sort=False).reset_index(drop=True)


def get_featurelist(source_dir: Path, exp_name: str):
    all_features_csv = source_dir / f"{exp_name}_all_features.csv"
    return pd.read_csv(all_features_csv)


############################
##      Prep
############################
def prep_featurefile(fname: Path, min_rt=0.8, **kwargs):
    """Take feature file and convert / get necessary data

    Args:
        fname (str or Path): Path to cppis type dataframe
        min_rt(float, optional): Filter RT less than this value. MUST BE >= 0
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
    # Validate filter
    if not min_rt >= 0.0:
        raise exc.InvalidFilter()

    rt_col = kwargs.get("rt_col", "rettime")
    if "cppis" in fname.name.lower():
        df = waters.cppis(fname)
    else:
        df = mzmine.feature_list(fname)

    return df[df[rt_col] > min_rt].copy()


def munge_featurelist(inp_spec: InputSpec, cond: str, out_dir: Path, config: Dict):
    """Pre-process data for specific condition before generating feature list.

    Args:
        inp_spec: InputSpec
        cond (str): condition ie ACE, ACED0, BLANK, etc
        out_dir (str or Path): Directory to put outputs in

    Returns:
        pd.DataFrame: DataFrame which has been pre-processed for IsoTracing
    """
    out_dir = Path(out_dir)
    replicates = inp_spec.get_feature_filepaths(cond)
    # All reps dataframe
    # Will be editted inplace along the way
    dfs = [prep_featurefile(Path(f)) for f in replicates]
    print(f"Combining DFs for {cond}")
    try:
        df = combine_dfs(dfs)
    except ValueError:
        print(f"No features found for {cond}")
        # create empty dataframe
        df = pd.DataFrame(columns=["precmz", "rettime", "precz", "sample"])
        df.to_csv(out_dir.joinpath("all_ions_sorted.csv"), index=False)
        averaged = pd.DataFrame(
            columns=["rettime", "precmz", "precz", "samples", "rep_count"]
        )
        averaged.to_csv(out_dir.joinpath("all_ions_averaged.csv"), index=False)
        return df

    df.to_csv(out_dir.joinpath("all_ions_sorted.csv"), index=False)

    # R-tree based replicate comparison
    print("Performing replicate comparison")
    averaged = dereplicator.replicate_compare(df, config=config)
    averaged.to_csv(out_dir.joinpath("all_ions_averaged.csv"), index=False)

    return averaged


def blank_subtract(blank_df, df, config, inplace=True):
    """Given a DF of blank feature list data and another feature list DF,
    remove blanks from DF

    Args:
        blank_df (pd.DataFrame): Blanks data frame
        df (pd.DataFrame): Other featurelist DF
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
    # Make sure indices are sequential integers 0,1,2,etc...
    df.reset_index(inplace=True, drop=True)
    # Keep a set of IDs to drop
    drop_me = dereplicator.find_overlap(df, blank_df, config=config)
    # Will return None if inplace=True, else return DataFrame
    return df.drop(drop_me, inplace=inplace)


############################
##      Scrape
############################
def prep_scan(
    fname: Path,
    inp_spec: InputSpec,
    min_intensity: int,
    min_rt: float = 0.8,
    **kwargs,
):
    """Drop unnecessary columns and add metadata to all scan DF

    Args:
        fname (Path): Path to all scan CSV
        inp_spec (InputSpec): Input specification
        min_intensity(int): Filter data less than this value. MUST BE >= 0
        min_rt(float, optional): Filter RT less than this value. MUST BE >= 0

    Returns:
        pd.DataFrame: Pre-processed all scan DF
    """
    # Validate filter
    if not min_rt >= 0.0 or not min_intensity >= 0:
        raise exc.InvalidFilter()

    # Some defaults with flexibility for kwargs
    rt_col = kwargs.get("rt_col", "rettime")

    if fname.suffix.lower() == ".mzml":
        df = waters.mzml(fname, min_intensity=min_intensity)
    elif "func001" in fname.name:
        df = waters.func001(fname, min_intensity=min_intensity)
    else:
        raise Exception(f"DATA MISSING? - {fname}")

    meta = inp_spec.get_filepath_info(fname)
    df["organism"] = meta["organism"]
    df["isotope"] = meta["isotope"]
    df["condition"] = f'{meta["condition"]}-{meta["replicate"]}'

    return df[df[rt_col] > min_rt].copy()


def blank_subtract(blank_df, df, config, inplace=True):
    """Given a DF of blank feature list data and another feature list DF,
    remove blanks from DF

    Args:
        blank_df (pd.DataFrame): Blanks data frame
        df (pd.DataFrame): Other featurelist DF
        inplace (bool, optional): Edit dataframe in place. Defaults to True.
    """
    # Make sure indices are sequential integers 0,1,2,etc...
    df.reset_index(inplace=True, drop=True)
    # Keep a set of IDs to drop
    drop_me = dereplicator.find_overlap(df, blank_df, config=config)
    # Will return None if inplace=True, else return DataFrame
    return df.drop(drop_me, inplace=inplace)


def get_scan_slice(
    df: pd.DataFrame, mz: float, low_scan: int, high_scan: int, mz_tol: float = 10.0
):
    """Return the indices of scan DF given mz, scan range

    Arguments:
        df (pd.DataFrame): scan DataFrame to slice
        mz (float): Mass to query
        low_scan (int): lowscan
        high_scan (int): highscan
        mz_tol(number): M/Z tolerance in PPM

    Returns:
        iterable: List-like of indices in scan DF
    """
    mz_tol = ppm_tolerance(mz, error=mz_tol)
    masks = (
        df["mz"] >= mz_tol[0],
        df["mz"] <= mz_tol[1],
        df["scanindex"] >= low_scan,
        df["scanindex"] <= high_scan,
        # [idx not in seen for idx in df.index],
    )
    func_slice = df[reduce(np.logical_and, masks)]
    indices = list(func_slice.index)
    return indices


def calc_exp(g):
    """Calculate the avgprecmz, Lowscan and highscan for an ExpId

    Args:
        g (pd.DataFrame): exp_id DataFrame

    Returns:
        tuple: (avgprecmz, lowscan, highscan)
    """
    return round(g.precmz.mean(), 4), g["lowscan"].min(), g["highscan"].max()


def isotope_slicer(
    df: pd.DataFrame,
    mz: float,
    low_scan: int,
    high_scan: int,
    min_scans: int = 5,
    mz_tol: float = 10.0,
):
    """Iteratively find all isotope data associated with a given mass.
    Continues until a slice has less than five datapoints.

    Args:
        df (pd.DataFrame): Func001 dataframe
        mz (float): Precursor mass to start scanning from
        low_scan (int): Low scan value in CPPIS
        high_scan (int): High scan value in CPPIS
        mz_tol (number): mz tolerance in PPM
    Returns:
        list: indices to mark after processing
    """
    results = []
    # Find base ion,
    # results.append(get_func_slice(df, mz, low_scan, high_scan, seen))
    results.append(get_scan_slice(df, mz, low_scan, high_scan, mz_tol=mz_tol))

    # find isotopes +/- 13C or 15N
    # Initialize while loop
    # Need to look forwards only because CPPIS only contains M0 peaks
    iso = df.iloc[0]["isotope"]
    if iso == 15:
        iso_func = n_isotope
    else:
        iso_func = c_isotope

    this_mz = mz
    while True:
        mn = iso_func(this_mz)
        this_mz = mn
        indices = get_scan_slice(df, mn, low_scan, high_scan, mz_tol=mz_tol)
        results.append(indices)
        # print(f"Found {len(indices)} features")
        # Stop condition
        if len(indices) < min_scans:
            break

    return results


def mark_scans(df, results):
    """Marks isotope analysis DF"""
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
                {"idx": i, "isotopomer": f"M{c}", "exp_id": e_id} for i in idc_filtered
            )
            # data.loc[idc_filtered, "isotopomer"] = f"M{c}"
            # data.loc[idc_filtered, "exp_id"] = e_id

    ddf = pd.DataFrame(data).set_index("idx")
    # return data[-data['exp_id'].isna()].copy()
    df1 = df.join(ddf)
    return df1[-df1["exp_id"].isna()].copy()


def add_scan_window(feature_df: pd.DataFrame, scan_df: pd.DataFrame, scanwindow: int):
    """For each feature find closest scan index to rettime in scan_df, then add window"""
    df = feature_df.copy()
    lowscan = []
    highscan = []
    for _, row in df.iterrows():
        closest = scan_df.iloc[(scan_df["rettime"] - row.rettime).abs().argsort()[0]]
        low, high = closest.scanindex - scanwindow, closest.scanindex + scanwindow
        lowscan.append(low)
        highscan.append(high)
    df["lowscan"] = lowscan
    df["highscan"] = highscan
    return df


############################
##      Analyze
############################
def calc_rep_stats(df, exp_id, iso, cond, min_scans=3):
    data = []

    isos = sorted(df.isotopomer.unique(), key=lambda x: int(x.strip("M")))
    if df.empty:
        return data
    # Calculate the slope data for each replicate
    for i in range(len(isos) - 1):
        g = df.drop_duplicates(["scanindex", "isotopomer"]).set_index("scanindex")
        mi = g[g.isotopomer == isos[i]]
        mj = g[g.isotopomer == isos[i + 1]]
        scans = np.intersect1d(mi.index, mj.index)
        # TODO: logging
        if len(scans) < min_scans:
            continue
        slope, intercept, _, _, std_err = linregress(
            mi.loc[scans, "intensity"].values, mj.loc[scans, "intensity"].values
        )

        data.append(
            {
                "exp_id": exp_id,
                "condition": cond,
                "isotope": iso,
                "isotopomer": f"{isos[i+1]}v{isos[i]}",
                "slope": slope,
                "intercept": intercept,
                "std_err": std_err,
                "scan_count": len(scans),
            }
        )
    return data


def aggregate_results(df):
    ic = df.groupby(["exp_id", "isotope", "isotopomer"])["slope"].agg(
        [("mean", "mean"), ("stdev", "std"), ("rep_count", "count")]
    )
    return pd.merge(
        df, ic, left_on=["exp_id", "isotope", "isotopomer"], right_index=True
    )


def confidence_test(t, p, alpha=0.05):
    # Simple p-test plus negative t-score (slope larger than M0)
    # Use this do define labelled/enriched
    return p <= alpha and t < 0.0


def generate_summary_df(features):
    temp = features[["exp_id", "rettime", "precmz", "precz"]].copy().set_index("exp_id")
    # Sort by ExpID
    temp["indexnumber"] = [int(i.split("_")[-1]) for i in temp.index]
    temp.sort_values("indexnumber", ascending=True, inplace=True)
    temp.drop("indexnumber", axis=1, inplace=True)
    # Average precmz, rettime
    data = [
        {
            "exp_id": idx,
            "rettime": temp.loc[idx, "rettime"].mean(),
            "precmz": temp.loc[idx, "precmz"].mean(),
            "precz": temp.loc[idx, "precz"].mean(),
        }
        for idx in temp.index.unique()
    ]
    return pd.DataFrame(data).set_index("exp_id")


def condition_stats(data_dir, idxs, cond):
    """Calculate summary for a given condition
    Takes as input the `all_isotope_conditions` data dir
    """
    res_csv = data_dir.glob(f"*{cond}.csv")
    df = pd.read_csv(next(res_csv))

    unlabelled = []
    label_list = []
    label_count = []
    for idx in idxs:
        slc = df[df["exp_id"] == idx]
        unlabelled.append(any(x in slc.isotope.unique() for x in (12, 14)))
        label_list.append(len(slc[slc.labelled == True]) > 0)
        label_count.append(isotope_count(slc))
    return {
        f"{cond}_unlabelled": unlabelled,
        f"{cond}_labelled": label_list,
        f"{cond}_count": label_count,
    }


def isotope_count(slc):
    # Takes slice of data with sig p-val and neg tstat
    # Sort by isotopomer ratio
    enriched = slc[(-slc.isotope.isin([12, 14])) & (slc.labelled == True)]
    if len(enriched) < 1:
        return 0
    return max(int(i.split("v")[0].strip("M")) for i in enriched.isotopomer)


def summarize_labels(data_dir, features, conditions):
    data_dir = Path(data_dir)
    print("Summarizing data")
    # New dataframe for data
    df = generate_summary_df(features)
    idxs = list(df.index)
    for cond in conditions:
        df = df.assign(**condition_stats(data_dir, idxs, cond))
    return df


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
    unique_isotopes = df["isotope"].unique()  # isotopes (U/L, 12/13, 14/15)
    nat_iso = unique_isotopes.min()  # unlabeled 12 or 14
    label_iso = unique_isotopes.max()  # labeled 13 or 15
    # Initialize results columns in case data is all null for a condition
    df["labelled"] = np.nan
    df["pval"] = np.nan
    df["tstat"] = np.nan
    exps = df.groupby("exp_id")
    for _, grp in exps:
        # if only 1 isotope condition (U/L), skip
        # or data is len 1
        if len(grp.isotope.unique()) < 2 or len(grp) == 1:
            continue
        nat_ratio = grp.loc[
            (grp.isotope == nat_iso) & (grp.isotopomer == "M1vM0"), "slope"
        ]
        # if only two instances of slope calculated, skip
        if len(nat_ratio) < 3:
            continue

        labelled_slc = grp[grp.isotope == label_iso]
        for idx in labelled_slc.isotopomer.unique():
            labelled_ratio = labelled_slc.loc[labelled_slc.isotopomer == idx, "slope"]
            if len(labelled_ratio) < 3:
                continue
            t, p = ttest_ind(
                nat_ratio.values,
                labelled_ratio.values,
                equal_var=False,
                nan_policy="omit",
            )
            labelled = confidence_test(t, p, alpha=0.05)
            indx = labelled_ratio.index.values
            df.loc[indx, "labelled"] = labelled
            df.loc[indx, "pval"] = p
            df.loc[indx, "tstat"] = t
    df.to_csv(out_dir.joinpath(f"iso_analysis_{cond}.csv"), index=False)
