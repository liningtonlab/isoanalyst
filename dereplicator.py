import json

import numpy as np
import pandas as pd
import rtree

CONFIG = {
    "ColsToMatch": ["RetTime", "PrecMz", "PrecZ"],
    "Tolerances": {
        "PrecMz": ["ppm", 10],
        "RetTime": ["window", 0.03],
        "PrecZ": ["window", 0.1],
    },
    "MinReps": 3,
}

def replicate_compare(df, config=CONFIG):
    """Take dataframe and perform Rtree comparison to replicate
    
    Args:
        df (pd.DataFrame): Data to work on
    
    Returns:
        pd.DataFrame: New dataframe with averaged data
    """
    gen_error_cols(df, config['Tolerances'])
    rects = get_rects(df, config)
    rtree = build_rtree(rects, config)
    con_comps = gen_con_comps(rtree, rects)
    file_col = []
    new_data = []
    for c in con_comps:
        if len(c) > 1:
            c_df = df.iloc[list(c)]
            unique_samples = set(c_df['Sample'].values)
            if len(unique_samples) >= config['MinReps']:
                new_data.append(collapse_data_rows(c_df, config['ColsToMatch']))
    return pd.DataFrame(new_data).round(4)


def gen_error_cols(df, errorinfo):
    """
    Uses the errorinfo dict to generate
    error windows for each of the columns.
    Mutates dataframe inplace for some memory conservation.
    possible error types are:
    * ppm - parts per million
    * perc - percentage
    * factor - a multiplier (ie 10 = 10x)
    * window - a standard fixed error window

    Args:
        df (pandas.DataFrame): input dataframe to calc error windows (modified in place)
        errorinfo (dict): dict of error information
    """

    for dcol, einfo in errorinfo.items():
        col = df[dcol]
        etype, evalue = einfo
        if etype == 'ppm':
            efunc = lambda x:x*(evalue*1e-6)
        if etype == 'perc':
            efunc = lambda x:x*(evalue/100)
        if etype == 'factor':
            efunc = lambda x:x*evalue
        if etype == 'window':
            efunc = lambda x:evalue
        if etype is None:
            efunc = lambda x:0
        errors = col.apply(efunc)
        df[f"{dcol}_low"] = df[dcol] - errors
        df[f"{dcol}_high"] = df[dcol] + errors


def gen_con_comps(rtree: rtree.index.Index, rects: np.ndarray) -> set:
    """
    Generate connected components subgraphs for a graph where nodes are hyperrectangles
    and edges are overlapping hyperrectangles. This is done using the rtree index and
    a depth first search.
    """
    seen = set()

    for i, _ in enumerate(rects):
        if i in seen:
            continue
        search_idxs = [i]
        c = {i}
        while search_idxs:
            search = search_idxs.pop()
            try:
                neighbors = set(rtree.intersection(rects[search]))
            except Exception as e:
                print(e)
                print(rects[search])
                raise e

            for n in neighbors - seen:  # set math
                c.add(n)
                search_idxs.append(n)
                seen.add(n)
        yield c


def build_rtree(rects: np.ndarray, config) -> rtree.index.Index:
    """
    Build RTree index for rectangles for fast range queries.
    df needs errors cols pre-calculated
    """
    dims = len(config["ColsToMatch"])
    p = rtree.index.Property()
    p.dimension = dims
    p.interleaved = False
    rgen = ((i, r, None) for i, r in enumerate(rects))
    idx = rtree.index.Index(rgen, properties=p)
    return idx


def get_rects(df: pd.DataFrame, config) -> np.ndarray:
    """
    Get the error portions of df
    """
    ecols = [f"{c}_low" for c in config["ColsToMatch"]]
    ecols = ecols + [f"{c}_high" for c in config["ColsToMatch"]]

    return df[ecols].values


def collapse_data_rows(df: pd.DataFrame, datacols: list, calc_bin_info: bool=False) -> dict:
    """
    Takes conncect component DF and return average of compared values
    as dict for appending to list for new DF construction
    """
    unique=set(df["Sample"].values)
    data = {k: df[k].mean() for k in datacols}
    if calc_bin_info:
        bin_info = {cn:[float(df[cn].min()),float(df[cn].max())] for cn in datacols}
        bin_info['n'] = len(df)
        data['bin_info'] = json.dumps(bin_info)

    data["Samples"] = "|".join(unique)
    data["LowScan"] = df["ScanLowRange"].min()
    data["HighScan"] = df["ScanHighRange"].max()
    data["rep_count"] = len(unique)
    data['RetTime'] = round(data['RetTime'], 3)
    return data


def find_overlap(df1, df2, config=CONFIG):
    """Find overlap of df1 in df2 and return set of iloc indices from df1
    
    Args:
        df1 (pd.DataFrame): Query DataFrame (Are any of DF1 rects in DF2)
        df2 (pd.DataFrame): Reference DataFrame
        config (dict, optional): Configuration dict. Defaults to CONFIG.

    Returns:
        set: Set of iloc indices of df1 
    """
    # Make temp copies of dataframes
    query = df1.copy()
    ref = df2.copy()
    gen_error_cols(ref, config['Tolerances'])
    gen_error_cols(query, config['Tolerances'])
    
    ref_rects = get_rects(ref, config)
    query_rects = get_rects(query, config)
    tree = build_rtree(query_rects, config)
    overlap = set()
    for i, r in enumerate(ref_rects):
        overlap.update(tree.intersection(r))
    return overlap
