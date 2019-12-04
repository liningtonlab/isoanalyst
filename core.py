#!/usr/bin/env python
# coding: utf-8
from pathlib import Path

import joblib
import pandas as pd

import dereplicator
from utils import (blank_subtract, calc_exp, combine_dfs, getfilenames_cppis,
                   isotope_slicer, munge_cppis, prep_func)


def cppis_masterlist(source_dir, conditions, exp_name):
    """
    From a directory takes a folder named 'CPPIS' and munges all the cppis.csv files first by conditions
    and then combines all files to create an mz masterlist of all features detected in experiment and basketed
    to align between conditions and replicates.
    """
    source_dir = Path(source_dir) # Make sure is Path object
    cppis_dir = source_dir.joinpath("CPPIS")
    print(f"Collecting files from {cppis_dir}")
    cppis_files = getfilenames_cppis(cppis_dir, conditions) # df containing paths, filenames and conditions

    print(f"Working on BLANKS")
    blank_dir = source_dir.joinpath("BLANKS")
    # make directory if not exists
    blank_dir.mkdir(parents=True, exist_ok=True)

    # peaks in blanks to subtract later
    blanks = munge_cppis(cppis_files, 'BLANK', blank_dir)

    # TODO: Re-do secondary condition preparation and processing
    # if sec_arg == 1:
    #     for p in primary:
    #         sec_dfs = [] #dfs to combine for each primary condition
    #         for s in secondary:
    #             cond = p+s
    #             out_dir = source_dir.joinpath(p, s) # create directory
    #             out_dir.mkdir(parents=True, exist_ok=True)

    #             all_collapsed = munge_cppis(cppis_files, cond, out_dir) #write files
    #             sec_dfs.append(all_collapsed) # add df to list to combine with other sec conds
    #         all_secondary = combine_dfs(sec_dfs)
    #         grouped_secondary = group_cppis(all_secondary, p, 'group_ID') # grouping all features for primary condition
    #         grouped_secondary.to_csv(source_dir.joinpath(p, f'{p}_All_ions.csv', index=False))
    #         prim_dfs.append(grouped_secondary) #add df to list to combine with other primary conds

    # else:
    def munge_condition(cond):
        out_dir = source_dir.joinpath(cond)
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Working on {cond}")
        all_collapsed = munge_cppis(cppis_files, cond, out_dir)
        return all_collapsed

    # Run pre-processing on conditions in separate processes
    prim_dfs = joblib.Parallel(n_jobs=-1)(joblib.delayed(munge_condition)(c) for c in conditions)

    all_primary = combine_dfs(prim_dfs)
    print("Substracting blanks")
    blank_subtract(blanks, all_primary) # blank subtraction on full dataset
    print("Grouping all features")
    all_primary.reset_index(inplace=True, drop=True)
    dereplicator.group_features(all_primary, exp_name, 'Exp_ID') # final grouping - Exp_ID used for func001 munging
    all_primary.to_csv(source_dir.joinpath(f'{exp_name}_All_features.csv'), index=False)
    return all_primary


def isotope_scraper(source_dir, conditions, master=None):
    '''This is the final function which will scrape all the isotope data for all ions in the
    '''
    print("Running isotope scraper")
    source_dir = Path(source_dir)
    # Find All feature CSV without knowing name
    # Enables passing DF instead of loading CSV
    if not isinstance(master, pd.DataFrame):
        all_features_csv = source_dir.glob('*_All_features.csv')
        master = pd.read_csv(next(all_features_csv)) # load mz masterlist (1 file)

    # TODO: Add flexibility
    func_files = list(source_dir.joinpath('func001').glob("*.csv"))

    # dir for all output files containing scan data for each cond
    out_dir = source_dir.joinpath('All_scan_data')
    out_dir.mkdir(parents=True, exist_ok=True)

    def run_isoslicer(c):
        # add replicate func files for current condition
        c_funcs = filter(lambda f: c in f.name, func_files)
        c_dfs = (prep_func(f) for f in c_funcs)

        # merged data frame of all func001s for the given condition 'c'
        print(f"Combining files for {c}")
        func_df = combine_dfs(c_dfs)

        # loops through all unique ions in the master list and iteratively adds
        # Set of seen indices in func_df, updated in isotope_slicer
        func_df['Isotopomer'] = None
        func_df['Exp_ID'] = None
        seen = set()
        len_uni = len(master.Exp_ID.unique())
        exps = master.groupby("Exp_ID")
        print(f"There are {len_uni} ions to look at for {c}")

        def do_slice(i, idx, g):
            nonlocal seen
            if i % 10 == 0:
                print(f"Working on {i}/{len_uni} for {c}")

            mz, low_scan, high_scan = calc_exp(g)
            #  function 'iso_slicer' slices relevant isotope data for a given mz and all its isotopomers
            #  within given scan range in given func file
            isotope_slicer(func_df, mz, low_scan, high_scan, idx, seen)

        # Will run slicing in multiple threads with shared memory for seen set
        joblib.Parallel(n_jobs=-1, require='sharedmem')(joblib.delayed(do_slice)(i, idx, g) for i, (idx, g) in enumerate(exps))

        # write iso_scan_data to a file for that condition
        func_df.to_csv(out_dir.joinpath(f'All_ions_{c}.csv'), index=False)

    # Run processing of each condition in separate process
    joblib.Parallel(n_jobs=-1)(joblib.delayed(run_isoslicer)(c) for c in conditions)
