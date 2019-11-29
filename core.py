#!/usr/bin/env python
# coding: utf-8
# import multiprocessing
import joblib
from pathlib import Path
import pandas as pd
from utils import (blank_subtract, combine_dfs, getfilenames_cppis,
                   group_cppis, munge_cppis, prep_func, isotope_slicer)
import time


def cppis_masterlist(source_dir, conditions):
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

    print("Setting up multi-processing...")
    futures = joblib.Parallel(8)(joblib.delayed(munge_condition)(c) for c in conditions)

    # Collect futures results into list
    prim_dfs = list(futures) 
    all_primary = combine_dfs(prim_dfs)
    print("Substracting blanks")
    blank_subtract(blanks, all_primary) # blank subtraction on full dataset 
    print("Grouping all features")
    group_cppis(all_primary, Exp_name, 'Exp_ID') # final grouping - Exp_ID used for func001 munging
    all_primary.to_csv(source_dir.joinpath(f'{Exp_name}_All_features.csv'), index=False)
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

    uniques = master.Exp_ID.unique() # list of unique features from mz masterlist

    # # Why create and write header for each condition?
    # cols = ['FunctionIndex', 'FunctionScanIndex', 'RT', 'MZ', 'mzFWHM', 'Counts',
    #    'Intensity', 'Organism', 'Isotope', 'Condition', 'Isotopomer', 'Exp_ID']
    # header = pd.DataFrame(columns = cols) # header to be used for data file for each condition

    # dir for all output files containing scan data for each cond
    out_dir = source_dir.joinpath('All_scan_data')
    out_dir.mkdir(parents=True, exist_ok=True)

    def run_isoslicer(c):
         # header.to_csv(out_dir.joinpath(f'All_ions_{c}.csv'), index = False) # make file to add iso data to
        
        # add replicate func files for current condition
        c_funcs = filter(lambda f: c in f.name, func_files)
        c_dfs = [prep_func(f) for f in c_funcs]

        # merged data frame of all func001s for the given condition 'c'
        func_df = combine_dfs(c_dfs) 
        
        # loops through all unique ions in the master list and iteratively adds
        # Set of seen indices in func_df, updated in isotope_slicer
        seen = set()
        print(f"There are {len(uniques)} ions to look at")
        for i, u in enumerate(uniques):
            if i % 10 == 0:
                print(f"Working on {i}/{len(uniques)}")

            exp = master[master.Exp_ID==u]
            # Don't bother when Z > 1 for now
            if exp.PrecZ.iloc[0] > 1:
                continue

            mz = round(exp.PrecMz.mean(), 4) # average those masses
            low_scan = exp['LowScan'].min() # get lowest scan
            high_scan = exp['HighScan'].max() # get highest scan
            
            #  function 'iso_slicer' slices relevant isotope data for a given mz and all its isotopomers 
            #  within given scan range in given func file
            start = time.time()
            isotope_slicer(func_df, mz, low_scan, high_scan, u, seen)
            end = time.time()
            print(f"IsoSlicer took {end-start} seconds... for {u}")
        # write iso_scan_data to a file for that condition
        func_df.to_csv(out_dir.joinpath(f'All_ions_{c}.csv'), index=False)

    for c in conditions:
        run_isoslicer(c)



if __name__ == "__main__":
    '''User provided arguments'''

    # Input directory

    # source_dir = '/Volumes/TimeMachineBackups/201908_RLUS1353/20190808/seven/'
    # source_dir = Path("./IsotopedataProcessing/J1").absolute()
    source_dir = Path("./IsotopedataProcessing/J3").absolute()


    # Experiment name or organism name, applies to whole dataset
    # for naming output files
    Exp_name = 'RLUS1353'


    # Experimental conditions
    # primary - positional argument always a list of labels used preferably any 3 letter code per label
    # secondary - optional argument for a nested condition within the labels
    #             ie days, media, additive, etc
    # Depending on which arguments are given, a list called 'conditions' is generated and used throughout 
    # all steps of the code including file names and directory/folder naming. 

    sec_arg = 0

    if sec_arg == 1:
        primary = ['SB','SA','GB','GA']
        secondary = ['N1','N2','N3','N4']
        conditions = []
        for p in primary:
            for s in secondary:
                conditions.append(p+s) 
                
    else:
        primary = ['ACE','PROP','MET','GLU']
        conditions = primary

    master = None # for continuing
    master = cppis_masterlist(source_dir, conditions)
    # isotope_scraper(source_dir, conditions, master=master)
