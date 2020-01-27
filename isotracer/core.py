#!/usr/bin/env python
# coding: utf-8
from pathlib import Path

import joblib
import pandas as pd
from scipy.stats import ttest_ind

import isotracer.dereplicator as dereplicator
import isotracer.utils as utils


def cppis_masterlist(source_dir, conditions, exp_name):
    """
    From a directory takes a folder named 'CPPIS' and munges all the cppis.csv files first by conditions
    and then combines all files to create an mz masterlist of all features detected in experiment and basketed
    to align between conditions and replicates.
    """
    source_dir = Path(source_dir) # Make sure is Path object
    cppis_dir = source_dir.joinpath("CPPIS")
    print(f"Collecting files from {cppis_dir}")
    cppis_files = utils.getfilenames_cppis(cppis_dir, conditions) # df containing paths, filenames and conditions

    print(f"Working on BLANKS")
    blank_dir = source_dir.joinpath("BLANKS")
    # make directory if not exists
    blank_dir.mkdir(parents=True, exist_ok=True)

    # peaks in blanks to subtract later
    blanks = utils.munge_cppis(cppis_files, 'BLANK', blank_dir)

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
        all_collapsed = utils.munge_cppis(cppis_files, cond, out_dir)
        return all_collapsed

    # Run pre-processing on conditions in separate processes
    prim_dfs = joblib.Parallel(n_jobs=-1)(joblib.delayed(munge_condition)(c) for c in conditions)

    all_primary = utils.combine_dfs(prim_dfs)
    print("Substracting blanks")
    utils.blank_subtract(blanks, all_primary) # blank subtraction on full dataset
    print("Grouping all features")
    all_primary.reset_index(inplace=True, drop=True)
    dereplicator.group_features(all_primary, exp_name, 'Exp_ID') # final grouping - Exp_ID used for func001 munging
    all_primary.to_csv(source_dir.joinpath(f'{exp_name}_All_features.csv'), index=False)
    return all_primary


def isotope_scraper(source_dir, conditions, exp_name, master=None, n_jobs=-1):
    '''This is the final function which will scrape all the isotope data for all ions in the
    '''
    print("Running isotope scraper")
    source_dir = Path(source_dir)
    # Enables passing DF instead of loading CSV
    if not isinstance(master, pd.DataFrame):
        master = utils.get_cppis_masterlist(source_dir)

    # TODO: Add flexibility
    func_files = list(source_dir.joinpath('func001').glob("*.csv"))

    # dir for all output files containing scan data for each cond
    out_dir = source_dir.joinpath('All_scan_data')
    out_dir.mkdir(parents=True, exist_ok=True)

    def run_isoslicer(cond, restart=False):
        outfile = out_dir.joinpath(f'All_ions_{cond}.csv')
        if not restart and outfile.exists():
            print(f"{cond} already processed")
            return 

        # add replicate func files for current condition
        c_funcs = filter(lambda f: cond in f.name, func_files)
        c_dfs = (utils.prep_func(f, exp_name) for f in c_funcs)

        # merged data frame of all func001s for the given condition 'c'
        print(f"Combining files for {cond}")
        func_df = utils.combine_dfs(c_dfs)

        # loops through all unique ions in the master list and iteratively adds
        # Set of seen indices in func_df, updated in isotope_slicer
        func_df['Isotopomer'] = None
        func_df['Exp_ID'] = None
        seen = set()
        len_uni = len(master.Exp_ID.unique())
        exps = master.groupby("Exp_ID")
        print(f"There are {len_uni} ions to look at for {cond}")

        # to_mark = {}
        # for i, (idx, g) in enumerate(exps):
        #     if i % 20 == 0 and i > 0:
        #         print(f"Working on {i}/{len_uni} for {cond}")
        #     mz, low_scan, high_scan = utils.calc_exp(g)
        #     #  function 'iso_slicer' slices relevant isotope data for a given mz and all its isotopomers
        #     #  within given scan range in given func file
        #     to_mark[idx] = utils.isotope_slicer(func_df, mz, low_scan, high_scan, seen)

        # # Parallelized solution
        def do_slice(i, idx, g):
            # nonlocal seen
            if i % 20 == 0 and i > 0:
                print(f"Working on {i}/{len_uni} for {cond}")

            mz, low_scan, high_scan = utils.calc_exp(g)
            #  function 'iso_slicer' slices relevant isotope data for a given mz and all its isotopomers
            #  within given scan range in given func file
            # return idx, utils.isotope_slicer(func_df, mz, low_scan, high_scan, seen)
            return idx, utils.isotope_slicer(func_df, mz, low_scan, high_scan)

        # Will run slicing in multiple threads with shared memory for seen set
        to_mark = [ (k, v) for k,v in
            # joblib.Parallel(n_jobs=-1, require='sharedmem')(joblib.delayed(do_slice)(i, idx, g) for i, (idx, g) in enumerate(exps))
            joblib.Parallel(n_jobs=-1, prefer="threads")(joblib.delayed(do_slice)(i, idx, g) for i, (idx, g) in enumerate(exps))
        ]

        print(f"Finished collecting ions for {cond}")
        print("Preparing output")
        res_df = utils.mark_func(func_df, to_mark)

        # write iso_scan_data to a file for that condition
        res_df.to_csv(outfile)

    # Run processing of each condition in separate process
    joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(run_isoslicer)(c) for c in conditions)
    # [run_isoslicer(c) for c in conditions]


def isotope_label_detector(source_dir, conditions, master=None, n_jobs=-1):
    source_dir = Path(source_dir)
    scan_dir = source_dir.joinpath("All_scan_data")
    slope_dir = source_dir.joinpath("All_slope_data")
    slope_dir.mkdir(parents=True, exist_ok=True)
    out_dir = source_dir.joinpath("All_isotope_analysis")
    out_dir.mkdir(parents=True, exist_ok=True)
    # Enables passing DF instead of loading CSV
    if not isinstance(master, pd.DataFrame):
        master = utils.get_cppis_masterlist(source_dir)

    def run_label_analysis(df, cond):
        print(f"Analyzing labels in {cond}")
        unique_isotopes = df['Isotope'].unique() # isotopes (U/L, 12/13, 14/15)
        nat_iso = unique_isotopes.min() # unlabeled 12 or 14
        label_iso = unique_isotopes.max() # labeled 13 or 15

        exps = df.groupby("Exp_ID")
        for _, grp in exps:
            # if only in one isotope condition (U/L), skip
            # if only one instance of slope calculated, skip
            if len(grp.Isotope.unique()) < 2 or grp.shape[0] == 1:
                continue
            nat_ratio = grp.loc[(grp.Isotope==nat_iso)&(grp.Isotopomer=="M0vM1"), "Slope"]
            if nat_ratio.shape[0] < 1:
                continue

            labelled_slc = grp[grp.Isotope==label_iso]
            for idx in labelled_slc.Isotopomer.unique():
                labelled_ratio = labelled_slc.loc[labelled_slc.Isotopomer==idx, "Slope"]
                t, p = ttest_ind(nat_ratio.values, labelled_ratio.values,
                                 equal_var=False, nan_policy ='omit')
                labelled = utils.conf_test(t, p, alpha=0.05)
                indx = labelled_ratio.index.values
                df.loc[indx, "labelled"] = labelled
                df.loc[indx, "pval"] = p
                df.loc[indx, "tstat"] = t

        df.to_csv(out_dir.joinpath(f"iso_analysis_{cond}.csv"), index=False)

    def run_label_detector(cond):
        out_file = slope_dir.joinpath(f"All_slope_data_{cond}.csv")
        print(f"Detecting labels in {cond}")
        fil = scan_dir.joinpath(f"All_ions_{cond}.csv")
        assert fil.exists()
        df = pd.read_csv(fil)
        grouped = df.groupby(["Exp_ID", "Isotope", "Condition"])
        data = []
        for (e_id, iso, c), g in grouped:
            res = utils.calc_rep_stats(g, e_id, iso, c)
            # print(data)
            if len(res) < 1:
                continue
            data.extend(res)
        res_df = utils.aggregate_results(pd.DataFrame(data))
        res_df.to_csv(out_file, index=False)
        run_label_analysis(res_df, cond)


    # Run processing of each condition in separate process
    joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(run_label_detector)(c) for c in conditions)
    # [run_label_detector(c) for c in conditions]
    sum_df = utils.summarize_labels(out_dir, master, conditions)
    sum_df.to_csv(source_dir.joinpath("data_summary.csv"))

