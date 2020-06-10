#!/usr/bin/env python
# coding: utf-8
from pathlib import Path

import joblib
import pandas as pd

import isotracer.dereplicator as dereplicator
import isotracer.utils as utils
from isotracer.config import CONFIG


def validate_input(source_dir, conditions):
    """Make sure required input directories exist

    Args:
        source_dir (str or Path): input directory
    """
    source_dir = Path(source_dir)
    assert source_dir.is_dir()
    cppis_dir = source_dir.joinpath("CPPIS")
    assert cppis_dir.is_dir()
    func_dir = source_dir.joinpath("func001")
    assert func_dir.is_dir()
    print("Validation successful!")


def cppis_masterlist(
    source_dir, conditions, exp_name, n_jobs=-1, config=None,
):
    """
    From a directory takes a folder named 'CPPIS' and munges all the cppis.csv files first by conditions
    and then combines all files to create an mz masterlist of all features detected in experiment and basketed
    to align between conditions and replicates.
    """
    if not config:
        config = CONFIG
    source_dir = Path(source_dir)  # Make sure is Path object
    cppis_dir = source_dir.joinpath("CPPIS")
    print(f"Collecting files from {cppis_dir}")
    cppis_files = utils.getfilenames_cppis(
        cppis_dir, conditions
    )  # df containing paths, filenames and conditions

    print(f"Working on BLANKS")
    blank_dir = source_dir.joinpath("BLANKS")
    # make directory if not exists
    blank_dir.mkdir(parents=True, exist_ok=True)

    # peaks in blanks to subtract later
    blanks = utils.munge_cppis(cppis_files, "BLANK", blank_dir, config=config)

    def munge_condition(cond):
        out_dir = source_dir.joinpath(cond)
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Working on {cond}")
        all_collapsed = utils.munge_cppis(cppis_files, cond, out_dir, config=config)
        return all_collapsed

    # Run pre-processing on conditions in separate processes
    prim_dfs = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(munge_condition)(c) for c in conditions
    )

    all_primary = utils.combine_dfs(prim_dfs)
    print("Substracting blanks")
    utils.blank_subtract(
        blanks, all_primary, config=config
    )  # blank subtraction on full dataset
    print("Grouping all features")
    all_primary.reset_index(inplace=True, drop=True)
    dereplicator.group_features(
        all_primary, exp_name, "Exp_ID", config=config
    )  # final grouping - Exp_ID used for func001 munging
    all_primary.to_csv(source_dir.joinpath(f"{exp_name}_All_features.csv"), index=False)
    return all_primary


def isotope_scraper(
    source_dir,
    conditions,
    exp_name,
    master=None,
    min_scans=5,
    n_jobs=-1,
    restart=False,
    config=None,
):
    """This is the final function which will scrape all the isotope data for all ions in the
    """
    if not config:
        config = CONFIG
    print("Running isotope scraper")
    source_dir = Path(source_dir)
    # Enables passing DF instead of loading CSV
    if not isinstance(master, pd.DataFrame):
        master = utils.get_cppis_masterlist(source_dir)

    # TODO: Add flexibility
    func_files = list(source_dir.joinpath("func001").glob("*.csv"))

    # dir for all output files containing scan data for each cond
    out_dir = source_dir.joinpath("All_scan_data")
    out_dir.mkdir(parents=True, exist_ok=True)

    def run_isoslicer(cond, restart=restart):
        outfile = out_dir.joinpath(f"All_ions_{cond}.csv")
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
        # disable N15 searching and use wider tolerance
        # iso = "15N" if cond in ("GLU", "GLN") else "13C"

        len_uni = len(master.Exp_ID.unique())
        exps = master.groupby("Exp_ID")
        print(f"There are {len_uni} ions to look at for {cond}")

        # Parallelized solution
        def do_slice(i, idx, g):
            # nonlocal seen
            if i % 20 == 0 and i > 0:
                print(f"Working on {i}/{len_uni} for {cond}")

            mz, low_scan, high_scan = utils.calc_exp(g)
            #  function 'iso_slicer' slices relevant isotope data for a given mz and all its isotopomers
            #  within given scan range in given func file
            return (
                idx,
                utils.isotope_slicer(
                    func_df,
                    mz,
                    low_scan,
                    high_scan,
                    min_scans=min_scans,
                    iso="C13",
                    mz_tol=utils.get_mz_tol(config),
                ),
            )

        to_mark = joblib.Parallel(n_jobs=-1, prefer="threads")(
            joblib.delayed(do_slice)(i, idx, g) for i, (idx, g) in enumerate(exps)
        )

        print(f"Finished collecting ions for {cond}")
        res_df = utils.mark_func(func_df, to_mark)

        # write iso_scan_data to a file for that condition
        res_df.to_csv(outfile)

    # Run processing of each condition in separate process
    joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(run_isoslicer)(c) for c in conditions)
    # [run_isoslicer(c) for c in conditions]


def isotope_label_detector(
    source_dir, conditions, exp_name, num_cond=1, master=None, n_jobs=-1
):
    source_dir = Path(source_dir)
    scan_dir = source_dir.joinpath("All_scan_data")
    slope_dir = source_dir.joinpath("All_slope_data")
    slope_dir.mkdir(parents=True, exist_ok=True)
    out_dir = source_dir.joinpath("All_isotope_analysis")
    out_dir.mkdir(parents=True, exist_ok=True)
    # Enables passing DF instead of loading CSV
    if not isinstance(master, pd.DataFrame):
        master = utils.get_cppis_masterlist(source_dir)

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
        utils.run_label_analysis(res_df, cond, out_dir)

    # Run processing of each condition in separate process
    joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(run_label_detector)(c) for c in conditions
    )
    # [run_label_detector(c) for c in conditions]
    sum_df = utils.summarize_labels(out_dir, master, conditions)
    sum_df.to_csv(source_dir.joinpath(f"{exp_name}_data_summary.csv"))
    filtered_df = utils.filter_summary(sum_df, num_cond)
    filtered_df.to_csv(source_dir.joinpath(f"{exp_name}_data_summary_filtered.csv"))
