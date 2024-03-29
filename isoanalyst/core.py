#!/usr/bin/env python
# coding: utf-8
from pathlib import Path
from typing import Dict, Optional

import joblib
import pandas as pd

import isoanalyst.dereplicator as dereplicator
import isoanalyst.utils as utils
from isoanalyst.config import CONFIG, get_mz_tol
from isoanalyst.input_spec import InputSpec


def validate_input(input_spec: InputSpec):
    """Make sure required input directories exist

    Args:
        input_spec (InputSpec): input file specification
    """
    assert input_spec.validate()


def generate_featurelist(
    input_spec: InputSpec,
    source_dir: Path,
    exp_name: str,
    config: Dict,
    n_jobs: int = -1,
    blank_remove: bool = True,
):
    """
    Create an mz ground truth list of all features detected in experiment and basketed
    to align between conditions and replicates. Optional blank subtraction.
    """
    print("Collecting feature list files")
    # define here so input_spec and config in scope and not needed as params
    def do_munge_featurelist(cond: str):
        out_dir = source_dir.joinpath(cond)
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Working on {cond}")
        all_collapsed = utils.munge_featurelist(
            inp_spec=input_spec, cond=cond, out_dir=out_dir, config=config
        )
        return all_collapsed

    # peaks in blanks to subtract later
    # This should not break in the event that there are no blanks
    # but the printing may be misleading
    if input_spec.get_feature_filepaths("blank"):
        blanks = do_munge_featurelist("blank")
    else:
        print("No blanks found")
        blank_remove = False

    # Run pre-processing on conditions in separate processes
    conditions = input_spec.get_conditions()
    cond_dfs = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(do_munge_featurelist)(c) for c in conditions
    )

    all_cond_df = utils.combine_dfs(cond_dfs)
    if blank_remove:
        print("Substracting blanks")
        utils.blank_subtract(blanks, all_cond_df, config=config)
    else:
        print("Not subtracting blanks")
    print("Grouping all features")
    all_cond_df.reset_index(inplace=True, drop=True)
    dereplicator.group_features(
        all_cond_df, exp_name, "exp_id", config=config
    )  # final grouping - exp_id used for scan munging
    all_cond_df.to_csv(source_dir.joinpath(f"{exp_name}_all_features.csv"), index=False)
    return all_cond_df


def isotope_scraper(
    input_spec: InputSpec,
    source_dir: Path,
    exp_name: str,
    config: Dict,
    min_scans: int,
    scanwindow: int,
    min_intensity: int,
    min_rt: float,
    n_jobs: int,
):
    """Scrape all the isotope data for all ions in the all scan data"""
    conditions = input_spec.get_isoconditions()
    print("Running isotope scraper")
    features = utils.get_featurelist(source_dir=source_dir, exp_name=exp_name)

    # dir for all output files containing scan data for each cond
    out_dir = source_dir.joinpath("all_scan_data")
    out_dir.mkdir(parents=True, exist_ok=True)

    def run_isoslicer(cond, features=features):
        outfile = out_dir.joinpath(f"all_ions_{cond}.csv")

        # add replicate scan files for current condition
        c_scans = input_spec.get_scan_filepaths(cond)
        scanfile = out_dir.joinpath(f"all_scans_{cond}.csv")
        if scanfile.exists():
            print(f"Loading all scan file - {scanfile}")
            scan_df = pd.read_csv(scanfile)
        else:
            c_dfs = (
                utils.prep_scan(
                    Path(f),
                    input_spec,
                    min_intensity=min_intensity,
                    min_rt=min_rt,
                )
                for f in c_scans
            )
            # merged data frame of all scans for the given condition
            print(f"Combining files for {cond}")
            scan_df = utils.combine_dfs(c_dfs)
            # Save checkpoint DF
            print(f"Saving all scan file - {scanfile}")
            scan_df.to_csv(scanfile, index=False)

        # If no scan ranges in feature_df, we need to add them
        if not "lowscan" in features.columns:
            # Make a copy for parallel safety
            print(f"Assigning scan range window to feature list for {cond}")
            features = utils.add_scan_window(features, scan_df, scanwindow=scanwindow)

        len_uni = len(features.exp_id.unique())
        exps = features.groupby("exp_id")
        print(f"There are {len_uni} ions to look at for {cond}")

        # Parallelized solution
        def do_slice(i, idx, g):
            # nonlocal seen
            if i % 20 == 0 and i > 0:
                print(f"Working on {i}/{len_uni} for {cond}")

            mz, low_scan, high_scan = utils.calc_exp(g)
            #  function 'iso_slicer' slices relevant isotope data for a given mz and all its isotopomers
            #  within given scan range in given scan file
            return (
                idx,
                utils.isotope_slicer(
                    scan_df,
                    mz,
                    low_scan,
                    high_scan,
                    min_scans=min_scans,
                    mz_tol=get_mz_tol(config),
                ),
            )

        to_mark = joblib.Parallel(n_jobs=-1, prefer="threads")(
            joblib.delayed(do_slice)(i, idx, g) for i, (idx, g) in enumerate(exps)
        )

        print(f"Finished collecting ions for {cond}")
        res_df = utils.mark_scans(scan_df, to_mark)

        # write iso_scan_data to a file for that condition
        res_df.to_csv(outfile)

    # Run processing of each condition in separate process
    joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(run_isoslicer)(c) for c in conditions)
    # [run_isoslicer(c) for c in conditions]


def isotope_label_detector(
    input_spec: InputSpec,
    source_dir: Path,
    exp_name: str,
    min_scans: int = 5,
    num_cond: int = 1,
    n_jobs: int = -1,
):
    scan_dir = source_dir.joinpath("all_scan_data")
    slope_dir = source_dir.joinpath("all_slope_data")
    slope_dir.mkdir(parents=True, exist_ok=True)
    out_dir = source_dir.joinpath("all_isotope_analysis")
    out_dir.mkdir(parents=True, exist_ok=True)
    features = utils.get_featurelist(source_dir=source_dir, exp_name=exp_name)
    conditions = input_spec.get_conditions()

    def run_label_detector(cond):
        out_file = slope_dir.joinpath(f"all_slope_data_{cond}.csv")
        print(f"Detecting labels in {cond}")
        
        #load both unlabeled and labeled scan data for cond
        scan_files = list(scan_dir.glob(f"all_ions_*{cond}.csv"))
        scan_dfs = []
        for s in scan_files:
            assert s.exists()
            s_df = pd.read_csv(s)
            scan_dfs.append(s_df)
        
        df = utils.combine_dfs(scan_dfs)
        
        #fil = scan_dir.joinpath(f"all_ions_{cond}.csv")
        #assert fil.exists()
        #df = pd.read_csv(fil)
        
        grouped = df.groupby(["exp_id", "isotope", "condition"])
        data = []
        for (e_id, iso, c), g in grouped:
            res = utils.calc_rep_stats(g, e_id, iso, c, min_scans=min_scans)
            # print(data)
            if len(res) < 1:
                continue
            data.extend(res)
        if data:
            agg_df = pd.DataFrame(data)
            res_df = utils.aggregate_results(agg_df)
            res_df.to_csv(out_file, index=False)
            utils.run_label_analysis(res_df, cond, out_dir)
        else:
            print(f"No labels to detect for {cond}")

    # Run processing of each condition in separate process
    joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(run_label_detector)(c) for c in conditions
    )
    # [run_label_detector(c) for c in conditions]
    sum_df = utils.summarize_labels(out_dir, features, conditions)
    sum_df.to_csv(source_dir.joinpath(f"{exp_name}_data_summary.csv"))
    filtered_df = utils.filter_summary(sum_df, num_cond)
    filtered_df.to_csv(source_dir.joinpath(f"{exp_name}_data_summary_filtered.csv"))
