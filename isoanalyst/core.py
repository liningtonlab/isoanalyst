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
    blanks = do_munge_featurelist("blank")

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
    )  # final grouping - exp_id used for func001 munging
    all_cond_df.to_csv(source_dir.joinpath(f"{exp_name}_all_features.csv"), index=False)
    return all_cond_df


def isotope_scraper(
    input_spec: InputSpec,
    source_dir: Path,
    exp_name: str,
    config: Dict,
    min_scans: int,
    min_intensity: int,
    min_rt: float,
    n_jobs: int,
    restart: bool = True,  # retry by default
):
    """Scrape all the isotope data for all ions in the all scan data"""
    conditions = input_spec.get_conditions()
    print("Running isotope scraper")
    features = utils.get_featurelist(source_dir=source_dir, exp_name=exp_name)

    # dir for all output files containing scan data for each cond
    out_dir = source_dir.joinpath("all_scan_data")
    out_dir.mkdir(parents=True, exist_ok=True)

    def run_isoslicer(cond, restart=restart):
        outfile = out_dir.joinpath(f"all_ions_{cond}.csv")
        if not restart and outfile.exists():
            print(f"{cond} already processed")
            return

        # add replicate func files for current condition
        c_funcs = input_spec.get_scan_filepaths(cond)
        c_dfs = (
            utils.prep_scan(
                Path(f), input_spec, min_intensity=min_intensity, min_rt=min_rt
            )
            for f in c_funcs
        )

        # merged data frame of all func001s for the given condition 'c'
        print(f"Combining files for {cond}")
        func_df = utils.combine_dfs(c_dfs)

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
            #  within given scan range in given func file
            return (
                idx,
                utils.isotope_slicer(
                    func_df,
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
        res_df = utils.mark_func(func_df, to_mark)

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
        fil = scan_dir.joinpath(f"all_ions_{cond}.csv")
        assert fil.exists()
        df = pd.read_csv(fil)
        grouped = df.groupby(["exp_id", "isotope", "condition"])
        data = []
        for (e_id, iso, c), g in grouped:
            res = utils.calc_rep_stats(g, e_id, iso, c, min_scans=min_scans)
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
    sum_df = utils.summarize_labels(out_dir, features, conditions)
    sum_df.to_csv(source_dir.joinpath(f"{exp_name}_data_summary.csv"))
    filtered_df = utils.filter_summary(sum_df, num_cond)
    filtered_df.to_csv(source_dir.joinpath(f"{exp_name}_data_summary_filtered.csv"))
