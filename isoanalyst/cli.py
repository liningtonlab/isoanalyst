import copy
import argparse
import logging
from pathlib import Path

import isoanalyst.core as core
from isoanalyst.config import CONFIG


def run_validate(args):
    print("Running validation...")
    try:
        core.validate_input(args.input_spec)
    except AssertionError as e:
        logger.error(e)


def run_prep(args):
    print("Running prep...")
    core.cppis_masterlist(
        source_dir=args.source_dir,
        conditions=args.conditions,
        exp_name=args.name,
        n_jobs=args.jobs,
        config=args.config,
    )


def run_scrape(args):
    print("Running scrape...")
    core.isotope_scraper(
        source_dir=args.source_dir,
        conditions=args.conditions,
        exp_name=args.name,
        n_jobs=args.jobs,
        restart=-args.retry,
        min_scans=args.minscans,
        config=args.config,
    )


def run_analyze(args):
    print("Running analyze...")
    core.isotope_label_detector(
        source_dir=args.source_dir,
        conditions=args.conditions,
        exp_name=args.name,
        num_cond=args.minconditions,
        min_scans=args.minscans,
        n_jobs=args.jobs,
    )


# Define Parser and Options
# Main parser
parser = argparse.ArgumentParser(
    prog="isoanalyst",
    description="Analyze MS data for isotopic labelling experiments.",
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "step",
    help="""Processing step:

    validate - Step 0 : validate input file structure
    prep - Step 1 : Prepare master list of ions including dereplication and optional blank removal
    scrape - Step 2 : Scrape all scan data for each of the ions
    analyze - Step 3 : Analyze all scan data for all of the data
    """,
    choices=("validate", "prep", "scrape", "analyze"),
)

parser.add_argument(
    "-i",
    "--input_specification"
    help="CSV file containing input specifications",
    type=Path,
)

parser.add_argument(
    "-j",
    "--jobs",
    type=int,
    default=-1,
    help="Number of jobs to run in parallel",
)

parser.add_argument(
    "--minscans",
    type=int,
    help="ONLY FOR SCRAPE STEP: Minumum number of scans for a real isotopomer (Default = 5)",
)

parser.add_argument(
    "--print_config",
    action="store_true",
    help="Print the configuration",
)

parser.add_argument(
    "--colstomatch",
    nargs="+",
    help="Column names to match in dereplication",
)

parser.add_argument("--mztol", type=float, help="PrecMz tolerance in PPM")

parser.add_argument("--rttol", type=float, help="RetTime tolerance in min")

parser.add_argument(
    "--minreps", type=int, help="Minium reps to consider in replication comparison"
)

parser.add_argument(
    "--minconditions",
    type=int,
    help="ONLY FOR ANALYZE STEP: Minimum number of conditions to output in filtered output",
    default=1,
)


def main():
    args = parser.parse_args()
    # Parse extra config options
    args.config = copy.deepcopy(CONFIG)
    if args.colstomatch:
        assert len(args.colstomatch) > 1
        args.config["ColsToMatch"] = args.colstomatch
    if args.mztol:
        assert args.mztol > 0
        args.config["Tolerances"]["PrecMz"][1] = args.mztol
    if args.rttol:
        assert args.rttol > 0
        args.config["Tolerances"]["RetTime"][1] = args.rttol
    if args.minreps:
        assert args.minreps > 0
        args.config["MinReps"] = args.minreps
    if args.print_config:
        print(args.config)
    # Default to the name of the source_dir
    if not args.name:
        args.name = args.source_dir.name
    # Append secondary named to primary
    if args.secondary:
        args.conditions = [f"{c}{s}" for c in args.conditions for s in args.secondary]
        # print("Detected secondary conditions")
        # print(args.conditions)
    if args.retry and not args.step == "scrape":
        logging.warning("'--retry' flag has no effect")
    if args.minscans and not args.step in ["scrape", "analyze"]:
        logging.warning("'--minscans' flag has no effect")
    if args.step in ["scrape", "analyze"] and not args.minscans:
        args.minscans = 5
    if args.step == "validate":
        run_validate(args)
    if args.step == "prep":
        run_prep(args)
    if args.step == "scrape":
        run_scrape(args)
    if args.step == "analyze":
        run_analyze(args)


if __name__ == "__main__":
    main()
