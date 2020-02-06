import argparse
import logging
from pathlib import Path

import isotracer.core as core

def run_validate(args):
    print("Running validation...")
    try:
        core.validate_input(args.source_dir, args.conditions)
    except AssertionError:
        print(f"source_dir: {args.source_dir} does not have the required structure")
        print("""
Please ensure directory format matches: \033[93m
    source_dir 
    ├── CPPIS 
    └── func001
    \033[0m """)

def run_prep(args):
    print("Running prep...")
    core.cppis_masterlist(
        source_dir=args.source_dir,
        conditions=args.conditions,
        exp_name=args.name,
        n_jobs=args.jobs,
    )

def run_scrape(args):
    print("Running scrape...")
    core.isotope_scraper(
        source_dir=args.source_dir,
        conditions=args.conditions,
        exp_name=args.name,
        n_jobs=args.jobs,
        restart=-args.retry,
    )

def run_analyze(args):
    print("Running analyze...")
    core.isotope_label_detector(
        source_dir=args.source_dir,
        conditions=args.conditions,
        n_jobs=args.jobs,
    )

# Define Parser and Options
# Main parser
parser = argparse.ArgumentParser(
    prog = "IsoTracer",
    description = "Analyze MS data for isotopic labelling experiments.",
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument(
    "step",
    help="""Processing step:
    
    validate - Step 0 : validate input file structure
    prep - Step 1 : Prepare master list of ions including dereplication and optional blank removal
    scrape - Step 2 : Scrape all scan data for each of the ions
    analyze - Step 3 : Analyze all scan data for all of the data
    """,
    choices=("validate", "prep", "scrape", "analyze")
)

parser.add_argument(
    "source_dir",
    help="Directory containing input data",
    type=Path,
)

parser.add_argument(
    "conditions",
    help="Conditions to consider (Minimum 1)",
    nargs="+",
)

parser.add_argument(
    "-n",
    "--name",
    help="Experiment name",
    default="IsoTracer",
)

parser.add_argument(
    "-s",
    "--secondary",
    help="Secondary conditions to consider (Minimum 1 if specified)",
    nargs="+",
)

parser.add_argument(
    "-j", 
    "--jobs", 
    type=int,
    default=-1,
    help="Number of jobs to run in parallel",
)

parser.add_argument(
    "--retry",
    action="store_true",
)


def main():
    args = parser.parse_args()
    if args.retry and not args.step == "scrape":
        logging.warn("'--retry' flag has no effect")
    if args.step == "validate":
        run_validate(args)
    if args.step == "prep":
        run_prep (args)
    if args.step == "scrape":
        run_scrape(args)
    if args.step == "analyze":
        run_analyze(args)


if __name__ == "__main__":
    main()
