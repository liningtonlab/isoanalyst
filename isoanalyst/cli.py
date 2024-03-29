import copy
import functools
import sys
from pathlib import Path
from typing import Optional

import click

from isoanalyst import VERSION
import isoanalyst.core as core
from isoanalyst.config import get_config, CONFIG
from isoanalyst.input_spec import InputSpec

## Helper functions
def common_options(f):
    options = [
        # click.option(
        #     "--configfile",
        #     help="Config file input",
        #     type=click.Path(exists=True),
        #     required=False,
        # ),
        click.option(
            "--name",
            "-n",
            help="Experiment name - used as output file directory",
            required=True,
        ),
        click.option(
            "--input_specification",
            "-i",
            help="Input specification filename",
            type=click.Path(exists=True),
            required=True,
        ),
    ]
    return functools.reduce(lambda x, opt: opt(x), options, f)


def gt_zero(ctx, param, value):
    try:
        assert value > 0
        return value
    except AssertionError:
        raise click.BadParameter(f"{param} must be greater than 0.")


def gte_zero(ctx, param, value):
    try:
        assert value >= 0
        return value
    except AssertionError:
        raise click.BadParameter(f"{param} must be greater than or equal to 0.")


############################
##      CLI Definition
############################
@click.group()
@click.version_option(version=VERSION)
def cli():
    """Isoanalyst CLI entrypoint

    Processing steps:

    validate - Step 0 : validate input file structure

    prep - Step 1 : Prepare ground truth list of features including dereplication and optional blank removal

    scrape - Step 2 : Scrape all scan data for each of the feature

    analyze - Step 3 : Analyze all scan data for all of the data
    """
    click.echo("Welcome to IsoAnalyst!")


############################
##      VALIDATE
############################
@cli.command("validate")
@common_options
def run_validate(
    name: str, input_specification: Path, configfile: Optional[Path] = None
):
    """
    Performs some simple checks on your input specification file
    """
    click.echo(f"Running validation for {name} on {input_specification}")
    spec = InputSpec.from_csv(input_specification)
    try:
        core.validate_input(spec)
        click.echo(click.style("Validation successful! ✅", fg="green"))

    except AssertionError as e:
        click.echo(click.style("Validation failed... ❌", fg="red"), err=True)
        sys.exit(1)


############################
##      PREP
############################
@cli.command("prep")
@common_options
@click.option(
    "--blank-remove/--no-blank-remove",
    default=True,
    show_default=True,
    help="Perform blank removal during feature aggregation (or not).",
)
@click.option(
    "--minreps",
    type=int,
    default=CONFIG["minreps"],
    show_default=True,
    callback=gt_zero,
    help="Minimum reps to consider in replication comparison",
)
@click.option(
    "--mztol",
    type=float,
    default=CONFIG["tolerances"]["precmz"][1],
    show_default=True,
    callback=gt_zero,
    help="M/Z tolerance in PPM",
)
@click.option(
    "--rttol",
    type=float,
    default=CONFIG["tolerances"]["rettime"][1],
    show_default=True,
    callback=gt_zero,
    help="rettime tolerance in min",
)
@click.option(
    "--jobs",
    "-j",
    type=int,
    default=-1,
    show_default=True,
    help="Maximum number of parallel processes",
)
def run_prep(
    name: str,
    input_specification: Path,
    jobs: int,
    blank_remove: bool,
    minreps: int,
    mztol: int,
    rttol: int,
    configfile: Optional[Path] = None,
):
    """Prepares the ground truth feature list"""
    click.echo(f"Running prep for {name} using {input_specification}")
    spec = InputSpec.from_csv(input_specification)
    config = get_config(configfile, minreps=minreps, mztol=mztol, rttol=rttol)
    source_dir = Path(name.replace(" ", "_"))
    core.generate_featurelist(
        input_spec=spec,
        source_dir=source_dir,
        exp_name=name,
        n_jobs=jobs,
        config=config,
        blank_remove=blank_remove,
    )


############################
##      SCRAPE
############################
@cli.command("scrape")
@common_options
@click.option(
    "--minscans",
    type=int,
    default=2,
    show_default=True,
    callback=gt_zero,
    help="Minimum number of scans",
)
@click.option(
    "--mztol",
    type=float,
    default=CONFIG["tolerances"]["precmz"][1],
    show_default=True,
    callback=gt_zero,
    help="M/Z tolerance in PPM",
)
@click.option(
    "--minintensity",
    type=int,
    default=0,
    show_default=True,
    callback=gte_zero,
    help="Minimum intensity threshold for data",
)
@click.option(
    "--minrt",
    type=float,
    default=0.8,
    show_default=True,
    callback=gt_zero,
    help="Ignores data before minimum RT (minutes)",
)
@click.option(
    "--scanwindow",
    type=int,
    default=10,
    show_default=True,
    callback=gt_zero,
    help="Number of scans to consider for isotope alignment. This only applies if you data is missing scan ranges.",
)
@click.option(
    "--jobs",
    "-j",
    type=int,
    default=-1,
    show_default=True,
    help="Maximum number of parallel processes",
)
def run_scrape(
    name: str,
    input_specification: Path,
    jobs: int,
    minscans: int,
    scanwindow: int,
    mztol: float,
    minintensity: int,
    minrt: float,
    configfile: Optional[Path] = None,
):
    """Collects relevant scan data for all members of the ground truth feature list"""
    click.echo(f"Running scrape for {name} using {input_specification}")
    spec = InputSpec.from_csv(input_specification)
    config = get_config(configfile, mztol=mztol)
    source_dir = Path(name.replace(" ", "_"))
    core.isotope_scraper(
        input_spec=spec,
        source_dir=source_dir,
        exp_name=name,
        n_jobs=jobs,
        min_scans=minscans,
        scanwindow=scanwindow,
        min_intensity=minintensity,
        min_rt=minrt,
        config=config,
    )


############################
##      ANALYZE
############################
@cli.command("analyze")
@common_options
@click.option(
    "--minconditions",
    type=int,
    default=1,
    show_default=True,
    callback=gt_zero,
    help="Minimum number of conditions to output in filtered output",
)
@click.option(
    "--minscans",
    type=int,
    default=5,
    show_default=True,
    callback=gt_zero,
    help="Minimum number of scans",
)
@click.option(
    "--jobs",
    "-j",
    type=int,
    default=-1,
    show_default=True,
    help="Maximum number of parallel processes",
)
def run_analyze(
    name: str,
    input_specification: Path,
    jobs: int,
    minscans: int,
    minconditions: int,
    configfile: Optional[Path] = None,
):
    """Performs Stable Isotope Labelling detecting and analysis"""
    click.echo(f"Running analysis for {name} using {input_specification}")
    spec = InputSpec.from_csv(input_specification)
    source_dir = Path(name.replace(" ", "_"))
    core.isotope_label_detector(
        input_spec=spec,
        source_dir=source_dir,
        exp_name=name,
        min_scans=minscans,
        num_cond=minconditions,
        n_jobs=jobs,
    )
