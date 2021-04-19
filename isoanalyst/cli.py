import copy
import functools
import sys
from pathlib import Path
from typing import Optional

import click

import isoanalyst.core as core
from isoanalyst.config import get_config, CONFIG
from isoanalyst.input_spec import InputSpec

## Helper functions
def common_options(f):
    options = [
        click.option(
            "--config_file",
            help="Config file input",
            type=click.Path(exists=True),
            required=False,
        ),
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


############################
##      CLI Definition
############################
@click.group()
@click.version_option(version="0.1.0")
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
def run_validate(name: str, input_specification: Path, config_file: Optional[Path]):
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
    help="precmz tolerance in PPM",
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
    config_file: Optional[Path],
    jobs: int,
    blank_remove: bool,
    minreps: int,
    mztol: int,
    rttol: int,
):

    click.echo(f"Running prep for {name} using {input_specification}")
    spec = InputSpec.from_csv(input_specification)
    config = get_config(config_file, minreps=minreps, mztol=mztol, rttol=rttol)
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
    "--minreps",
    type=int,
    default=CONFIG["minreps"],
    show_default=True,
    callback=gt_zero,
    help="Minimum reps to consider in replication comparison",
)
@click.option(
    "--minscans",
    type=int,
    default=2,
    show_default=True,
    callback=gt_zero,
    help="Minumum number of scans",
)
@click.option(
    "--mztol",
    type=float,
    default=CONFIG["tolerances"]["precmz"][1],
    show_default=True,
    callback=gt_zero,
    help="precmz tolerance in PPM",
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
def run_scrape(
    name: str,
    input_specification: Path,
    config_file: Optional[Path],
    jobs: int,
    minscans: int,
    minreps: int,
    mztol: float,
    rttol: float,
):
    click.echo(f"Running scrape for {name} using {input_specification}")
    spec = InputSpec.from_csv(input_specification)
    config = get_config(config_file, minreps=minreps, mztol=mztol, rttol=rttol)
    source_dir = Path(name.replace(" ", "_"))
    core.isotope_scraper(
        input_spec=spec,
        source_dir=source_dir,
        exp_name=name,
        n_jobs=jobs,
        min_scans=minscans,
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
    help="Minumum number of scans",
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
    config_file: Optional[Path],
    jobs: int,
    minscans: int,
    minconditions: int,
):
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
