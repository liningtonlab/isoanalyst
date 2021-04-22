# IsoAnalyst

IsoAnalyst is a mass spectrometry (MS) metabolomics data analysis program designed to determine
the number of stable isotopically labeled (SIL) tracers incorporated into metabolites across parallel
SIL tracer experiments. Unlabeled control samples are required for each SIL tracer used, and separate 
pre-processed datasets from both labeled and unlabeled samples are required input. IsoAnalyst compares
isotopologue distributions between MS features in the unlabeled and labeled datasets and generates a 
summary file containing the number of heavy isotopes incorporated into every MS feature in each SIL condition.

Full program documentation is available at [liningtonlab.github.io/isoanalyst/](https://liningtonlab.github.io/isoanalyst/).

## Installation

**Create a conda environment with dependencies**

`conda env create -f environment.yml`

This should install a Python 3.8+ environment with all the necessary packages. 

You then need to activate the virtual environment.

`conda activate isoanalyst`

__Note:__ This program uses the `rtree` library, which requires bindings
for `libspatialindex`. This has been known to cause problems on Windows installs. `rtree` installed 
from PyPi (`pip install rtree`) should work out of the box. The `environment.yml` installation specification 
has been configured to do this for you. If you receive an error on running
`isoanalyst --version` this library is the most likely culprit. You can try
to install it manually again with `conda install -c conda-forge rtree` (this has been found to be 
required on macOS).

Then to install the IsoAnalyst CLI program run

`python setup.py install`

To test whether the program has been properly installed, you can run

`isoanalyst --version`


## Running the program

The `isoanalyst` program is split into four sub-applications which are intended to be run in sequence.
To run all steps in a simple, reproducible manner, the [Snakemake pipeline](#snakemake-pipeline) is recommended.

Two common parameters are required for all steps of the program:

`-n` or `--name` is the experiment name and is used as the output path and for marking mass features.

`-i` or `--input_specification` is the path the [input specification file](#input-specification).


#### Steps

1. Validation - `isoanalyst validate -n EXPNAME -i input_spec.csv` - validates your file paths, and input specification.
2. Prep - `isoanalyst prep -n EXPNAME -i input_spec.csv` - generated the ground truth feature list for isotope label detection.
3. Scrape - `isoanalyst scrape -n EXPNAME -i input_spec.csv` - collects all scan data relevant to the ground truth feature list.
4. Analyze - `isoanalyst analyze -n EXPNAME -i input_spec.csv` - performs SIL analysis and generates summary reports.

Type `isoanalyst --help` to see complete details on how the CLI works. (On Windows, replace `isoanalyst` with `isoanalyst.exe`.)

The CLI options are explained in greater detail [in the documentation @ liningtonlab.github.io/isoanalyst/](https://liningtonlab.github.io/isoanalyst/#cli-options)


### Input Specification

To remove the error prone inference from the program, a simple CSV input specification has been devised for input files and related parameters.
[Below is a simple example](#example-input-specification).

__All columns are required__. The columns are defined in the following manner

- `filepath` 
  - This is the path the input file, either complete or relative to where the `isoanalyst` program is run from.
- `organism`
  - The name of the organism for that particular experiment. This may commonly be all the same value.
- `type`
  - The type of input file: One of feature list __`f`__, or all scan __`s`__.
- `element`
  - The element involved in SIL detected (currently the program only support C and N) in the given condition.
- `isotope`
  - The isotope being detected in the given condition.
- `condition`
  - The name of the experimental SIL condition (or "BLANK").
- `replicate`
  - The replicate number for a given condition.

All fields are required for all conditions __except__ `isotope` and `element` for the "BLANK" condition.

Full documentation on input specification is available [in the documentation @ liningtonlab.github.io/isoanalyst/](https://liningtonlab.github.io/isoanalyst/#input-specification).


#### Example Input Specification

| filepath                                                                          | organism | type | element | isotope | condition | replicate |
| --------------------------------------------------------------------------------- | -------- | ---- | ------- | ------- | --------- | --------- |
| feature_lists/20180409_RLUS135312ACED0-1_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | ACE       | 1         |
| feature_lists/20180409_RLUS135312ACED0-2_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | ACE       | 2         |
| feature_lists/20180409_RLUS135312ACED0-3_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | ACE       | 3         |
| feature_lists/20180409_RLUS135312ACED0-4_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | ACE       | 4         |
| feature_lists/20180409_RLUS135314GLUD0-1_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | GLU       | 1         |
| feature_lists/20180409_RLUS135314GLUD0-2_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | GLU       | 2         |
| feature_lists/20180409_RLUS135314GLUD0-3_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | GLU       | 3         |
| feature_lists/20180409_RLUS135314GLUD0-4_seen.mzml_chromatograms_deconvoluted.csv | RLUS1353 | f    | C       | 12      | GLU       | 4         |
| mzmls/20180409_RLUS135312ACED0-1_seen.mzml                                        | RLUS1353 | s    | C       | 12      | ACE       | 1         |
| mzmls/20180409_RLUS135312ACED0-2_seen.mzml                                        | RLUS1353 | s    | C       | 12      | ACE       | 2         |
| mzmls/20180409_RLUS135312ACED0-3_seen.mzml                                        | RLUS1353 | s    | C       | 12      | ACE       | 3         |
| mzmls/20180409_RLUS135312ACED0-4_seen.mzml                                        | RLUS1353 | s    | C       | 12      | ACE       | 4         |
| mzmls/20180409_RLUS135314GLUD0-1_seen.mzml                                        | RLUS1353 | s    | C       | 12      | GLU       | 1         |
| mzmls/20180409_RLUS135314GLUD0-2_seen.mzml                                        | RLUS1353 | s    | C       | 12      | GLU       | 2         |
| mzmls/20180409_RLUS135314GLUD0-3_seen.mzml                                        | RLUS1353 | s    | C       | 12      | GLU       | 3         |
| mzmls/20180409_RLUS135314GLUD0-4_seen.mzml                                        | RLUS1353 | s    | C       | 12      | GLU       | 4         |
| mzmls/20180410_RLUS135313ACED0-1_seen.mzml                                        | RLUS1353 | s    | C       | 13      | ACE       | 1         |
| mzmls/20180410_RLUS135313ACED0-2_seen.mzml                                        | RLUS1353 | s    | C       | 13      | ACE       | 2         |
| mzmls/20180410_RLUS135313ACED0-3_seen.mzml                                        | RLUS1353 | s    | C       | 13      | ACE       | 3         |
| mzmls/20180410_RLUS135313ACED0-4_seen.mzml                                        | RLUS1353 | s    | C       | 13      | ACE       | 4         |
| mzmls/20180410_RLUS135315GLUD0-1_seen.mzml                                        | RLUS1353 | s    | N       | 15      | GLU       | 1         |
| mzmls/20180410_RLUS135315GLUD0-2_seen.mzml                                        | RLUS1353 | s    | N       | 15      | GLU       | 2         |
| mzmls/20180410_RLUS135315GLUD0-3_seen.mzml                                        | RLUS1353 | s    | N       | 15      | GLU       | 3         |
| mzmls/20180410_RLUS135315GLUD0-4_seen.mzml                                        | RLUS1353 | s    | N       | 15      | GLU       | 4         |


### Snakemake pipeline

To run the full IsoAnalyst pipeline, [snakemake](https://snakemake.readthedocs.io/en/stable/) is recommended for much better reproducibility. 
An example workflow including a `Snakefile` is available in the `example_workflow` directory of this repository.

Running snakemake is as simple as navigating to the root directory for your data (`cd example_workflow` in this example),
and running `snakemake -j1`. This will produce several `.done` file for each step of the pipeline.
If you change tolerances or configuration and would like to re-run the pipeline, either remove these `.done` files (`rm *.done`)
OR use `snakemake -j1 -F` to force snakemake to re-run all steps.

If you would like to use snakemake in your own work, following these steps:

1. Copy the `example_workflow/Snakefile` to the directory just below where you would like to run the program from: i.e. `example_workflow` in the example above.
2. Edit the `Snakefile` to contain the required parameters - experiment name, input specification, tolerances, etc.
3. Navigate to the directory with your Snakefile and run `snakemake -n` to make sure everything is ready to run
4. Run analysis by running `snakemake -j1`


## Data Requirements

See the [complete documentation @ liningtonlab.github.io/isoanalyst/](https://liningtonlab.github.io/isoanalyst/#data-requirements) for details.
