## Welcome to IsoAnalyst

IsoAnalyst is a mass spectrometry (MS) metabolomics data analysis program designed to determine
the number of stable isotopically labeled (SIL) tracers incorporated into metabolites across parallel
SIL tracer experiments. Unlabeled control samples are required for each SIL tracer used, and separate 
pre-processed datasets from both labeled and unlabeled samples are required input. IsoAnalyst compares
isotopologue distributions between MS features in the unlabeled and labeled datasets and generates a 
summary file containing the number of heavy isotopes incorporated into every MS feature in each SIL condition.


## Installation

**Create a conda environment with dependencies**

`conda env create -f environment.yml`

This should install a Python 3.8+ environment with all the necessary packages. 

You then need to activate the virtual environment.

`conda activate isoanalyst`

Then to install the IsoAnalyst CLI program run

`python setup.py install`

To test whether the program has been properly installed, you can run

`isoanalyst --version`

__Note:__ This program uses the `rtree` library, which requires bindings
for `libspatialindex`. This has been known to cause problems on Windows installs. `rtree` installed 
from PyPi (`pip install rtree`) should work out of the box. The `environment.yml` installation specification 
has been configured to do this for you. If you receive an error on running
`isoanalyst --version` this library is the most likely culprit. You can try
to install it manually again with `conda install -c conda-forge rtree` (this has been found to be 
required on macOS).


## CLI Options

The `isoanalyst` program is split into four sub-applications which are intended to be run in sequence.
To run all steps in a simple, reproducible manner, the [Snakemake pipeline](https://github.com/liningtonlab/isoanalyst#snakemake-pipeline) 
is recommended.

Both the root app and each sub-app have their own CLI help context, which is accessible through the `--help` flag.

##### Root app help

```
[jvansan@cpu ~]$ isoanalyst --help
Usage: isoanalyst [OPTIONS] COMMAND [ARGS]...

  Isoanalyst CLI entrypoint

  Processing steps:

  validate - Step 0 : validate input file structure

  prep - Step 1 : Prepare ground truth list of features including
  dereplication and optional blank removal

  scrape - Step 2 : Scrape all scan data for each of the feature

  analyze - Step 3 : Analyze all scan data for all of the data

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  analyze   Performs Stable Isotope Labelling detecting and analysis
  prep      Prepares the ground truth feature list
  scrape    Collects relevant scan data for all members of the ground truth...
  validate  Performs some simple checks on your input specification file
```

Two common parameters are __required__ for all steps of the program:

`-n` or `--name` is the experiment name and is used as the output path and for marking mass features.

`-i` or `--input_specification` is the path the [input specification file](#input-specification).

The number of parallel jobs is available for all steps except the validate step.

`-j` or `--jobs`

The number of parallel jobs to run. A job for each can condition will be run in parallel.
Defaults to -1, indicating all available detected CPU cores/threads.

#### Validate

The `isoanalyst validate` step performs some simple checks on your input specification file, such as making
sure all the required information is present, and that there are no missing specified input files.

##### Validate step help

```
[jvansan@cpu ~]$ isoanalyst validate --help
Welcome to IsoAnalyst!
Usage: isoanalyst validate [OPTIONS]

  Performs some simple checks on your input specification file

Options:
  -i, --input_specification PATH  Input specification filename  [required]
  -n, --name TEXT                 Experiment name - used as output file
                                  directory  [required]

  --help                          Show this message and exit.
```

There are no additional options for the validate step.


#### Prep

The `isoanalyst prep` step collects all your feature lists to produce a single ground truth feature lists of
unlabelled features.

This step performs a replicate comparison on each condition, and can remove features present from a "BLANK"
condition. All the features from the conditions are then aggregated into the ground truth feature list, which 
will be output as `NAME/NAME_all_features.csv` (where `NAME` is the specified experiment name).

Individual lists of ions from each condition (`COND`) are all stored in  `NAME/COND/all_ions_sorted.csv` for all the
listed features and `NAME/COND/all_ions_averaged.csv` for all the replicate compared features.

##### Prep step help

```
[jvansan@cpu ~]$ isoanalyst prep --help
Welcome to IsoAnalyst!
Usage: isoanalyst prep [OPTIONS]

  Prepares the ground truth feature list

Options:
  -i, --input_specification PATH  Input specification filename  [required]
  -n, --name TEXT                 Experiment name - used as output file
                                  directory  [required]

  --blank-remove / --no-blank-remove
                                  Perform blank removal during feature
                                  aggregation (or not).  [default: True]

  --minreps INTEGER               Minimum reps to consider in replication
                                  comparison  [default: 3]

  --mztol FLOAT                   M/Z tolerance in PPM  [default: 10.0]
  --rttol FLOAT                   rettime tolerance in min  [default: 0.03]
  -j, --jobs INTEGER              Maximum number of parallel processes
                                  [default: -1]

  --help                          Show this message and exit.
```

##### OPTIONS

- `--blank-remove / --no-blank-remove`

Perform blank removal on the ground truth feature list. This will only work if
any "BLANK" condition input feature files are present in the input specification.
`True` == `--blank-remove` by default as blanks are highly recommended.

- `--minreps`

The minimum number of replicates a feature must be present in a single condition
to be considered a real feature. Default value is 3.

- `--mztol`

The M/Z tolerance for replicate comparison in PPM. Defaults to 10.0.

- `--rttol`

The retention time tolerance for replicate comparison in minutes. Defaults to 0.03.


#### Scrape

The `isoanalyst scrape` step collects all the scan by scan data for each of the features in the
ground truth feature list.

This step performs collects all scans for each of the ground truth features from the prepared
labelled and unlabelled all scan data. It collects all the scans for each replicate individually
to make the analysis step possible.

Paired down all scan data for each condition after import and optional intensity tolerance filtering (`COND`) 
are stored in `NAME/all_scan_data/all_scans_COND.csv`,and extracted scans are stored in `NAME/all_scan_data/all_ions_COND.csv`.

##### Scrape step help

```
[jvansan@cpu ~]$ isoanalyst scrape --help
Welcome to IsoAnalyst!
Usage: isoanalyst scrape [OPTIONS]

  Collects relevant scan data for all members of the ground truth feature
  list

Options:
  -i, --input_specification PATH  Input specification filename  [required]
  -n, --name TEXT                 Experiment name - used as output file
                                  directory  [required]

  --minscans INTEGER              Minimum number of scans  [default: 2]
  --mztol FLOAT                   M/Z tolerance in PPM  [default: 10.0]
  --minintensity INTEGER          Minimum intensity threshold for data
                                  [default: 0]

  --minrt FLOAT                   Ignores data before minimum RT (minutes)
                                  [default: 0.8]

  --scanwindow INTEGER            Number of scans to consider for isotope
                                  alignment. This only applies if you data is
                                  missing scan ranges.  [default: 10]

  -j, --jobs INTEGER              Maximum number of parallel processes
                                  [default: -1]

  --help                          Show this message and exit.
```

##### OPTIONS

- `--minintensity`

The minimum intensity threshold for scan data. This is applied during scan data import.
While you may be able to detect more labelled peaks with a lower intensity, this can
drastically slow down processing. Defaults to 0, i.e. no filtering.

- `--minrt`

The minimum retention time value for considering scan data. In our analysis we ignored the
first 0.8 minutes, thus this is the default value.

- `--minscans`

The minimum number of scans to be considered for collection. Less than two scans will
does not work during analysis, thus the Default is 2.

- `--mztol`

The M/Z tolerance for scan scrapping in PPM. Defaults to 10.0.

- `--scanwindow`

The number of scans to consider on either side of the central scan of a ground truth feature.
This is only applicable if you import data without a scan range for features, such as when
importing from MzMine feature lists.

#### Analyze

The `isoanalyst analyze` step uses the ground truth feature list and collected scans to
detect stable isotope labelling and summarize the extent of labelling of each feature in
each of the specified conditions.

For each condition (`COND`), the algorithm computes slope data between each of the isotopomers, output
in `NAME/all_slope_data/all_slope_data_COND.csv`. These slopes are then used to determine the extend of
labelling, output in `NAME/all_isotope_analysis/iso_analysis_COND.csv`. These results are aggregated in
`NAME/NAME_data_summary.csv` and a filtered version (based on the `--minconditions` flag) in
`NAME/NAME_data_summary_filtered.csv`.

##### Analyze step help

```
[jvansan@cpu ~]$ isoanalyst analyze --help
Welcome to IsoAnalyst!
Usage: isoanalyst analyze [OPTIONS]

  Performs Stable Isotope Labelling detecting and analysis

Options:
  -i, --input_specification PATH  Input specification filename  [required]
  -n, --name TEXT                 Experiment name - used as output file
                                  directory  [required]

  --minconditions INTEGER         Minimum number of conditions to output in
                                  filtered output  [default: 1]

  --minscans INTEGER              Minimum number of scans  [default: 5]
  -j, --jobs INTEGER              Maximum number of parallel processes
                                  [default: -1]

  --help                          Show this message and exit.
```

##### OPTIONS

- `--minconditions`

The minumum number of conditions to be considered for annotation in the final filtered output file.

- `--minscans`

The minimum number of scans to be considered for slope analysis used in isotope labelling detection.
This number must be greater than or equal to the number of `minscans` from the scrape step.


## Input Specification

To remove the error prone inference from the program, a simple CSV input specification has been devised for input files and related parameters.
[Below is a simple example](#example-input-specification) which matches the [`example_workflow` in the repository](https://github.com/liningtonlab/isoanalyst/tree/main/example_workflow).

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

##### Example Input Specification

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

##### Notes

Although this is discussed in greater details in the [Data Requirements Section](#data-requirements), here are some important
guidelines for organizing your input specification.

1. Blank data should only be present as a feature list
2. Include only unlabelled data as features lists. Doing otherwise is likely to corrupt your analysis.
3. Detection of the `M0` peak in the unlabelled condition should always be done with the 12/13 C isotope. This is why in the `example_workflow` input specification, the isotope = 12 and element = C for the GLU unlabelled condition, even though the labelling is done with 15N (isotope = 15, element = N for the GLU labelled condition).


## Data Requirements

#### 1.Scan Data

Files containing centroided scan data are required for all unlabeled controls and labeled samples. Four technical replicates are recommended for each condition.

Generic CSV tabular inputs minimally with `["FunctionScanIndex", "RT", "MZ", "Intensity"]` column headers. These are compatible msExpress func001 files.

Scan data can be imported as mzMLs from [msConvert](http://proteowizard.sourceforge.net/tools/msconvert.html). Standard GNPS settings recommended for MSconvert. 


#### 2. Feature Lists

Feature lists are required for unlabeled control samples only. Four technical replicates are recommended for every feedstock condition used. Three additional feature lists are required for solvent injection blanks. 

Generic CSV tabular inputs minimally with `["Sample", "PrecMz", "PrecZ", "PrecIntensity", "RetTime", "ScanLowRange", "ScanHighRange"]` column headers. These are compatible msExpress CPPIS files.

Feature lists may be made using the standard MZmine 2 workflow. Parameters for mass detection, chromatogram building, and deconvolution are dependent on the MS instrument used and data quality (ie signal to noise). The isotope peak grouper function must be used in order to deistope all features in the list. Every feature in this list should be identified by the monosisotopic m/z value. 

Export .csv files from MZmine 2 with the columns `["row m/z", "row retention time"]`

Peak intensity is not required for feature lists; m/z and retention time values are used to get the isotopologue peak intensities from the scan data. The program will then set a window of 10 scans on both sides of the center scan index detected during all scan importing.


#### Data Collection Recommendation

Not supported for DDA.

Either DIA with interleaved scans or MS1 only.


### Development

Unit tests run on the [`pytest` framework](https://docs.pytest.org/en/6.2.x/).
Be sure to install pytest in your development environment (`pip instll pytest` or `conda install pytest`) 
and make sure all tests pass before committing code.

Running pytest is as simple as running `pytest` from the root directory of the code repository.

To install the `isoanalyst` toolchain in a development environment, during installation use

`pytest setup.py develop` 
