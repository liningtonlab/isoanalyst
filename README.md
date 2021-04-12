# isoanalyst

# Installation and Dependecies

**Create a conda environment with dependencies**

`conda env create -f environment.yml`

This should install a Python 3.8+ environment with all the necessary packages. 

__Note:__ This program uses the `rtree` library, which requires bindings
for `libspatial`. This has been known to cause problems on Windows installs. `rtree` installed from PyPi (`pip install rtree`) should work out of the box. The `environment.yml` has been configured to do this for you.

Then to install the CLI program run

`python setup.py install`

Alternatively, if any changes to the code (i.e.) development are to be done, install the CLI in development mode (I think this has some performance implications but allows for hot-changes in the code)

`python setup.py develop`

### Running the program

Type `isoanalyst -h` to see details on how the CLI works.

### Suggested file structure

```
ROOT
├── Snakefile
└── data
    ├── CPPIS
    └── func001
```

Also, _BLANK_ data should only be present in the CPPIS directory, or else your analysis will get corrupted.

### Snakemake pipeline

The snakemake pipeline allows for much better reproducibility. I have commented in the file to show some details, but to run the pipeline the steps are as follows:

1. Copy the `workflow/Snakefile` to the directory just below where your data is located: i.e. `ROOT` in the example above
2. Edit the `Snakefile` to contain the required parameters - datadir name, experiment name, conditions, tolerances, etc.
3. Run `snakemake -n` from the `ROOT` directory to make sure everything is ready to run
4. Run analysis by running `snakemake --cores 1` from the `ROOT` directory and wait for analysis to complete (or fail).

If you want to re-run an analysis you can either type `snakemake --cores 1 -f` or remove all the `*.done` files.

### Data Requirements

#### 1. Scan Data

Generic CSV tabular inputs minimally with `["ScanIndex", "RetTime", "Mz", "Intensity"]` column headers.

Can import func0001 CSVs from msExpress.

Can import mzMLs from [msConvert](http://proteowizard.sourceforge.net/tools/msconvert.html).

#### 2. Feature lists

Generic CSV tabular inputs minimally with `["ScanIndex", "RetTime", "Mz", "Intensity"]` column headers.

### Development and Testing

To run the tests, you must first install dev dependencies `conda install pytest`.