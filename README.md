# isoanalyst

# Installation and Dependecies

**Create a conda environment with dependencies**

`conda create -n isoanalyst -c bioconda -c conda-forge snakemake-minimal rtree scipy joblib pandas`

This should install a Python 3.8 environment with all the necessary packages.

Then to install the CLI program run

`python setup.py install`

Alternatively, if any changes to the code (i.e.) development are to be done, install the CLI in development mode (I think this has some performance implications but allows for hot-changes in the code)

`python setup.py develop`

#### Running the program

Type `isoanalyst -h` to see details on how the CLI works.

#### Suggested file structure

```
ROOT
├── Snakefile
└── data
    ├── CPPIS
    └── func001
```

Also, _BLANK_ data should only be present in the CPPIS directory, or else your analysis will get corrupted.

#### Snakemake pipeline

The snakemake pipeline allows for much better reproducibility. I have commented in the file to show some details, but to run the pipeline the steps are as follows:

1. Copy the `workflow/Snakefile` to the directory just below where your data is located: i.e. `ROOT` in the example above
2. Edit the `Snakefile` to contain the required parameters - datadir name, experiment name, conditions, tolerances, etc.
3. Run `snakemake -n` from the `ROOT` directory to make sure everything is ready to run
4. Run analysis by running `snakemake --cores 1` from the `ROOT` directory and wait for analysis to complete (or fail).

If you want to re-run an analysis you can either type `snakemake --cores 1 -f` or remove all the `*.done` files.
