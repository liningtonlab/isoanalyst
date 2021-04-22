# Example Workflow

Inside this directory is a fully functional working example of how to use IsoAnalyst
including using `snakemake` to run the full pipeline.

The input data contained herein as specified in the `input_spec.csv`, but are composed
of mzML files in the `mzmls` directory and CSV feature list files in the `feature_lists`
directory.

The data were prepared by taking a selection of about 10 features from one of the full publication
expirements, and extracting all the mass features within 20 Da M/Z and 0.06 min retention time
from the original all scan data in mzML files.
These data were then composed back into new (much smaller) mzML files using the [psims](https://github.com/mobiusklein/psims)
library.

These mzML files were then re-processed using MzMine (v2.53) to prepare feature lists using our
recommended generalized approach and exported as CSVs.'

Example output files are in the `EXAMPLE` directory, or you can install IsoAnalyst yourself
and re-run the analysis by running `snakemake -j1` from this directory.