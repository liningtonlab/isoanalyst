## Welcome to IsoAnalyst


## CLI Options

## Input Specification

## Data Requirements

#### 1. Feature lists

Generic CSV tabular inputs minimally with `["Sample", "PrecMz", "PrecZ", "PrecIntensity", "RetTime", "ScanLowRange", "ScanHighRange"]` column headers. These are compatible msExpress CPPIS files.

Working on MzMine2 import pipeline.

MzMine outputs only require the `["row m/z", "row retention time"]` columns. 

__Question:__ Intensity? Use MzMine Peak Area?

The program will then set a window of 10 scans on both sides of the center scan index detected during all scan importing.


#### 2. Scan Data

Generic CSV tabular inputs minimally with `["FunctionScanIndex", "RT", "MZ", "Intensity"]` column headers.

Can import func0001 CSVs from msExpress.

Can import mzMLs from [msConvert](http://proteowizard.sourceforge.net/tools/msconvert.html).

### Development and Testing

To run the tests, you must first install dev dependencies `conda install pytest`.


#### Data Collection Recommendation

Not supported for DDA.

Either DIA with interleaved scans or MS1 only.