import pandas as pd
import numpy as np
import scipy.stats
from sklearn.linear_model import LinearRegression

def scan_data_QC(df,c,e):
    # if only one scan is detected
    if isinstance(df, pd.Series):
        r = 'not enough scans'
        to_log(e,c,r)
        return False

    # if ion only detected in unlabeled or labeled
    elif len(df['Isotope'].unique()) < 2:
        r = 'not in both isotope conditions'
        to_log(e,c,r)
        return False

    return True


# should probably establish headers for this log file somewhere in the function
def to_log(e,c,r):
    """Custom logging by Cat

    Args:
        e (str): Experiment ID
        c (str): Condition
        r (str): Reason / message
    """
    log = 'log_file.txt'
    with open(log, 'a') as logfile:
        logfile.write(e+','+c+','+r+'\n')



# will run this function per every replicate condition ie 12ACED0-1
#  compare all Mi v Mj combinations and return a table with statistics
def MID_ratio_comparison(ion_rep):
    Mn = ion_rep['Isotopomer'].unique().tolist() # list of unique isotopomers M0,M1,M2,etc
    Exp_ID = ion_rep.index[0]
    Condition = ion_rep['Condition'].iloc[0]
    Isotope = ion_rep['Isotope'].iloc[0]

    ion_rep = ion_rep.drop_duplicates(['FunctionScanIndex','Isotopomer'])
    # set index to scan number
    new_ix = ion_rep.set_index(['FunctionScanIndex'])
    # slice index containing only relevant scan #s from intersect
    #scan_sliced = new_ix.ix[scan_intsct]
    # slice out only the relevant data columns - Isotopomer and Intensity
    data_sliced = new_ix[['Intensity','Isotopomer']]
    # create hierarchical index
    stacked = data_sliced.set_index(['Isotopomer'], append = True)
    # unstack so columns are MO/M1 and the values are all intensity
    unstacked = stacked.unstack()
    # Eliminate 'Intenisty' from column name to simplify slicing
    final_munged = unstacked['Intensity']
    cols = pd.Index(['Exp_ID','Condition','Isotope','Isotopomer','Slope','Y_intercept','R_squared','Scan_count'])
    header = pd.DataFrame(columns = cols)


    for m in range(len(Mn)):
        if Mn[m] != Mn[-1]:
            Mi_scans = ion_rep.loc[ion_rep['Isotopomer'] == Mn[m]]['FunctionScanIndex'] # get M0 scans
            Mj_scans = ion_rep.loc[ion_rep['Isotopomer'] == Mn[m+1]]['FunctionScanIndex'] # get M1 scans

            # intersect the scan #s so that only scans containing data are compared
            scan_intsct = np.intersect1d(Mi_scans, Mj_scans)
            MiMj_scans = final_munged.ix[scan_intsct]
            scan_count = len(MiMj_scans)
            if scan_count <= 2:
                r = 'scan intersect too small for' + Mn[m] + 'v' + Mn[m+1]
                to_log(Exp_ID,Condition,r)

            elif scan_count >= 3:

                model = LinearRegression()
                one = MiMj_scans.loc[:,[Mn[m]]]  # M0 data slice (x)
                two = MiMj_scans.loc[:,[Mn[m+1]]]  # M1 data slice (y)
                model.fit(one,two)  # fit the model - must do this first to use some other LinReg commands

                slope = float(model.coef_) # slope
                y_int = float(model.intercept_) # y int
                R2 = model.score(one,two) # R2 value- clearly not a good peak

                Isotopomer = Mn[m] + 'v' + Mn[m+1]

                data = {'Exp_ID' : Exp_ID,
                       'Condition' : Condition,
                       'Isotope' : Isotope,
                       'Isotopomer' : Isotopomer,
                       'Slope' : slope,
                       'Y_intercept' : y_int,
                       'R_squared' : R2,
                       'Scan_count' : scan_count}


                header = header.append(data, ignore_index=True)

    return header


