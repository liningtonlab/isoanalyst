from pathlib import Path

from core import cppis_masterlist, isotope_scraper


def main():
    '''User provided arguments'''

    # Input directory

    # source_dir = '/Volumes/TimeMachineBackups/201908_RLUS1353/20190808/seven/'
    # source_dir = Path("./IsotopedataProcessing/J1").absolute()
    source_dir = Path("./IsotopedataProcessing/J2").absolute()


    # Experiment name or organism name, applies to whole dataset
    # for naming output files
    exp_name = 'RLUS1353'


    # Experimental conditions
    # primary - positional argument always a list of labels used preferably any 3 letter code per label
    # secondary - optional argument for a nested condition within the labels
    #             ie days, media, additive, etc
    # Depending on which arguments are given, a list called 'conditions' is generated and used throughout
    # all steps of the code including file names and directory/folder naming.

    sec_arg = 0

    if sec_arg == 1:
        primary = ['SB','SA','GB','GA']
        secondary = ['N1','N2','N3','N4']
        conditions = []
        for p in primary:
            for s in secondary:
                conditions.append(p+s)

    else:
        primary = ['ACE','PROP','MET','GLU']
        conditions = primary

    master = None # for continuing
    master = cppis_masterlist(source_dir, conditions, exp_name)
    isotope_scraper(source_dir, conditions, master=master)

if __name__ == "__main__":
    main()
