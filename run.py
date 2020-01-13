from pathlib import Path

from isotracer import cppis_masterlist, isotope_label_detector, isotope_scraper

# TODO: Need to address RuntimeWarning of invalid Degrees of Freedom in ttest_ind

def main():
    '''User provided arguments'''
    # TODO: Generate CLI for code

    # Input directory
    source_dir = Path("./isotope_experiments_data/SERY").absolute()

    # Experiment name or organism name, applies to whole dataset
    # for naming output files
    exp_name = 'SERY'

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
    isotope_scraper(source_dir, conditions, exp_name, master=master)
    isotope_label_detector(source_dir, conditions, master=master)


if __name__ == "__main__":
    main()
