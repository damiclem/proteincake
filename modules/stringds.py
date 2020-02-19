######################
### STRING DATASET ###
######################


# Dependencies
import pandas as pd
import numpy as np
import gzip


# Load sting dataset as Pandas DataFrame object
def load(path, keep=None, sep=' '):
    # Load dataset
    string_ds = pd.read_csv('data/string.txt.gz', compression='gzip', header=0, sep=sep)
    # Subset pandas dataframe: keep only some rows
    if keep is not None:
        string_ds = string_ds[string_ds.protein1.isin(keep) | string_ds.protein2.isin(keep)]
    # Return retrieved dataset
    return string_ds


# Unit testing
if __name__ == '__main__':

    #print(load('data/string.txt.gz').head())
    string_ds = load('data/string.txt.gz', keep=['9606.ENSP00000380460',
                                                 '9606.ENSP00000292431',
                                                 '9606.ENSP00000322804'])
    print(string_ds)
