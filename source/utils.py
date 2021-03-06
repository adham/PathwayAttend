"""
Utility functions

Adham Beyki
Deakin UNiversity
odinay@gmail.com
2018-11-15
"""

import numpy as np
import pandas as pd


def prep_genexp_data(genexp_path, pam50_path, sample_filter_path):
    """ Reads the data and returns a dataframe of gene expresession and
    labels of samples.
    """

    # read data files
    genexp_df = pd.read_csv(genexp_path)
    genexp_df.rename(columns={'Unnamed: 0': 'GENE-ID'}, inplace=True)
    pam50_df = pd.read_csv(pam50_path)
    sample_filter = pd.read_csv(sample_filter_path, header=None)


    # filter samples to keep PAM50 samples 
    col_idxs = sample_filter.values.astype(bool).flatten()
    col_idxs = np.hstack([np.array(True), col_idxs])    # insert a True at the beginning to keep GENE-ID column
    cols = genexp_df.columns[col_idxs]
    genexp_df = transpose_df(genexp_df[cols])    

    # remove gene ENSG00000179546 as it has NULL values for some samples
    idxs = genexp_df.isna().sum() != 0
    cols = genexp_df.columns[idxs].tolist()
    genexp_df.drop(cols, axis=1, inplace=True)

    # check if filtering the samples and their order is correct
    assert (genexp_df['Sample-ID'].apply(lambda s: '-'.join(s.split('-')[:4])[:-1]) == pam50_df["Sample ID"]).sum() == len(pam50_df)

    # only consider LumA and LumB
    idxs = pam50_df['PAM50'].apply(lambda s: s in ['LumA', 'LumB'])
    data_df = genexp_df[idxs].copy()
    data_df['PAM50'] = pam50_df[idxs]['PAM50']

    data_df.reset_index(inplace=True, drop=True)
    
    return data_df