import json
import gzip
import copy
import argparse

import numpy as np
import pandas as pd

"""
Functions that create the datasets (target and background) for PDB network
"""
def pdb_target_dataset(original_dataset_path, do_dataset_path,
                       mapping_dataset_path, human_dataset_path):

    """
    1. original_dataset_path: the path of the target dataset, a.k.a the dataset of protein tha passed our model
    2. go_dataset_path: the path of the full go dataset (entry_ac, go_id)
    3. mapping_dataset_path: the path of the file mapping from pdb to uniprot
    4. human_dataset_path: the path of the full human dataset
    5. col_name_entry: the name of the column containing the proteins uniprot id
    6. col_name_pdb: the name of the column containing the proteins pdb id
    """

    ### 1. Load the mapping
    mapping_df = pd.read_table(mapping_dataset_path, header=1)
    mapping_df.columns = [col.lower() for col in mapping_df.columns]
    ### 2. Load the original dataset
    original_df = pd.read_table(original_dataset_path)
    ### 3. Load the do dataset
    background_df = pd.read_table(do_dataset_path, dtype=str, index_col=[0])
    ### 4. Load the human dataset
    human_df = pd.read_table(human_dataset_path)
    # 4.1 Take out the entry_ac that have a pdb_id
    protein_with_pdb = human_df.entry_ac[human_df.pdb_ids.isna() == False]
    # 4.2 Take out from  the background dataset the protein without a pdb
    background_df = background_df[background_df.entry_ac.isin(protein_with_pdb)]
    ### 5. Get all original proteins with a pdb and get all the other proteins which shares the same pdb
    # 5.1. Get a dataset with key (uniprot_id, pdb_id)
    values = []
    for n in range(original_df.shape[1]):
        key = original_df.loc[n, 'entry_ac']
        value = original_df.loc[n, 'pdb_ids']
        if type(value) == str:
            pdb_ids = value.split(';')[:-1]
            for ids in pdb_ids:
                values.append([key, ids.lower()])
    pdb_original = pd.DataFrame(values, columns=['entry_ac', 'pdb_ids'])
    # 5.2 Merge the new dataset with the mapping df to get all the proteins with that pdb id
    target_dataset = pd.merge(pdb_original, mapping_df, left_on='pdb_ids', right_on='pdb', how='left')
    ### 6. Get the GO of every pdb_in our target_dataset
    target_dataset = background_df[background_df.entry_ac.isin(target_dataset['sp_primary'])]
    ###

    return target_dataset, background_df

if __name__ == '__main__':
    # 1. Define arguments
    parser = argparse.ArgumentParser()
    ### Paths of the two dataframe that must be compared
    parser.add_argument('--original_df_path',     type=str,   default='data/results/jackhmmer.tsv')
    parser.add_argument('--do_df_path',           type=str,   default='data/do/do.csv')
    parser.add_argument('--mapping_df_path',      type=str,   default='data/pdb/pdb_chain_uniprot.tsv')
    parser.add_argument('--human_df_path',        type=str,   default='data/human.csv')
    ### Out directory of the results and of the WordCloud
    parser.add_argument('--out_path_target',      type=str,   default='data/pdb/pdb_target_do.csv')
    parser.add_argument('--out_path_background',  type=str,   default='data/pdb/pdb_background_do.csv')


    # 2. Define dictionary of args
    args = parser.parse_args()

    # 3. Get the target and background datasets
    df_target, df_background = pdb_target_dataset(original_dataset_path=args.original_df_path,
                                                  do_dataset_path=args.do_df_path,
                                                  mapping_dataset_path=args.mapping_df_path,
                                                  human_dataset_path=args.human_df_path)

    # 4. Save the Dataframe
    df_target.to_csv(args.out_path_target, sep='\t')
    df_background.to_csv(args.out_path_background, sep='\t')
