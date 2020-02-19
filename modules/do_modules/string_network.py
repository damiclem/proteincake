import json
import gzip
import copy
import argparse

import numpy as np
import pandas as pd

def load(path, keep=None, sep=' '):
    # Load dataset
    string_ds = pd.read_csv(path, compression='gzip', header=0, sep=sep)
    # Subset pandas dataframe: keep only some rows
    if keep is not None:
        string_ds = string_ds[string_ds.protein1.isin(keep) | string_ds.protein2.isin(keep)]
    # Return retrieved dataset
    return string_ds

if __name__ == '__main__':
    # 1. Define arguments
    parser = argparse.ArgumentParser()
    ### Paths of the two dataframe that must be compared
    parser.add_argument('--original_path',        type=str,   default='data/results/jackhmmer.tsv')
    parser.add_argument('--do_path',              type=str,   default='data/do/do.csv')
    parser.add_argument('--string_path',          type=str,   default='data/string/string.txt.gz')
    parser.add_argument('--human_path',           type=str,   default='data/human.csv')
    ### Out directory of the results and of the WordCloud
    parser.add_argument('--out_path_target',      type=str,   default='data/string/string_target_do.csv')
    parser.add_argument('--out_path_background',  type=str,   default='data/string/string_background_do.csv')
    ### Names of the columns of the input dataframe containing the GO_id
    parser.add_argument('--col_name_entry_ac',    type=str,   default='entry_ac')
    ### Minimum Combined Score of the interaction
    parser.add_argument('--min_score',            type=int,   default=700)

    # 2. Define dictionary of args
    args = parser.parse_args()

    # 3. Load the required datasets
    # Human
    human =  pd.read_csv(args.human_path, sep='\t')
    human = human[human.string_id.isna() == False] #Remove proteins with no STRING
    human.string_id = human.string_id.map(lambda x: str(x).replace(';', '').strip())
    # Original
    original = pd.read_csv(args.original_path, sep='\t')
    # Do
    do = pd.read_csv(args.do_path, sep='\t', dtype=str)
    # String
    string = load(args.string_path)

    # 4. Merge the STRING interactors
    original_string_ids = set([i[:-1] for i in original.string_id.tolist() if type(i) == str])
    # Get direct interactors
    original_interaction = string[string.protein1.isin(original_string_ids)]
    # Filter by score
    original_interaction = original_interaction[original_interaction.combined_score > args.min_score]
    # Define interactors ids
    interactors_string_ids = set(original_interaction.protein2.tolist())
    # Define union of the two sets
    all_string_ids = original_string_ids | interactors_string_ids
    # Get all proteins in original dataset, plus direct interactors
    original = human[human.string_id.isin(all_string_ids)]

    # 5. Get target and background and save it
    # String target DO dataset
    string_target_do = do[do.entry_ac.isin(original.entry_ac)]
    string_target_do.to_csv(args.out_path_target, sep='\t')

    # String background DO dataset
    string_background_do = do[do.entry_ac.isin(human.entry_ac)]
    string_background_do.to_csv(args.out_path_background, sep='\t')
