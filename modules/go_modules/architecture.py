import numpy as np
import pandas as pd
import argparse

"""
Function that extract the proteins belonging to a certain architecture and maps each protein to its
GO annotation. Accept in input the architecture of interest and two dataframes. The first one maps e
very protein to its architecture and requires columns 'entry_ac' and 'architecture', the second one
requires columns 'entry_ac' and 'go'.
"""
def select_architecture(arch, arch_df, go_df):
    assert arch in arch_df.architecture.values, 'Architecture not found'
    assert set(arch_df.entry_ac).issubset(set(go_df.entry_ac)), 'Architecture entries do not match'

    # 1. Select proteins that present the input architecture
    entries = arch_df.entry_ac[arch_df.architecture == arch]
    # 2. Retrieve GO for each protein
    return go_df.loc[go_df.entry_ac.isin(entries)]



"""
Loops over the architectures observed in 'original_arch' and apply 'select_architecture' to
each of them saving the outputs in BASE_PATH folder.
"""
def main(original_arch_path, ds_gene_ontology_path, base_path):

    # Load file that maps each protein to its architecture
    original_arch = pd.read_csv(original_arch_path, sep='\t', index_col=0)
    # Load file that maps each protein to its GO annotations
    ds_gene_ontology = pd.read_csv(ds_gene_ontology_path, sep='\t', dtype=str)

    # List all the observed architectures
    arch_list = list(set(original_arch.architecture))

    # Save tsv files for each architecture
    for arch in arch_list:
        go_arch_df = select_architecture(arch=arch, arch_df=original_arch, go_df=ds_gene_ontology)
        go_arch_df.to_csv(base_path+arch+'_arch.csv', sep='\t')



if __name__ == '__main__':

    # Define arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('--original_path',     type=str,   default='data/architecture/go_architectures/architecture_background.csv')
    parser.add_argument('--go_path',           type=str,   default='data/go/go.csv')
    parser.add_argument('--out_path',          type=str,   default='data/architecture/go_architectures/')

    # Define dictionary of args
    args = parser.parse_args()

    main(args.original_path, args.go_path, args.out_path)
