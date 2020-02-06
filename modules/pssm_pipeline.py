"""
Pipeline PSSM
"""
import os
import pandas as pd
import argparse
import parser_lib as prs

def pssm_pipeline(pfam_path, model_path, out_results_path, df_out_path,
                  # outfmt=6,
                  num_iteration=3, evalue=0.05):
    """
    Requires: os, pandas as pd
    Required files: the MSA fasta file and the test set fasta file

    Input file:
    1. pfam_path:           the path of the test set(file required). str
    2. model_path:          the path where we save the PSSM model (file generated). str
    3. out_results_path:    the path where we save the PSSM results (file generated). str
    4. df_out_path:         the path where we save the DataFrame (file generated). str
    """
    # 1. Build database
    os.system('makeblastdb -dbtype prot -in {} -parse_seqids'.format(pfam_path))
    print('Database created')
    # 2. Run psiblast
    os.system('psiblast -in_pssm {} -db {} -num_iterations {} -evalue {} -out {}'.
              format(model_path, pfam_path, num_iteration, evalue, out_results_path))
    print('Psiblast done')
    #3. Create dataframe
    df = prs.pssm_parser(path=out_results_path)
    df.to_csv(df_out_path)
    print('Results saved')


### MAIN
if __name__ == '__main__':

    # 1. Define the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--pfam_path',          type=str,   default='data/all_pfam.fasta')
    parser.add_argument('--model_path',         type=str,   default='models/model.pssm')
    parser.add_argument('--out_results_path',   type=str,   default='results/PSSM_results.txt')
    parser.add_argument('--df_out_path',        type=str,   default='results/PSSM_df.csv')
    # parser.add_argument('--outfmt', type=int, default=6)
    parser.add_argument('--num_iteration',      type=int,   default=1)
    parser.add_argument('--evalue',             type=float, default=0.001)

    # 2. Generate dictionary of args
    args = parser.parse_args()

    # 3. Run the pipeline
    pssm_pipeline(args.pfam_path, args.model_path, args.out_results_path,
                  args.df_out_path,
                  # args.outfmt,
                  args.num_iteration, args.evalue)
