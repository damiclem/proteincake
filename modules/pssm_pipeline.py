"""
Pipeline PSSM
"""
import os
import pandas as pd
import argparse
import parser as prs

def pssm_pipeline(pfam_path, model_path, out_results_path, df_out_path,
                  # outfmt=6,
                  num_iteration=3, evalue=0.05):
    """
    Requires: os, pandas as pd
    Required files: the MSA fasta file and the test set fasta file

    Input file:
    1. msa_path:       the path of the msa.fasta file (file required). str
    2. test_set_path:  the path of the test_set (file required). str
    3. hmm_model_path: the path where we save the HMM model (file generated). str
    4. out_file_path:  the path where we save the first result table (file generated). str
    5. domtblout_path: the path where we save the second result table (file generated). str
    6. df_out_path:    the path where to save the final DataFrame (file generated). str
    7. evalue:         the maximum e-value that a protein/domain need to have in order to be inserted in the final DataFrame.
                       Default is 0.05. float
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



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pfam_path', type=str, default='data/all_pfam.fasta')
    parser.add_argument('--model_path', type=str, default='models/model.pssm')
    parser.add_argument('--out_results_path', type=str, default='results/PSSM_results.txt')
    parser.add_argument('--df_out_path', type=str, default='results/PSSM_df.csv')
    # parser.add_argument('--outfmt', type=int, default=6)
    parser.add_argument('--num_iteration', type=int, default=3)
    parser.add_argument('--evalue', type=float, default=0.005)

    args = parser.parse_args()

    pssm_pipeline(args.pfam_path, args.model_path, args.out_results_path,
                  args.df_out_path,
                  # args.outfmt,
                  args.num_iteration, args.evalue)
