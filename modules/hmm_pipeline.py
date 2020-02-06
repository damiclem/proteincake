"""
Pipeline HMM
"""
import os
import pandas as pd
import argparse
import parser_lib as prs

"""
Pipeline that start from a MSA file, build an HMM model and test it using a test set. Requires MSA
"""
def hmm_pipeline(msa_path, test_set_path, hmm_model_path, out_file_path, domtblout_path, df_out_path, evalue=0.05):
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
    ### 1. Build the model
    os.system('hmmbuild {} {}'.format(hmm_model_path, msa_path))
    print("Model has been built")
    ### 2. Test the model
    os.system('hmmsearch --domtblout {} -o {} {} {}'.
             format(domtblout_path, out_file_path, hmm_model_path, test_set_path))
    print("Model has been tested")
    ### 3. Get dataframe
    df = prs.hmm_parser(path=domtblout_path)
    df = df[df['i-Evalue'] < evalue]
    df.to_csv(df_out_path)
    print("Done!")



if __name__ == '__main__':

    # 1.Define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa_path',       type=str,   default='data/msa.fasta')
    parser.add_argument('--test_set_path',  type=str,   default='data/all_pfam.fasta')
    parser.add_argument('--hmm_model_path', type=str,   default='models/model.hmm')
    parser.add_argument('--out_file_path',  type=str,   default='results/HMM_results')
    parser.add_argument('--domtblout_path', type=str,   default='results/HMM_results_domtblout')
    parser.add_argument('--df_out_path',    type=str,   default='results/HMM_df.csv')
    parser.add_argument('--evalue',         type=float, default=0.005)

    # 2.Define dictionary of args
    args = parser.parse_args()

    # 3.Run Pipeline
    hmm_pipeline(args.msa_path, args.test_set_path, args.hmm_model_path, args.out_file_path,
                 args.domtblout_path, args.df_out_path, args.evalue)
