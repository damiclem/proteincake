import pandas as pd
import argparse
import parser_lib

"""
Get the union of the passed proteins of the HMM and PSSM (Get proteins that passed PSSM or HMM).
"""
def ensamble(pssm_path, hmm_path, true_positive_path, out_path):
    # Get PSSM results
    df_pssm = pd.read_csv(pssm_path, index_col = [0])
    index_pssm = set(df_pssm.index)
    # Get HMM results
    df_hmm = pd.read_csv(hmm_path, index_col = [0, 1])
    index_hmm = set([index[0] for index in df_hmm.index])
    # Get true_positive
    true_positive_prot = parser_lib.fasta_to_dict(path=true_positive_path)
    true_positive = set([key.split('|')[1]for key in true_positive_prot.keys()])


    # Get the union of the two results
    ensamble_positive = index_pssm.union(index_hmm)
    ensamble_df = pd.DataFrame()
    PSSM = [x in index_pssm for x in ensamble_positive]
    HMM = [x in index_hmm for x in ensamble_positive]
    TP = [x in true_positive for x in ensamble_positive]
    ensamble_df['PSSM'] = PSSM
    ensamble_df['HMM'] = HMM
    ensamble_df['TP'] = TP
    ensamble_df.index = ensamble_positive
    ensamble_df.to_csv(out_path)


if __name__ == '__main__':
    # 1. Define the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--pssm_df_path',   type=str,   default='results/PSSM_df.csv')
    parser.add_argument('--hmm_df_path',    type=str,   default='results/HMM_df.csv')
    parser.add_argument('--true_positive',  type=str,   default='data/positive_pfam.fasta')
    parser.add_argument('--out_path',       type=str,   default='results/ensamble_df.csv')

    # 2. Generate dictionary of args
    args = parser.parse_args()

    # 3. Run the pipeline
    ensamble(args.pssm_df_path, args.hmm_df_path, args.true_positive, args.out_path)
