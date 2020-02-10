###########
### HMM ###
###########


# Dependencies
import sys
import subprocess
import pandas as pd
import argparse


# Parse HMM string result to Pandas DataFrame
def parse(result):
    """
    Input:
    1. PSSM result (string)
    Output:
    1. PSSM result (Pandas DataFrame)
    """
    # Turn result from string to list (rows) of lists (columns)
    result = [row.split() for row in result.split('\n')
              if row != '' and row[0] != '#']
    # Define column names
    columns = ['target_db', 'target_name', 'target_gene', 'target_accession', 'target_len',
               'query_name', 'query_accession', 'query_len',
               'full_seq_evalue', 'full_seq_score', 'full_seq_bias',
               'dom_nr', 'dom_of', 'dom_c_evalue', 'dom_i_evalue', 'dom_score', 'dom_bias',
               'hmm_from', 'hmm_to', 'align_from', 'align_to', 'env_from', 'env_to',
               'acc', 'description']
    # Define number of columns and number of rows
    num_rows, num_cols = len(result), len(columns)
    # Update each row
    for i in range(num_rows):
        # First column must be splitted
        result[i] = result[i][0].split('|') + result[i][1:]
        # Last column must be grouped together: it is the description
        result[i][num_cols] = ' '.join(result[i][num_cols:])
        result[i] = result[i][:num_cols]
    # Turn result table into DataFrame Object and return it
    return pd.DataFrame(result, columns=columns)


# Run HMM script, return result
def run(msa_path, test_path, model_path):
    """
    Creates and tests an HMM model
    Input:
    1. test_path:   path to test set (.fasta formatted)
    1. msa_path:    path to multiple seuqence alignment file (.fasta formatted)
    3. model_path:  path where we save the HMM model
    4. out_path:    path where to save output
    7. evalue:      maximum e-value threshold (default 0.05)
    Output:
    1. Text output retrieved from psiblast
    """
    # Build the model
    subprocess.run(['hmmbuild', model_path, msa_path],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # Test the model
    out = subprocess.run(['hmmsearch',
                          '-o', '/dev/null',  # Output is silenced
                          '--domtblout', '/dev/stdout',  # Domtblout redirected to stdout
                          model_path, test_path],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         check=True)
    # Turn bytes output to string
    return out.stdout.decode('utf-8')


if __name__ == '__main__':

    # 1. Define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa_path',       type=str,   default='data/msa.fasta')
    parser.add_argument('--test_path',      type=str,   default='data/human.fasta')
    parser.add_argument('--model_path',     type=str,   default='models/model.hmm')
    parser.add_argument('--out_path',       type=str)
    parser.add_argument('--evalue',         type=float, default=0.005)

    # 2. Define dictionary of args
    args = parser.parse_args()

    # 3. Run HMM model test
    try:
        # Get hmm result
        hmm_result = run(args.msa_path, args.test_path, args.model_path,
                         out_path=args.out_path, evalue=args.evalue)
        hmm_result = parse(hmm_result)
        # Case out path has been set: write to file
        if args.out_path:
            # Write to file
            hmm_result.to_csv(args.out_path, index=False, sep='\t')
        # Case out path has not been set: write to console
        else:
            print(hmm_result.to_string())

    # Catch eventual errors
    except subprocess.CalledProcessError as e:
        # Exit with error (last character is newline)
        sys.exit(e.stderr.decode('utf-8'))
