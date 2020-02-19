###########
### HMM ###
###########


# Dependencies
import sys
import subprocess
import tempfile
import pandas as pd
import numpy as np
import argparse


# Available algorithms
JACKHMMER = 'jackhmmer'
HMMSEARCH = 'hmmsearch'


# Parse HMM string result to Pandas DataFrame
def parse(result, raw=False):
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
    # Turn result table into DataFrame Object
    result = pd.DataFrame(result, columns=columns)
    # Format dataframe: return predicted sequence accession, start, end, 'evalue'
    if not raw:
        # Return subset of dataset columns
        result = result[['target_name', 'align_from', 'align_to', 'dom_i_evalue']]
        # Map columns
        result.columns = ['entry_ac', 'seq_start', 'seq_end', 'e_value']
    # Slice the returned dataset in order to retrieve only
    return result


# Fit HMM model (use hmmbuild)
def fit(msa_path, model_path):
    # Build the model
    out = subprocess.run(['hmmbuild', model_path, msa_path],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         check=True)
    return out.stdout.decode('utf-8')


# Evaluate HMM model, return result
def test(algorithm=HMMSEARCH, *args, **kwargs):
    """
    Test an HMM model
    Input:
        1. algorithm:       which algorithm to be used in evaluation
        2. *args, **kwargs: other arguments passed to inner functions
    Output:
        1. text output retrieved from psiblast
    """
    # Test the given model
    if algorithm == HMMSEARCH:
        out = test_hmmsearch(*args, **kwargs)
    elif algorithm == JACKHMMER:
        out = test_jackhmmer(*args, **kwargs)
    # Turn bytes output to string
    return out.stdout.decode('utf-8')


# Evaluate HMM model using hmmsearch (requires fit)
def test_hmmsearch(model_path, test_path):
    """
    Input:
        1. model_path:  path to fitted HMM model
        2. test_path:   path to test set (.fasta formatted)
    Output:
        1. text output retrieved from psiblast
    """
    # Test the given model
    return subprocess.run([HMMSEARCH,
                          '-o', '/dev/null',  # Output is silenced
                          '--domtblout', '/dev/stdout',  # Domtblout redirected to stdout
                          model_path, test_path],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          check=True)

# Evaluate HMM model using jackhmmes (does not require fit)
def test_jackhmmer(seq_path, test_path):
    """
    Input:
        1. seq_in:      input sequence
        2. test_path:   path to test set (.fasta formatted)
    Output:
        1. text output retrieved from psiblast
    """
    # Test the given model
    return subprocess.run([JACKHMMER, '-o', '/dev/null', '--domtblout', '/dev/stdout', seq_path, test_path],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)


if __name__ == '__main__':

    # 1. Define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--algorithm',      type=str,   default=HMMSEARCH)
    parser.add_argument('--fit',            type=str,   default=True)
    parser.add_argument('--seq_path',       type=str,   default='data/domain.fasta')
    parser.add_argument('--msa_path',       type=str,   default='data/msa.fasta')
    parser.add_argument('--test_path',      type=str,   default='data/human.fasta')
    parser.add_argument('--model_path',     type=str,   default='models/model.hmm')
    parser.add_argument('--out_path',       type=str)
    parser.add_argument('--e_value',        type=float, default=0.001)

    # 2. Define dictionary of args
    args = parser.parse_args()

    # 3. Run HMM model test
    try:

        # Define output container
        hmm_out = None

        # Case HMMSEARCH algorithm
        if args.algorithm == HMMSEARCH:
            # Fit the model if required
            if args.fit:
                fit(msa_path=args.msa_path, model_path=args.model_path)
            # Evaluate the model
            hmm_out = test(algorithm=HMMSEARCH, model_path=args.model_path, test_path=args.test_path)


        # Case JACKHMMER algorithm
        elif args.algorithm == JACKHMMER:
            hmm_out = test(algorithm=JACKHMMER, seq_path=args.seq_path, test_path=args.test_path)

        # Error: no algorithm has been chosen
        else:
            sys.exit('Error: no valid algorithm has been chosen')

        # Parse output to pandas DataFrame object
        hmm_out = parse(hmm_out)

        # Filter on e-value
        hmm_out = hmm_out[hmm_out.e_value.astype(np.float) < float(args.e_value)]

        # Case out path has been set: write to file
        if args.out_path:
            # Write to file
            hmm_out.to_csv(args.out_path, index=False, sep='\t')
        # Case out path has not been set: write to console
        else:
            print(hmm_out.to_string())

    # Catch eventual errors
    except subprocess.CalledProcessError as e:
        # Exit with error (last character is newline)
        sys.exit(e.stderr.decode('utf-8'))
