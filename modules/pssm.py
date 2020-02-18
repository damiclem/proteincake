############
### PSSM ###
############


# Dependencies
import sys
import subprocess
import pandas as pd
import argparse


# Parse PSSM string result to Pandas DataFrame
def parse(result, raw=False):
    """
    Input:
    - PSSM result, formatted as string
    Output:
    - PSSM result, formatted as Pandas DataFrame
    """
    # Create DataFrame
    result = pd.DataFrame([row.split('\t') for row in result.split('\n') if row != ''][:-1])
    # Set column names
    result.columns = ['query_acc_ver', 'subject_acc_ver', 'perc_identity',
                      'align_length', 'mismatches', 'gap_opens', 'q_start', 'q_end',
                      's_start', 's_end', 'evalue', 'bit_score']
    # Format output if required
    if not raw:
        # Get subset of columns: sequence accession, start, end, evalue
        result = result[['subject_acc_ver', 's_start', 's_end', 'evalue']]
        # Map columns
        result.columns = ['entry_ac', 'seq_start', 'seq_end', 'e_value']
    # Return DataFrame object
    return result


# Fit the model (generate pssm)
def fit(blast_path, msa_path, model_path=None):
    # Run PSSM creration
    pssm = subprocess.run(['psiblast', '-subject', blast_path, '-in_msa', msa_path, '-out_pssm', model_path],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # Return output
    return pssm.stdout.decode('utf-8')


# Evaluate the model against a test set
def test(model_path, test_path, out_fmt=6, num_iterations=3, e_value=0.05):
    """
    Creates a blast db, the runs psi-blast with a given PSSM
    Errors are handled by catching SubprocessError exception
    More info about subprocess module: https://docs.python.org/3/library/subprocess.html
    By default, output is retrieved in tabular format (outfmt=6), which has the following columns:
    query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q.start, q.end, s.start, s.end, evalue, bit score
    Input:
        1. blast_path:      path to input blast file (FASTA format)
        2. model_path:      path to output model
        3. test_path:       path to test set (FASTA format)
        4. out_fmt:         output format (default 6, tabular)
        5. num_iterations:  number of psi-blast iterations
        6. e_value:         e-value threshold
    Output:
        1. text output retrieved from psiblast
    """
    # Build database
    subprocess.run(['makeblastdb', '-dbtype', 'prot', '-in', test_path, '-parse_seqids'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # Run psi-blast
    out = subprocess.run(['psiblast', '-in_pssm', model_path, '-db', test_path,
                          '-num_iterations', str(num_iterations),
                          '-evalue', str(e_value), '-outfmt', str(out_fmt)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # Turn tabular output as string
    return out.stdout.decode('utf-8')


### MAIN
if __name__ == '__main__':

    # 1. Define the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--fit',            type=bool,  default=True)
    parser.add_argument('--blast_path',     type=str,   default='data/blast.fasta')
    parser.add_argument('--msa_path',       type=str,   default='data/msa.edited.fasta')
    parser.add_argument('--model_path',     type=str,   default='models/model.pssm')
    parser.add_argument('--test_path',      type=str,   default='data/human.fasta')
    parser.add_argument('--out_path',       type=str)
    parser.add_argument('--out_fmt',        type=int,   default=6)
    parser.add_argument('--num_iterations', type=int,   default=1)
    parser.add_argument('--e_value',        type=float, default=0.001)

    # 2. Generate dictionary of args
    args = parser.parse_args()

    # 3. Run psi-blast
    try:

        # Fit model, if requested
        if args.fit:
            # Create new pssm using msa as input
            fit(args.blast_path, args.msa_path, args.model_path)

        # Run psi-blast with given PSSM
        psi_blast = test(args.model_path, args.test_path, out_fmt=args.out_fmt,
                         num_iterations=args.num_iterations, e_value=args.e_value)
        # Create dataset
        psi_blast = parse(psi_blast)

        # Case out path has been set: write to file
        if args.out_path:
            # Write to file
            psi_blast.to_csv(args.out_path, index=False, sep='\t')
        # Case out path has not been set: write to console
        else:
            print(psi_blast.to_string())

    # Catch eventual errors
    except subprocess.CalledProcessError as e:
        # Exit with error (last character is newline)
        sys.exit(e.stderr.decode('utf-8').strip())
