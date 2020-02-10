############
### PSSM ###
############

# Dependencies
import sys
import subprocess
import pandas as pd
import argparse


# Parse PSSM string result to Pandas DataFrame
def parse(result):
    """
    Input:
    - PSSM result, formatted as string
    Output:
    - PSSM result, formatted as Pandas DataFrame
    """
    # Create DataFrame
    result = pd.DataFrame([row.split('\t') for row in result.split('\n') if row != ''])
    # Set column names
    result.columns = ['query_acc_ver', 'subject_acc_ver', 'perc_identity',
                      'align_length', 'mismatches', 'gap_opens', 'q_start', 'q_end',
                      's_start', 's_end', 'evalue', 'bit_score']
    # Return DataFrame object
    return result


# Run PSSM script, returns result
def run(test_path, model_path, outfmt=6, num_iterations=3, evalue=0.05):
    """
    Creates a blast db, the runs psi-blast with a given PSSM
    Errors are handled by catching SubprocessError exception
    More info about subprocess module: https://docs.python.org/3/library/subprocess.html
    By default, output is retrieved in tabular format (outfmt=6), which has the following columns:
    query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q.start, q.end, s.start, s.end, evalue, bit score
    Input:
    1. test_path:       Path to multiple sequence alignment, .fasta format
    2. model_path:      Path to model newly created model
    3. outfmt:          Output format (default tabular)
    4. num_iterations:  Number of psi-blast iterations
    5. evalue:          E-value threshold
    Output:
    1. Text output retrieved from psiblast
    """
    # Build database
    subprocess.run(['makeblastdb', '-dbtype', 'prot', '-in', test_path, '-parse_seqids'],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # Run psi-blast
    out = subprocess.run(['psiblast', '-in_pssm', model_path, '-db', test_path,
                          '-num_iterations', str(num_iterations),
                          '-evalue', str(evalue), '-outfmt', str(outfmt)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # Turn tabular output as string
    return out.stdout.decode('utf-8')


### MAIN
if __name__ == '__main__':

    # 1. Define the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--test_path',       type=str,   default='data/human.fasta')
    parser.add_argument('--model_path',      type=str,   default='models/model.pssm')
    parser.add_argument('--out_path',        type=str,   required=False)
    parser.add_argument('--outfmt',          type=int,   default=6)
    parser.add_argument('--num_iterations',  type=int,   default=1)
    parser.add_argument('--evalue',          type=float, default=0.001)

    # 2. Generate dictionary of args
    args = parser.parse_args()

    # 3. Run psi-blast
    try:
        # Run psi-blast with given PSSM
        pssm_result = run(args.test_path, args.model_path, outfmt=args.outfmt,
                          num_iterations=args.num_iterations, evalue=args.evalue)
        # Create dataset
        pssm_result = parse(pssm_result)
        # Case out path has been set: write to file
        if args.out_path:
            # Write to file
            pssm_result.to_csv(args.out_path, index=False, sep='\t')
        # Case out path has not been set: write to console
        else:
            print(pssm_result.to_string())
    # Catch eventual errors
    except subprocess.CalledProcessError as e:
        # Exit with error (last character is newline)
        sys.exit(e.stderr.decode('utf-8')[:-1])
