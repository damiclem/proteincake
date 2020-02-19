#########################
### TMALIGN ALGORITHM ###
#########################


# Dependencies
import numpy as np
import pandas as pd
import re
import subprocess
import argparse
import sys


# Constants
SCRIPT_PATH = r'./resources/TMalign'  # Defualt path to script


# Parse result
def parse(out):

    # Define chains container
    chains = dict()
    # Get chains data
    for i in range(2):
        # Define current chain container
        chains.setdefault(i, dict())
        # Get chain name
        chains[i]['name'] = str(re.search(
            r'Name of Chain_{:d}:(.+?)\n'.format(i+1), out
        ).group(1))
        # Get chain length
        chains[i]['len'] = int(re.search(
            r'Length of Chain_{:d}: (\d+)?'.format(i+1), out
        ).group(1))

    # Get TM score for each chain
    tm_score = re.findall(
        r'TM-score= ([\d.]+)?', out
    )
    # Store TM scores separately
    chains[0]['tm_score'] = float(tm_score[0])
    chains[1]['tm_score'] = float(tm_score[1])

    # Get normalized TM score
    tm_score = float(tm_score[2])

    # Get aligned residues length
    align_len = int(re.search(
        r'Aligned length=[ ]+(\d+)?', out
    ).group(1))
    # Get root mean squared deviation (RMSD)
    rmsd = float(
        re.search(r'RMSD=[ ]+([\d.]+)?', out
    ).group(1))
    # Get n. identical / n. aligned
    seq_id = float(re.search(
        r'Seq_ID=n_identical/n_aligned=[ ]+([\d.]+)?', out
    ).group(1))

    # Return retrieved data
    return chains, align_len, rmsd, tm_score, seq_id


# Execute alignment
def align(pdb1_path, pdb2_path,script_path=SCRIPT_PATH):
    # Run script, store result
    out = subprocess.run([script_path, '-a', 'T', pdb1_path, pdb2_path],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         check=True)
    # Parse alignment result
    return out.stdout.decode('utf-8')


# Execute multiple pairwise alignment
def multi_align(pdb_paths, triangular=True, script_path=SCRIPT_PATH):
    # Define index (get only file name)
    index = [re.search(r'/(\w+)?[.\w]+$', pdb_path).group(1) for pdb_path in pdb_paths]
    # Define number of PDBs
    n = len(pdb_paths)
    # Define a squared matrix containing scores for each PDB pair
    rmsd_mat = np.zeros(shape=(n, n), dtype=np.float)  # RMSD
    tm_mat =  rmsd_mat.copy()  # TM score

    # Fill matrix
    for i in range(n):  # Rows
        for j in range(n):  # Columns
            # Check in triangular
            if triangular and j <= i: continue
            # Execute pairwise alignment
            align_out = align(pdb_paths[i], pdb_paths[j], script_path=script_path)
            # Parse alignment results
            chains, _, rmsd, tm_score, _ = parse(align_out)
            # Store scores
            tm_mat[i, j] = tm_mat[j, i] = 1 - tm_score
            rmsd_mat[i, j] = rmsd_mat[j, i] = rmsd

    # Return filled matrices
    return (pd.DataFrame(data=rmsd_mat, index=index, columns=index),
            pd.DataFrame(data=tm_mat, index=index, columns=index))


# Execute (multiple) pairwise structural alignment
if __name__ == '__main__':

    # 1. Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_paths', type=str, nargs='+', required=True)
    parser.add_argument('--script_path', type=str, default=SCRIPT_PATH)
    parser.add_argument('--out_type', type=str, default='rmsd')
    parser.add_argument('--out_path', type=str)
    args = parser.parse_args()

    # 2. Run alignment
    try:

        # Get rmsd and tm score
        rmsd, tm_score = multi_align(args.pdb_paths)
        # Get output scores
        score = tm_score if args.out_type == 'tmscore' else rmsd

        # If output path is defined, write output to file
        if args.out_path:
            score.to_csv(args.out_path, sep='\t')
        # Output path not defined: print out results
        else:
            print(score.to_string())

    # 3. Catch eventual errors
    except subprocess.CalledProcessError as e:
        # Exit with error (last character is newline)
        sys.exit(e.stderr.decode('utf-8').strip())
