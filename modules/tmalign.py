#########################
### TMALIGN ALGORITHM ###
#########################


# Dependencies
import numpy as np
import pandas as pd
import re
import subprocess


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
    tm_scores = re.findall(
        r'TM-score= ([\d.]+)?', out
    )
    # Store TM scores separately
    chains[0]['tm_score'] = float(tm_scores[0])
    chains[1]['tm_score'] = float(tm_scores[1])

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
    return chains, align_len, rmsd, seq_id


# Execute alignment
def align(pdb1_path, pdb2_path, script_path=SCRIPT_PATH):
    # Run script, store result
    out = subprocess.run([script_path, pdb1_path, pdb2_path],
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
    rmsd_scores = np.zeros(shape=(n, n), dtype=np.float)  # RMSD
    tm_scores =  rmsd_scores.copy()  # TM score

    # Fill matrix
    for i in range(n):  # Rows
        for j in range(n):  # Columns
            # Check in triangular
            if triangular and j <= i: continue
            # Execute pairwise alignment
            align_out = align(pdb_paths[i], pdb_paths[j], script_path=script_path)
            # Parse alignment results
            chains, _, rmsd, _ = parse(align_out)
            # Sotre scores
            tm_scores[i, j] = chains[0]['tm_score']
            rmsd_scores[i, j] = rmsd

    # Return filled matrices
    return (pd.DataFrame(data=rmsd_scores, index=index, columns=index),
            pd.DataFrame(data=tm_scores, index=index, columns=index))
