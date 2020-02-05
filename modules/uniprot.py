#################################
### Handle UniProt search API ###
#################################

# Dependencies
import requests as req
import time
import json
import xmltodict


# Constants
BASE_URL = r'https://www.uniprot.org'


# Define function for getting a single protein
def get_protein(protein_id, database='uniref', params={}):
    # Make API request
    res = req.get('/'.join([BASE_URL, database, str(protein_id) + '.fasta']),
                  headers={'Accept': 'text/x-fasta'}, 
                  params=params)
    # Return results
    return res.status_code == 200, res.text, res


"""
# Define function for retrieving a protein sequence from uniprot
def get_proteins(proteins_id, params={}, batch_size=100, accept='text/x-fasta'):
    # Result container
    out = ''
    # Define paramteres
    params={**{'offset': 0, 'size': -1},  # Default parameters
            **params}  # User defined params
    # Define total number of prpteins
    proteins_num = len(proteins_id)
    # Loop through each protein id
    for i in range(0, proteins_num, batch_size):
        # Get batch
        proteins_id_ = proteins_id[i:min(i+batch_size, proteins_num)]
        # Add list of protein ids to parameters
        params['accession'] = ','.join(proteins_id_)
        print(params['accession'])
        print()
        # Make request to API
        res = req.get('/'.join([BASE_URL, 'proteins']),
                      headers={'Accept': accept},
                      params=params)
        # Check query result
        if res.status_code != 200: break
        # Concatenate output
        out = out + res.text    
    # Return results
    return res.status_code == 200, out, res
"""
    