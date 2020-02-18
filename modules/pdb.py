###########################
### PDB FILES RETRIEVAL ###
###########################


# Dependencies
import requests
import gzip
import re
from os.path import isfile


# Constants
BASE_URL = r'http://files.rcsb.org'  # Path to PDB
# http://files.rcsb.org/download/1eg3.pdb.gz

# Retrieve PDB file, given a pdb_id
def download(pdb_id, format='pdb', out_path=None, compressed=True):
    # Define path to pdb
    pdb_file = '.'.join([pdb_id, format] + (['gz'] if compressed else []))
    # Redefine output path (if directory, add output file)
    if not isfile(out_path):
        out_path = '/'.join([re.sub(r'[/]$', '', out_path), '.'.join([pdb_id, format])])
    # Make download rquest
    response = requests.get('/'.join([BASE_URL, 'download', pdb_file]))
    # Check response status
    if response.status_code == 200:
        # Get file content (as bytes)
        content = response.content
        # If content is compressed: decompress it
        content = gzip.decompress(content) if compressed else content
        # Write to file
        with open(out_path, 'wb') as out_file:
            # Save uncompressed file
            out_file.write(content)
    # Return if response statusis 200 OK, downloaded file path, response
    return response.status_code == 200, out_path, response


# Unit testing
if __name__ == '__main__':

    # Retrieve a known file
    download('1eg3', out_path='data/', format='pdb', compressed=True)
