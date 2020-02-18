#################################
### Handle UniProt search API ###
#################################

# Dependencies
import requests as req
import time
import json
import tempfile
import xmltodict
import gzip


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


# Define function for retrieving mapping entries batches
def map_ids(ids, from_db='NF90', to_db='NF90', out_format='fasta', batch_size=100, params={}):
    """
    Maps batches of entries from <from_db> to <to_db>, formatted as <out_format> format
    Note that GET action cannot handle too many data, hence query must be splitted in batches
    More information about this api can be found at https://www.uniprot.org/help/api_idmapping
    Input:
    1. ids:     batch of some database IDs (e.g. Pfam, UniProtKB, ...)
    2. format:  output format
    3. params:  extra parameters to be used in query
    Output:
    1. Mapping <from> -> <to>, retrieved as specified by <format>
    """
    # Define output container
    out = ''
    # Define number of proteins id
    num_ids = len(ids)
    # Loop through each batch
    for i in range(0, num_ids, batch_size):
        # Define ids batch
        batch_ids = ids[i:min((i+batch_size), num_ids)]
        # Make API request
        res = req.get('/'.join([BASE_URL, 'uploadlists']),
                      params={**params, **{
                        'from': from_db,
                        'to': to_db,
                        'format': out_format,
                        'query': ' '.join(batch_ids)
                      }})
        # Check result
        if res.status_code != 200: break
        # Store results
        out = out + res.text
    # Return results
    return res.status_code == 200, out, res


# Define function for making a generic query
def make_query(query, params={}):
    """
    Executes a generic query in UniProt, by using official APIs.
    UniProt's official query API info: https://www.uniprot.org/help/api_queries
    Retrievable data columns: https://www.uniprot.org/help/uniprotkb_column_names
    More in depth documentation for database(...) parameter: https://www.uniprot.org/docs/dbxref
    Input:
        1. query:       search query, formatted with UniProt Standards
        2. params:      other search parameters, which are concatenated after query
    Output:
        1. status:      did the request return a 200 OK status code? 1|0
        2. result:      result of the query (reliable only if status==1)
        3. response:    response object, useful for debug purposes
    """

    # Define params
    params = {**{
        'sort': 'score',
        'format': 'tab',
        'compress': 'no',
        'columns': 'id'
    }, **params}
    # Add query
    params.setdefault('query', query)
    # Make query
    response = req.get('/'.join([BASE_URL, 'uniprot']),
                        headers={'Accept': 'text/plain'},
                        params=params)
    # Get result
    result = response.text
    status = (response.status_code == 200)
    # If compressed: uncompress
    if params.get('compress', 'no') == 'yes' and status:
        result = str(gzip.decompress(response.content), 'utf-8')
    # Get result
    return status, result, response


# Unit testing
if __name__ == '__main__':

    # Test multiple sequences retrieval trhough fasta
    status, result, response = map_ids(['UniRef90_P43582', 'UniRef90_J8Q8J2', 'UniRef90_J5RH20'])
    # Check status
    assert status, 'Error while retrieving fasta result'
    # Check fasta format
    assert sum([1 for c in result if c == '>']) == 3, 'Error: format is not fasta'

    # Test query results
    status, result, response = make_query(query='insulin', params={
        'compress': 'yes',
        'columns': ','.join(['id', 'entry name', 'reviewed', 'protein names', 'genes', 'organism', 'length']),
        'sort': 'score',
        'format': 'tab'
    })
    # Check status
    assert status, 'Error while retireving query result'
    # Get result header
    rows = result.split('\n')
    header = rows[0].split('\t')
    assert header[0] == 'Entry', 'Error, first header item does not match the expected one'
    assert header[-1] == 'Length', 'Error, last header item does not match the expected one'
