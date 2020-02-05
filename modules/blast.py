################################
### Handle BLAST web service ###
################################

# Dependencies
import requests as req
import json
import xmltodict


# Constants
BASE_URL = r'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast'


# Parse xml jpb result
def parse_xml_job_result(job_result):
    # Parse xml to ordered dictionary
    job_result = xmltodict.parse(job_result)
    # Retrieve list of hits
    hits_in = list(job_result['EBIApplicationResult']['SequenceSimilaritySearchResult']['hits']['hit'])
    hits_out = list() # New list with simplifiued hits objects
    # Turn each hit into a simple vocabulary
    for i in range(len(hits_in)):
        # Get hit attributes
        hit_attr = {'database': str(hits_in[i]['@database']),
                    'id': str(hits_in[i]['@id']), 
                    'ac': str(hits_in[i]['@ac']), 
                    'description': str(hits_in[i]['@description'])}
        # Get hit content
        hits_in[i] = hits_in[i]['alignments']['alignment']
        # Check if hit is a list, if not, turn it into a list of one
        if type(hits_in[i]) != list:
            # Create a list of one
            hits_in[i] = [hits_in[i]]
        # Loop through each hit
        for j in range(len(hits_in[i])):
            # Get hit object
            hit = hits_in[i][j]
            # Simplify hit vocabulary: numeric values
            hits_out.append({**{'score': int(hit['score']),
                                'bits': float(hit['bits']),
                                'expectation': float(hit['expectation']),
                                'identity': float(hit['identity']),
                                'positives': float(hit['positives']),
                                'gaps': int(hit['gaps']),
                                'strand': str(hit['strand']),
                                'pattern_seq': str(hit['pattern']),
                                'match_seq': str(hit['matchSeq']['#text']),
                                'match_start': int(hit['matchSeq']['@start']),
                                'match_end': int(hit['matchSeq']['@end'])}, 
                             # Add hit attributes
                             **hit_attr})
    # Return hits
    return hits_out

        
# Make parameters request
def get_parameters():
    # Retrieve web service results
    res = req.get('/'.join([BASE_URL, 'parameters']), 
                  headers={'Content-type': 'application/json',
                           'Accept': 'application/json'})
    # Return tuple (is status code '200 OK'?, result content, result object)
    return res.status_code == 200, res.json().get('parameters'), res


# Submit job and retrieve its id
def run_job(email, sequence, params={}):
    # Merge user-specified and default parameters
    params = {**{'program': 'blastp',
                 'stype': 'protein',
                 'database': 'uniprotkb'},
              **params}
    # Update mandatory parameters
    params['email'] = email
    params['sequence'] = sequence
    # Make request
    res = req.post('/'.join([BASE_URL, 'run']), 
                           headers={'Content-type': 'application/x-www-form-urlencoded',
                                    'Accept': 'text/plain'}, 
                           params=params)
    # Return result
    return res.status_code == 200, res.text, res


# Retrieve job status
def get_job_status(job_id):
    # Make request
    res = req.get('/'.join([BASE_URL, 'status', job_id]),
                  headers={'Accept': 'text/plain'})
    # Return results
    return res.status_code == 200, res.text, res


# Retrieve job result
def get_job_result(job_id, result_type='xml', result_parse=parse_xml_job_result):
    # Make request
    res = req.get('/'.join([BASE_URL, 'result', job_id, result_type]),
                  headers={'Accept': 'application/xml'})
    # Return results
    return res.status_code == 200, result_parse(res.text), res