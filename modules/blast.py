################################
### Handle BLAST web service ###
################################

# Dependencies
import requests as req
import pandas as pd
import json
import xmltodict


# Constants
BASE_URL = r'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast'


# Parse xml job result
def parse_xml(job_result):
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


# Parse tabular result
def parse_tbl(job_result, sep=','):
    return pd.DataFrame(
        # Set data, skip comments and empty rows
        data=[row.split(sep) for row in job_result.split('\n') if row != '' and row[0] != '#'],
        # Define column names
        columns={
            'query_ac': str(), 'subj_ac': str(), 'perc_identity': float(),
            'align_len': int(), 'mismatches': int(), 'gap_opens': int(),
            'q_start': int(), 'q_end': int(), 's_start': int(), 's_end':int(),
            'e_value': float(), 'bit_score': float()
        })


# Make parameters request
def get_parameters():
    # Retrieve web service results
    res = req.get('/'.join([BASE_URL, 'parameters']),
                  headers={'Content-type': 'application/json',
                           'Accept': 'application/json'})
    # Return tuple (is status code '200 OK'?, result content, result object)
    return res.status_code == 200, res.json().get('parameters'), res


# Submit job and retrieve its id
def run_job(email, sequence, program='blastp', matrix='BLOSUM62',
            alignments=1000, scores=1000, evalue='1e-3', filter=True,
            seqrange='START-END', gapalign=True, align=6, stype='protein',
            database='uniprotkb', params={}):
    """
    Run BLAST job on EBI web service, then retrieve the id of the started job
    Official documentation for BLAST web api on EBI: https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?pageId=68167377
    Input:
        1.  email:      a valid email address, used to track the user
        2.  sequence:   the actual sequence on which the BLAST must be run
        3.  program:    blast program to be used
        4.  matrix:     the score matrix used to run BLAST
        5.  alignments: maximum number of alignments to retrieve
        6.  scores:     maximum number of match score summaries reported in the result output
        7.  evalue:     evalue threshold (as string)
        8.  filter:     filter low complexity regions
        9.  seqrange:   specify a range or section of the input sequence to use in the search
        10. gapalign:   specify wether to align gaps or not
        11. align:      output format (default 10, tabular)
        12. stype:      sequence type (protein|RNA|DNA)
        13. database:   database on which BLAST must be run
        14. params:     other params for the request
    Output:
        1. Is the status of the response? 1|0
        2. Response text
        3. Response object (useful only for debug)
    """
    # Merge user-specified and default parameters
    params = {**params, **{
        # Mandatory parameters
        'email':        email,
        'sequence':     sequence,
        'program':      program,
        'matrix':       matrix,
        'alignments':   alignments,
        'scores':       scores,
        'exp':          evalue,
        'filter':       'T' if filter else 'F',
        'seqrange':     seqrange,
        'gapalign':     gapalign,
        'align':        align,
        'stype':        stype,
        'database':     database,
    }}
    # Make request
    res = req.post('/'.join([BASE_URL, 'run']),
                   headers={
                    'Content-type': 'application/x-www-form-urlencoded',
                    'Accept': 'text/plain'
                   },
                   data=params)
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
def get_job_result(job_id, result_type='out', result_parse=lambda x: parse_tbl(x, sep='\t')):
    # Make request
    res = req.get('/'.join([BASE_URL, 'result', job_id, result_type]))
    # Return results
    return res.status_code == 200, result_parse(res.text), res
