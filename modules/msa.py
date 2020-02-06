#########################################################
### Handle Multiple Sequence Alignment using EBI apis ###
#########################################################


# Dependencies
import requests as req
import time


# Constants
BASE_URL = r'https://www.ebi.ac.uk/Tools/services/rest'
CLUSTALO = 'clustalo'
MUSCLE = 'muscle'


# Send multiple sequence alignment request to EBI apis
def run_job(email, sequence, algorithm=CLUSTALO, params={}):
    # Set mandatory parameters
    params['email'] = email
    params['sequence'] = sequence
    # Make api request to clustal omega endpoint
    response = req.post('/'.join([BASE_URL, algorithm, 'run']), headers={
        'Content-type': 'application/x-www-form-urlencoded',
        'Accept': 'text/plain',
    }, params=params)
    # Retrieve result
    return response.status_code == 200, response.text, response


# Retrieve job status
def get_job_status(job_id, algorithm=CLUSTALO):
    # Make request
    response = req.get('/'.join([BASE_URL, algorithm, 'status', job_id]),
                  headers={'Accept': 'text/plain'})
    # Return results
    return response.status_code == 200, response.text, response


# Retrieve job result
def get_job_result(job_id, algorithm=CLUSTALO, result_type='aln-fasta', result_parse=lambda x: x):
    # Make request
    response = req.get('/'.join([BASE_URL, algorithm, 'result', job_id, result_type]),
                  headers={'Accept': 'text/x-fasta'})
    # Return results
    return response.status_code == 200, result_parse(response.text), response


# Wrapper for Clustal Omega multiple sequence alignment algorithm
def run_clustalo(email, sequence, params={}):
    return run_job(email=email, sequence=sequence, algorithm=CLUSTALO, params={**{
        # Default parameters
        'outfmt': 'fa',
        'order': 'input',
        'stype': 'protein'
    }, **params})


# Wrapper for MUSCLE omega multiple sequence alignment algorithm
def run_muscle(email, sequence, params={}):
    return run_job(email=email, sequence=sequence, algorithm=MUSLCE, params={**{
        # Default parameters
        'format': 'fasta'
    }, **params})


# Test
if __name__ == '__main__':
    # Test clustal omega multiple sequence alignment
    status, job_id, response = run_clustalo(
        email='damiano.clementel@studenti.unipd.it',
        sequence=""">UniRef90_P43582 WW domain-containing protein WWM1 n=9 Tax=Saccharomyces TaxID=4930 RepID=WWM1_YEAST
MAQSKSNPPQVPSGWKAVFDDEYQTWYYVDLSTNSSQWEPPRGTTWPRPKGPPPGVNNEK
SSRQQADQAPPPYSSQSTPQVQAGAQAQQPRYYQPQQPQYPQYPQQQRYYPQQAPMPAAA
PQQAYYGTAPSTSKGSGHGGAMMGGLLGVGAGLLGGAMLEHAFDDHNYDGPDTVVVENNY
YGDDAGGSDGGFDDAGGFDGGFDDGFDGSDF
>UniRef90_J8Q8J2 Wwm1p n=1 Tax=Saccharomyces arboricola (strain H-6 / AS 2.3317 / CBS 10644) TaxID=1160507 RepID=J8Q8J2_SACAR
MAQSKSNPPQVPSGWKAVFDDEYQTWFYVDLSTNNSQWEPPKGASFPRPKGPPPAANNEK
TSRQQGDQAPPPYSAQSRTQPQPQAQQAQQGRYYQPQQPQYPQQPQQQSYYPQQVPMAAA
AAPQQGYYGATPTAAKSSGRSGAMMGGLLGVGAGLLGGAMLEHAFDDHSHGGPGPVVENN
YYGDDNGGGFGGPFDGGFDGEFDGGFDGGDF
>UniRef90_J5RH20 WWM1-like protein n=2 Tax=Saccharomyces TaxID=4930 RepID=J5RH20_SACK1
MAQSKGNPPQVPSGWKAVFDDEYQTWFYVNLSTNSSQWEPPKGTTWPRPKGPPPGVNNEK
SSRQEVDQAPPPYSSQSRAQPQAPAQQTRYYQPQQSQYPQQPQQQRYYQQQAPMAAAAPQ
QAYYGTTPSAAKSSGHGGAMMGGLLGVGAGLLGGAMLEHAFDDHDNYDQPDNVVVENNYY
GDDGGFDGGFDGGFDGGDF""",
        params={}
    )
    # Check status
    assert status, 'Error: job not started'

    # Test status retrieval
    while True:
        # Check response status
        status, job_status, response = get_job_status(job_id, algorithm='clustalo')
        # Check if response status is 200 OK
        assert status, 'Error: cannot retrieve job status'
        # Check if job has finished, then go on
        if job_status == 'FINISHED': break
        # Add delay
        time.sleep(3)

    # Check job result retrieval (fasta)
    status, job_result, response = get_job_result(job_id, algorithm='clustalo')
    # Check response status
    assert status, 'Error: cannot retrieve job result'
    # Check correctedness of job result
    assert sum([1 for c in job_result if c == '>']) == 3, 'Error: retrieverd fasta file is not coherent'
