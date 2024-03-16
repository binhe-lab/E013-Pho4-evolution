"""
this code is written largely with the help of ChatGPT
it's not working yet, as the http request continues giving 500 or 404 errors
switched to s4pred for local predictions
author: Bin He
date: 2024-03-15
"""

import sys
import os
import requests
import pdb

def psipred_prediction(name, sequence):
    url = "http://bioinf.cs.ucl.ac.uk/psipred/api/submission/"
    data = {
        "job": "psipred",
        "submission_name": name,
        "email": 'bin-he@uiowa.edu',
        "input_data": sequence
    }

    response = requests.post(url, data=data)
    if response.status_code == 200:
        job_id = response.json()["id"]
        return job_id
    else:
        print(f"Error submitting sequence: {response.text}")
        return None

def retrieve_psipred_result(job_id):
    url = f"http://bioinf.cs.ucl.ac.uk/psipred/api/result/{job_id}"
    response = requests.get(url)
    if response.status_code == 200:
        result = response.json()
        return result
    else:
        print(f"Error retrieving result for job {job_id}: {response.text}")
        return None

def main(filename):
    pdb.set_trace()
    sequence_name = os.path.splitext(os.path.basename(filename))[0]
    with open(filename, 'r') as file:
        sequence = ''.join(line.strip() for line in file.readlines() if not line.startswith('>'))
    print(f"Submitting sequence {sequence_name}...")
    job_id = psipred_prediction(sequence_name, sequence)
    if job_id:
        print(f"Job ID for sequence {sequence_name}: {job_id}")
        print(f"Waiting for result for sequence {sequence_name}...")
        result = None
        while not result:
            result = retrieve_psipred_result(job_id)
        output_filename = f"../data/psipred/{sequence_name}.ss"
        with open(output_filename, 'w') as output_file:
            output_file.write(result['ss'])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 20240315-psiprd-query.py <filename>")
        sys.exit(1)
    filename = sys.argv[1]
    main(filename)
