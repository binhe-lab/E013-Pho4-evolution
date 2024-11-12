"""
title: batch prediction using iupred3 for multiple sequence fasta
date: 2024-11-11
comment:
    this script is written with help from ChatGPT, modified by Bin He
    it assumes Biopython is installed
"""

import os
import sys
import subprocess
from Bio import SeqIO

# Check if the user provided a file name as an argument
if len(sys.argv) != 3:
    print("Usage: python 20241111-run-iupred3-batch.py <fasta_file> <outdir>")
    sys.exit(1)

# Paths to the input FASTA file and the `iupred3.py` script
fasta_file = sys.argv[1]
output_dir = sys.argv[2]
iupred_script = "./iupred3.py"  # Replace with the path to the `iupred3.py` script
temp_file = "./tmp.fasta"        # temporary file used as input to `iupred3.py`

# Create a directory to store the output files
os.makedirs(output_dir, exist_ok=True)

# Modes to run
modes = ['long', 'short']

# Iterate through each sequence in the FASTA file
for record in SeqIO.parse(fasta_file, "fasta"):
    sequence_name = record.id
    sequence = str(record.seq)
    # write the sequence to a temporary file as input for iupred3.py
    with open(temp_file, 'w') as tempfile:
        tempfile.write(f">{sequence_name}\n{sequence}\n")

    for mode in modes:
        # Prepare output file path
        output_file = os.path.join(output_dir, f"{sequence_name}-{mode}.dis")

        # Run the `iupred3.py` script for each sequence in each mode and save the output
        with open(output_file, "w") as outfile:
            process = subprocess.run(
                ["python", iupred_script, temp_file, mode],
                stdout=outfile,
                stderr=subprocess.PIPE
            )

        # Check for errors
        if process.returncode != 0:
            print(f"Error processing sequence {sequence_name}: {process.stderr.decode()}")

print("Processing complete. Results are saved in:", output_dir)

