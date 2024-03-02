""" merge Sanger sequences
Input: fasta file, whether sequences are forward or reverse
Output: assembled sequences

Bin He
6 september 2019
"""

import re
from datetime import datetime
from Bio import SeqIO

# open file for storing the output
out = open("../output/"+datetime.now().strftime("%Y%m%d")+"-motif-count.txt","w")

# open file for read only
for record in SeqIO.parse("../data/20190628-promoter-1500bp-upstream.fasta", "fasta"):
    # extract record id information
    desc = record.description
    geneID = re.search(r"FEATURE NAME=(.*)\]", desc).group(1)
    geneName = re.search(r"for gene: (.*);", desc).group(1)
    # count pattern matches
    seq  = str(record.seq)
    matches = re.finditer(pat, seq)
    for m in matches:
        out.write(geneID+"\t"+geneName+"\t"+str(m.start()-len(seq))+"\t"+m.group()+"\n")
