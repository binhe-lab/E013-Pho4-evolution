"""
title: run PADDLE in batch mode for AED test constructs
author: Bin He
date: 2024-11-11
updated: 2024-11-11 9:30pm
note: this script is largely adapted from the ipython notebook from the PADDLE authors
note2: updated to work on extra sequences. see below
note3: updated to use [0,12] as ylim
"""

import numpy as np
import pandas as pd
import tensorflow as tf # tested on version 2.2.0
import matplotlib
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['font.sans-serif'] = 'Arial'
from matplotlib import pyplot as plt
from Bio import SeqIO
from datetime import date
import paddle

# Set Numpy to display floats with 3 decimal places
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

########
# Run predictions using PADDLE. This network requires predicted secondary
# structure (from PSIPRED) and predicted disorder (from IUPRED2, in both the
# short and long modes) as input in addition to the protein sequence. This
# model is the most accurate and should be used for predicted ADs in wild-type
# proteins.
########

# Load PADDLE models. There are 10 in total, and their individual predictions
# are averaged to obtain the final result.
pad = paddle.PADDLE()

# Specify the fasta file
fasta_file = 'data/fasta/20241111-AED-test-Y1H-extra-constructs.fa'

# Date prefix
d = date.today().strftime("%Y%m%d") # add a date prefix to the filename

# Iterate through the AED construct sequences
for record in SeqIO.parse(fasta_file, "fasta"):
    # parse name and sequence from fasta record
    sequence_name = record.id
    prot = str(record.seq)
    
    # set up the input and output file names
    ss_file = 'data/s4pred-local/' + sequence_name + '.ss2'
    dis_long_file = 'data/iupred3/' + sequence_name + '-long.dis'
    dis_short_file = 'data/iupred3/' + sequence_name + '-short.dis'
    output_file = 'output/' + d + '-' + sequence_name + '-PADDLE-prediction.csv'
    img_file = 'output/img/' + d + '-' + sequence_name + '-PADDLE-z-score-0-12.png'    
    
    # Load pre-computed secondary structure predicted by s4pred
    # this is run without using BLAST, this speeds up secondary structure prediction.
    seq, helix, coil = paddle.read_ss2(ss_file)
    assert prot == seq

    # Load pre-computed disorder predicted by IUPRED2, in both
    # the short and long modes.
    prot, dis_short = paddle.read_iupred(dis_short_file)
    assert prot == seq
    prot, dis_long = paddle.read_iupred(dis_long_file)
    assert prot == seq

    # Run predictions on all 53 amino acid long tiles across the protein.
    # This function requires matching protein sequence and secondary structure scores.
    # Returns a Numpy array of size (protein_length-52) which gives the
    # predicted activation Z-score for the 53aa tiles starting at positions
    # 1, 2, 3, ..., protein_length-52.
    # High-strength ADs can be called by finding >=5 consecutive positions with Z-score > 6.
    # Medium-strength ADs can be called by finding >=5 consecutive positions with Z-score > 4.
    preds = pad.predict_protein(prot, helix, coil, dis_short, dis_long)
    
    # Plot predicted Z-scores based on the tiles' position in the protein.
    plt.figure(figsize=(8,2.5))
    plt.plot(np.arange(len(preds)) + (53+1)/2, preds)
    plt.xlabel('Center position of 53aa tile')
    plt.ylabel('Predicted Z-score')
    plt.title(sequence_name)
    plt.xlim(0,350)
    plt.ylim(0,12)
    plt.savefig(img_file, transparent = True, dpi = 300)
    plt.show()

    # Z-scores can be converted to fold-activation. This is only useful
    # as a rough reference for activation in S. cerevisiae, as the
    # fold-activation will vary between different experimental conditions.
    pred_act = paddle.score2act(preds)
    
    #plt.figure(figsize=(8,2.5))
    #plt.semilogy(np.arange(len(pred_act)) + (53+1)/2, pred_act)
    #plt.xlabel('Center position of 53aa tile')
    #plt.ylabel('Predicted activation')
    #plt.title(name + 'Pho4')
    #plt.xlim(0,len(prot))
    #plt.show()

    # Export both the Z-score and the fold-activation data
    pred_df = pd.DataFrame({'Pos': np.arange(len(pred_act)) + (53+1)/2, 'Z_score': preds, 'Fold': pred_act})
    pred_df.to_csv(output_file, index = False)
