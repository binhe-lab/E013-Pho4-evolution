# Goal
Use PADDLE (Sanborn et al 2021 eLife) to predict activation potential in ScPho4 and CgPho4.

# Install
The PADDLE scripts and models were downloaded from <https://github.com/asanborn/PADDLE> on 2024-03-05

To install `tensorflow` on my M1 mac, I used the following lines

```unix
conda create -n E013 python=3.10.0 tensorflow
conda deactivate; conda activate E013
mamba install -c conda-forge jupyterlab
pip install matplotlib pandas
```

# Analysis
## 2024-03-05
The first batch of predictions were made using the ipython notebook provided by the author in the github repo. `PSIPRED` was run using the [PSIPRED Workbench](bioinf.cs.ucl.ac.uk/psipred/) server. `IUPRED2` disorder predictions were made using the [IUPRED server](https://iupred2.elte.hu/download_new). These results were then manually coded into the ipython notebook, adding one cell per sequence. Results were saved as CSV and PNG files.

## 2024-11-11
Here, the goal is to respond to reviewer's question, i.e., whether the Y1H constructs designed to test the AEDs may have inadvertantly created new ADs by synthetically fusing regions.

To automate the analysis, I followed the following steps:
- Use `s4pred` to predict the secondary structure locally, without BLAST
- Use `iupred3` to predict disorder locally. Wrote a wrapper script to process multiple sequences in a fasta.
- Wrote a wrapper script to run PADDLE on a multiple fasta input.

I used `s4pred` and `iupred3` locally to make the secondary structure and disorder predictions for the construct sequences (`data/fasta/20241104-AED-test-Y1H-constructs.fa` and `data/fasta/20241111-AED-test-Y1H-extra-constructs.fa`; the latter contains several additional sequences Emily is working on). Details can be found in the respective `script` and `data` folders. I then wrote a wrapper script to run PADDLE on all the sequences (`20241111-run-PADDLE-AED.py`). I ran the script as follows

```unix
python3 20241111-run-PADDLE-AED.py >20241111-run-PADDLE-AED.log 2>20241111-run-PADDLE-AED.err
```
