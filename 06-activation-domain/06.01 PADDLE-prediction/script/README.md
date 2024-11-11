this folder contains auxilary script files. the main scripts for running PADDLE is in the parent folder

## 2024-03-15, install `s4pred`
followed [instructions](https://github.com/psipred/s4pred?tab=readme-ov-file) and installed the software in my local computer (not inside the git repo)

to make predictions, I copied the fasta file to the `s4pred` folder and used the following line to run the prediction

```unix
python run_model.py -s --save-files --outdir ./tmp/ ./20201121-merge-pho4-orthologs-aa.fa
```

The result files were saved in the `tmp` folder
