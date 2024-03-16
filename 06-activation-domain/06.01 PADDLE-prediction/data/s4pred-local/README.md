following [S4PRED instruction](https://github.com/psipred/s4pred?tab=readme-ov-file) and ran the analyses locally using the following command

```unix
python run_model.py --save-files --outdir ./tmp/ --save-by-idx ./20201121-merge-pho4-orthologs-aa.fa
```

installed `biopython` and `pytorch` using
`mamba install biopython=1.78 pytorch torchvision -c pytorch`
