more details can be found in `../script/README.md`

following [S4PRED instruction](https://github.com/psipred/s4pred?tab=readme-ov-file) and ran the analyses locally.

installed `biopython` and `pytorch` using
`mamba install biopython=1.78 pytorch torchvision -c pytorch`

## 2024-03-15, Pho4 orthologs

```unix
python run_model.py --save-files --outdir ./tmp/ --save-by-idx ./20201121-merge-pho4-orthologs-aa.fa
```

## 2024-11-11, AED constructs

```unix
python run_model.py --save-files --outdir ./20241111-AED-constructs --save-by-idx ./20241104-AED-test-Y1H-constructs.fa
```
