the dated files are run and downloaded from [IUPred2A](https://iupred2a.elte.hu/) server. The undated files are generated using `awk`, separating the two merged files into individual sequence results. the command used is below:

```unix
# for the long disorder prediction
grep -v '^#' 20240316-merge-pho4-orthologs-aa-long.result | awk 'BEGIN {FS="\t"; OFS="\t"} /^>/ {file=substr($1, 2, 4) "Pho4-long.dis"; next} { print $1,$2,$3 > file}'
# for the short disorder prediction
grep -v '^#' 20240316-merge-pho4-orthologs-aa-short.result | awk 'BEGIN {FS="\t"; OFS="\t"} /^>/ {file=substr($1, 2, 4) "Pho4-short.dis"; next} { print $1,$2,$3 > file}'
```
