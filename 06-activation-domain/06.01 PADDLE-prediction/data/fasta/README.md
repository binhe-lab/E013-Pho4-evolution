the individual fasta sequence files were "split" from the "merge" file using `awk`

`awk '/^>/ {file="tmp/" substr($1, 2, 4) "Pho4.fa"} { print > file }' 20201121-merge-pho4-orthologs-aa.fa'
