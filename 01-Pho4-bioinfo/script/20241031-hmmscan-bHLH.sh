# title: use hmmscan to identify the region of bHLH in Pho4 orthologs
# author: Bin He
# 
hmmscan -o 20241031-hmmscan.log --domtblout ../output/20241031-Pho4-ortholog-bHLH-hmmscan.txt --domE 1e-3 ../input/domain-hmm-profile/PF00010.hmm ../input/Pho4-alignment/20241031-Pho4-16sps.aa.fa
