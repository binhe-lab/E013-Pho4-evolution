---
title: "Count CgPho4 binding motif flanking nucleotide frequency"
output: html_notebook
author: Bin He
date: 17 mai 2017
---

input: tab-delimited spreadsheet containing motif and flanking sequences under the CgPho4 ChIP peaks

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 

## Background & Goal

I want to know the nucleotide (A/T/C/G) frequencies at the flanking nucleotides surrounding the E-box motif. 

In particular, we know that for ScPho4, the 5' T-CACGTG motifs are bound at much lower rate than the A/C/G-CACGTG ones, because Cbf1, another bHLH TF, binds to the 5' T version with high affinity. For CgPho4, we know that it binds to nearly all ScPho4 bound sites, and 50% more. The question is, do CgPho4-specific sites have a different sequence bias compared to the ScPho4-bound ones? 

To answer this question, I first prepared a table with each row corresponding to an E-box motif underneath a CgPho4 ChIP peak. In that row, I included the DNA sequences of the E-box motif "CACGTG" (invariable) as well as the 5bp flanking sequences on each side.

Now the goal is to parse that spreadsheet and output the count matrix at each flanking nucleotide position (we name them after Sebastian Maerkl's 2007 paper, i.e. position -3 ~ 3 represent the core motif, and the 5' flanking positions are called -4, -5, etc. and vice versa)

## Analysis 


```r
library(data.table)
library(GenomicRanges)
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist, unsplit
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:data.table':
## 
##     shift
```

```
## Loading required package: GenomeInfoDb
```


```r
myFreq <- function( seq, a=6, b=11 ){
  # goal: count frequencies at flanking nucleotides
  # input: a vector of sequences pre-aligned, with the 6bp E-box in the middle and 5 flanking nucleotides on each side
  # parameter: a indicates the position in the sequence that corresponds to the beginning of the E-box, while b for the end
  # output: a matrix with four rows and ten columns. Each row represents one nucleotide (a/t/c/g) and each column for one of the ten flanking nucleotides
  split <- strsplit(seq, split = "") # a list of vectors of single letters
  factor <- factor( unlist(split), levels = c("a","t","c","g") )
  mat <- matrix( as.numeric(factor), nrow = length(seq), byrow = T )
  mat <- mat[,-c(a:b)]
  c <- apply( mat, 2, tabulate )
  rownames(c) <- c("a","t","c","g")
  colnames(c) <- as.character(c(-8:-4,4:8))
  freq <- round(c / length(seq) * 100, 0)
  return( freq )
}
```

### CgPho4 binding sites in _S. cerevisiae_


```r
# load data
chip <- fread("./input/20170107_ChIP_stats_in_Scer.txt")
chip <- chip[-23] # remove an overlapping range
motif.seq <- fread("./table_output/20170517_CgPho4_motif_seq_under_peak_in_Scer.txt", skip = "#")

# map columns from chip to motif.seq using GRanges
## build GRanges objects
chip.gr <- with(chip, GRanges( seqnames = CHR, ranges = IRanges( start = START, end = END ) ))

motif.gr <- with(motif.seq, GRanges( seqnames = SEQ_NAME, ranges = IRanges( start = START, end = END ) ))

## find overlaps
map <- findOverlaps( motif.gr, chip.gr, type = "within")

## merge
bound <- chip[subjectHits(map), .(ScPho4.bound, CgPho4.bound)]
motif.tab <- cbind(motif.seq, bound)

## Calculate flanking nucleotides frequency
freq <- list( 
  all.sc4 = myFreq( motif.tab[CgPho4.bound & ScPho4.bound, SEQUENCE] ), 
  all.cg4 = myFreq( motif.tab[,SEQUENCE] ),
  cg4.only = myFreq( motif.tab[CgPho4.bound & !ScPho4.bound, SEQUENCE] )
  )
print( freq )
```

```
## $all.sc4
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 32 18 32 28 16  6 14 25 28 28
## t 31 22 24 33 14 18 35 31 24 26
## c 15 35 15 21 32 44 24 21 22 22
## g 22 25 29 18 38 33 27 24 26 24
## 
## $all.cg4
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 30 24 32 32 13 10 18 20 28 27
## t 30 22 24 26 13 15 30 34 24 30
## c 19 27 17 19 33 47 25 24 23 21
## g 22 27 27 24 42 28 27 21 25 23
## 
## $cg4.only
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 26 34 32 38  6 18 24 12 28 24
## t 28 22 26 14 10 10 20 40 24 36
## c 24 14 20 14 34 52 28 30 24 18
## g 22 30 22 34 50 20 28 18 24 22
```
> CgPho4-specific sites have even lower "A/T" ratio at the -4 and +4 sites compared to sites bound by both ScPho4 and CgPho4.


```r
# first I need to import Xu's mmc2 data
xu.mmc2 <- read.table("./input/20160903_mmc2_Zhou_2011_data_in_Sac03.csv",header = TRUE, as.is = TRUE, sep = "\t")
names(xu.mmc2)[32:33] <- paste(c("Pho4","Pho2"),"Recruitment.no.vs.hiPi",sep = ".")
## remove rows with no ChIP data
xu.mmc2 <- subset(xu.mmc2, Alignability != 0)

# then I need to construct a location set based on the motif locations, expanding each to include the E-box and 5bp flanking sequences on both ends
xu.gr <- with( xu.mmc2, GRanges(seqnames = CHR, ranges = IRanges(Location-2-5,Location+3+5), Bound4 = (Pho4.binding.No.Pi == 1)) ) # Location is set on N~4~, and the range is now expanded to include the whole E-box as well as 5bp flanking nucleotides

## convert the GRanges object to a data.frame and write the result into a text file
## allMotif <- data.frame(SEQ_NAME = seqnames(xu.gr), START = start(xu.gr), END = end(xu.gr)) # convert to data.frame
## write.table(allMotif, "./table_output/20170518_all_CACGTG_motifs_in_Scer.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## In MochiView: import the file generated above and use the export->Locations function to export the sequences based on the locations
## The output is then processed in Vim to make the flanking nucleotides appear as lower cases, to be compatible with myFreq()

# Re-import the sequence set
allMotif.seq <- fread("./input/20170518_all_Ebox_with_5bp_flanking.txt")
## merge into xu.mmc2
allMotif.tab <- data.table(xu.mmc2, SEQUENCE = allMotif.seq$SEQ)
freq$all.motif <- myFreq(allMotif.seq$SEQ)
freq$nucleosome.low.25pc <- myFreq(allMotif.tab[Nucleosome.High.Pi <= quantile(Nucleosome.High.Pi, prob = 0.25), SEQUENCE])
## print result
print(freq)
```

```
## $all.sc4
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 32 18 32 28 16  6 14 25 28 28
## t 31 22 24 33 14 18 35 31 24 26
## c 15 35 15 21 32 44 24 21 22 22
## g 22 25 29 18 38 33 27 24 26 24
## 
## $all.cg4
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 30 24 32 32 13 10 18 20 28 27
## t 30 22 24 26 13 15 30 34 24 30
## c 19 27 17 19 33 47 25 24 23 21
## g 22 27 27 24 42 28 27 21 25 23
## 
## $cg4.only
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 26 34 32 38  6 18 24 12 28 24
## t 28 22 26 14 10 10 20 40 24 36
## c 24 14 20 14 34 52 28 30 24 18
## g 22 30 22 34 50 20 28 18 24 22
## 
## $all.motif
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 31 31 30 32 22 24 33 27 29 30
## t 31 28 29 29 27 23 29 30 29 31
## c 20 19 20 18 26 28 19 20 20 19
## g 19 22 21 21 25 26 19 23 22 20
## 
## $nucleosome.low.25pc
##   -8 -7 -6 -5 -4  4  5  6  7  8
## a 33 30 26 32 25 25 32 28 30 27
## t 34 36 31 26 28 24 32 28 28 32
## c 16 17 18 16 24 20 16 21 19 21
## g 18 17 24 26 22 30 19 24 24 20
```
> 

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file).
