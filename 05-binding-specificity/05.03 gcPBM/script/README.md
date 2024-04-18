# Script file expanation
`normalize_agilent_array_alexa488_deBruijn_plus_custom.pl` script for processing the scanned images. modified by Yuning from Raluca's lab, 2024-04-17
`modules/masliner_fmedian.txt`: for combining multiple scans using linear regression. see `PBM_analysis_details.doc` in the `../../05.02 universal PBM/doc` folder.

# Data preprocessing notes
2024-04-18, HB

Following the Analysis_details documentation above, and also based on the recommendation of Wei (see `../input/20230121-gcPBM-scan-raw/README.md`, we will

1. Use the high concentrations for both Pho4, namely chamber 2-4 for CgPho4 (2uM) and 4-4 for ScPho4 (2 uM).
1. For chamber 2-4, we will combine the "550_30" and "550_10" scans.
1. For chamber 4-4, we will use "550_5" scan.

## Combine scans using Masliner
Based on the PBM_analysis_details documentation above, we used the following command
```unix
perl modules/masliner_fmedian.txt \
  -g1 ../input/20230121-gcPBM-scan-raw/20230121_cgPho4_yPho4_550_10_2-4.gpr \
  -g2 ../input/20230121-gcPBM-scan-raw/20230121_cgPho4_yPho4_550_30_2-4.gpr \
  -o ../input/masliner-processed/20230121_cgPho4_yPho4_550_10+30_2-4.gpr \
  -ll 2000 -lh 50000 -f1 488 -f2 488
```

## Detrending
For ScPho4, using chamber 4-4, with 550_5 scan
```unix
perl normalize_agilent_array_alexa488_deBruijn_plus_custom.pl \
  -i ../input/20230121-gcPBM-scan-raw/20230121_cgPho4_yPho4_550_5_4-4.gpr \
  -s ../input/20230125_8x60k_pho4_orthologos_analysis.txt \
  -o ../input/normalized-intensities/20240418-ScHigh
```

Output
```
Array format recognized as '8x60K'.
Reading Input File "../input/20230121-gcPBM-scan-raw/20230121_cgPho4_yPho4_550_5_4-4.gpr".
NOT normalizing Alexa488 Intensities by Cy3.
Overlall Alexa488_median: 1399
Calculating Neighborhood Median for Alexa488. Radius = 7.
Printing to output files for "../input/normalized-intensities/20240418-ScHigh".
```

For CgPho4, using chamber 2-4, with the 550_30 and 550_10 combined output
```unix
perl normalize_agilent_array_alexa488_deBruijn_plus_custom.pl \
  -i ../input/masliner-processed/20230121_cgPho4_yPho4_550_10+30_2-4.gpr \
  -s ../input/20230125_8x60k_pho4_orthologos_analysis.txt \
  -o ../input/normalized-intensities/20240418-CgHigh
```

Output:
```
Array format recognized as '8x60K'.
Reading Input File "../input/masliner-processed/20230121_cgPho4_yPho4_550_10+30_2-4.gpr".
NOT normalizing Alexa488 Intensities by Cy3.
Overlall Alexa488_median: 1505.5
Calculating Neighborhood Median for Alexa488. Radius = 7.
```

## Remove ^M carriage returns
The output files contain "^M" that break some of the lines. To remove them, use the following command

`tr -d '^M' < input > output`: to type ^M, actually type Ctrl_V Ctrl_M. replace "input" and "output" with the input and output file names
