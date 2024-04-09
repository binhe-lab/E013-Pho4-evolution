Below are the README file written by Jinye. The microscopy images are stored on RDSS. The output from the ImageJ analysis are stored in this folder for quantitative analysis. The sub-folder contains JY's processed data, which we will not use for analysis. Note that in some of the analyses, she subtracted the autofluorescence background using the average of five nonfluorescent cells. Our current analysis uses data before the background subtraction, because 1) not all replicates have the non-fluorescent cells measured in the same experiment; 2) subtraction of autofluorescence results in negative numbers for some cells. Neither are serious issues. I (Bin) plan to analyze the subtracted data to confirm the results.

A note on the data file headers:
- Pair_No: id for cell nested within "Group"
- Area: area in pixels
- Mean: mean fluorescence intensity, = RawIntDen / Area
- IntDen: Mean * area in scaled units, ignore for this analysis
- Median: median fluorescence intensity
- RawIntDen: sum of pixel values in the selected area
- Group: genotype
- Locus: area selected

Data collector&analyzer: Jinye Liang
Archived date: 20240228

Experiment setup: Pho4-mNeON localization in mid-log phase S.cerevisiae strains grown in SC complete medium; Leica epi-fluorescence microscope (Bright Field, GFP, DAPI)

Protocol: Manually select the DAPI stained region to indicate nucleus localization, see https://docs.google.com/document/d/1gsf_96zyUpWXnpNM9d6yda8O_A2H_LCB0lZvZ_rq_lI/edit?usp=sharing

Folder description: two biological replicates with Pho2: d013024 and d020624; two biological replicates w/o Pho2: d021424 and d021624

Raw image & Leica setting folder localization: rdss_bhe2\jinye-liang\data\microscopy (4 folders based on the date)

Description: This folder includes the raw image, Fiji parameters for whole cell (wc), cytoplasm and nuclei (nc), Fiji parameter output files(Mean and Median gray value, area, Integrated Density = mean gray value * area and etc.), ratio data calculated based on the output files and R Markdown files for analyzing the data. plot_nuc_ratio is the R repository.
