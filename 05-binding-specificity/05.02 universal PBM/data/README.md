this folder contains the intermediate results from the uPBM experiments. See `../docs/PBM_analysis_details.doc` for details. an excerpt from that tutorial explains the files in this folder:

> The program outputs four files: 
> (1) a “raw data” file, containing a summary of the un-normalized probe intensities, flags, and sequences extracted from the GPR files; 
> (2) an “all data” file, containing all the information above plus the calculated values for expected Cy3, observed/expected Cy3, Cy3-normalized Alexa488, and spatially-detrended Alexa488; 
> (3) a “deBruijn” file, containing a ranked list of the normalized intensities and sequences of all combinatorial ‘all k-mer’ probes (with control spots removed); and 
> (4) a “regression” file, containing the coefficients determined from the linear regression over Cy3 probe intensities and sequences, including the R2 value indicating the quality of the fit.  The most important of these is the “deBruijn” file, which can be input directly to the sequence analysis program (described below) or used by other motif finding algorithms.  The other output files are mainly intended for users more familiar with microarray analysis who wish to apply additional normalization or filtering steps.
