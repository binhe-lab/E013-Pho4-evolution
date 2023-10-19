QC for the flow data using scripts I developed before. The HTML output files are large and not included in the repository. The markdown files should contain some of the output. The main steps take are:

1. a rectangular gate to remove non-cell events.
1. a single cell gate to remove doublets.
1. a clustering gate to identify a homogeneous population.
1. calculate a cell-size normalized fluorescence value (not optimal for this purpose) and include in the exported dataset.
