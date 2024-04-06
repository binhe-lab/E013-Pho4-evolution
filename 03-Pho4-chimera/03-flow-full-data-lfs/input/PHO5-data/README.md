QC for the flow data using scripts I developed before. The HTML output files are large and not included in the repository. The markdown files should contain some of the output. The main steps take are:

1. a rectangular gate to remove non-cell events.
1. a single cell gate to remove doublets.
1. a clustering gate to identify a homogeneous population.
1. calculate a cell-size normalized fluorescence value (not optimal for this purpose) and include in the exported dataset.

The data files in this folder comes from multiple days and can be grouped into sets. Below is a table explaining how the data are grouped

| Dates                   | Experiment                                      | Comment                                                      |
| ----------------------- | ----------------------------------------------- | ------------------------------------------------------------ |
| 2023.02.08 - 2023.03.31 | Main set of chimeras with 5 region split        | Contains an incomplete 6 region split design, where Lindsey tested some region 4 splits; also contains some alternative breakpoint designs. This is the main dataset for the paper. |
| 2023.12.28 - 2023.12.29 | A complete set with alternative 4/5 breakpoints | In the main set, there were a few chimeras with different break points. Among them, those with an alternative, asymmetric 4/5 breakpoints rescued some of the non-functional chimera in the main set. In this experiment, the goal is to make a complete set with that alternative 4/5 breakpoints and decide which one to use for the paper. |
| 2024.01.22 - 2024.01.29 | Sc45 and Cg4ext test                            | Most of the chimeras in the main set with P2ID:Cg + DBD:Sc are non functional. To test the hypothesis that DBD:Sc requires P2ID:Sc to be active, Lindsey constructed a few of those non-functional chimeras, replacing DBD:Sc with a combined P2ID-DBD:Sc. In addition, three chimeras were made to test a different idea related to the same non functional set. This time, the test is to see if fusing a small stretch of residues that are predicted to be an extention of the basic region and part of the first alpha helix only in CgPho4 to DBD:Sc would rescue those chimeras. |

