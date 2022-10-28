# peakscan
This is the first chunck for microsatellite analysis and uses the library Fragman developed by Giovanny Covarrubias-Pazaran et al.

# Packages and software required
Fragman 1.0.9 - Peak scoring functions<br>
tidyr 1.1.4 - Tidyverse best practices<br>
dplyr 1.0.7 - Tidyverse best practices<br>
magrittr 2.0.1 - Use of L->R piping of functions<br>
qpdf 1.1 - PDF manipulations<br>

# Data Required
  .fsa files
  maker_info.R
  ladder_info.R

# Algorithm proposed by Jessica Blanton et al.
Pseudocode:
1) Read fsa binary files (new)
2) Check if the fsa directory has fsa files from the same batch: source("check_fsa_v_batch.R")
3) Retrieve metadata from fsa files into a table and a file called fsa_info.txt: source("get_fsa_metadata.R")
4) Retrieve the eletrophoretic data and show plots from each channel: source("storing_inds_rev3.R")
5) Associate the actual dyes to the channels ("associate_dye_names.R")

# Proposed modifications
1) Save the .fsa files from each marker set in different folder (ie. fsa_set1, fsa_set2...)
2) Structure markers sets and markers info in the source maker_info.R
3) Read .fsa files for particular marker set and analyse all markers at once
5) Use the markers_info to create panels for peak scoring
6) Report .fsa files with bad ladder
7) Check and edit manually the electropherograms for size matching with expected allele sizes
8) Export the results in a data.frame including sample_ID, peak sizes and peak scores (heights) 
