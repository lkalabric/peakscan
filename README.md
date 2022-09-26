# peakscan
Peak scan analysis for microsatellite allelic discrimination and scoring

# Algorithm adopted from Walter Blank's Excel plugins
Pseudocode:
1) Read fsa binary files (new)
2) Check if the fsa directory has fsa files from the same batch: source("check_fsa_v_batch.R")
3) Retrieve metadata from fsa files into a table and a file called fsa_info.txt: source("get_fsa_metadata.R")
4) Retrieve the eletrophoretic data and show plots from each channel: source("storing_inds_rev3.R")
5) Associate the actual dyes to the channels ("associate_dye_names.R")
