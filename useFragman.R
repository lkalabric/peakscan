```{r, chunk_header}
# Script name: useFragman.R
# Author: Luciano Kalabric
# Purpose: Literally implement package 'Fragman' to peak scan fsa files
# Creation date: Sept 23 2022
# Last update: Sept 24 2022
```
# Check required packages and autoinstall ======================================
list.of.packages <- c("Fragman","pacman")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Simpler version to check and autointall packages in R ========================
# if (!require('Fragman')) install.packages('Fragman'); library('Fragman')

# Set the working directory ====================================================
# setwd("~/GitHub/peakscan/peak-analysis/")
setwd("C:/Users/kalab/OneDrive/Documentos/GitHub/peakscan/peak-analysis")
# Original method

# Loads package Fragman ========================================================
library('Fragman')

# Use of package Fragman =======================================================
# The core of the package Fragman and the workflow of the fragment analysis rely in 4 functions:
# 1) storing.inds (function in charge of reading the FSA or txt(CQS) files and storing them with a list structure)
fsa_dir <- "~/GitHub/peakscan/fsa"
fsa_data <- storing.inds(fsa_dir)
# my_samples <- storing.inds(fsa_dir, channels=5, fourier=TRUE, 
#                              saturated=TRUE, lets.pullup=TRUE, 
#                                plotting=FALSE, rawPlot=FALSE,
#                                  llength=3000, ulength=80000)
# JB method
source("storing_inds.R")
fsa_data <- storing_inds(fsa_dir, channels = 5, rawPlot = TRUE, fourier = TRUE, saturated = TRUE, lets.pullup = FALSE)

# 2) ladder.info.attach (uses the information read from the FSA files and a vector containing the
#                       ladder information (DNA size of the fragments) and matches the peaks from the channel where
#                       the ladder was run with the DNA sizes for all samples. Then loads such information in the R
#                       environment for the use of posterior functions)

# 3) overview2 (create friendly plots for any number of individuals specified and can be used to
#              design panels (overview2) for posterior automatic scoring (like licensed software does), or make
#              Fragman-package 3
#              manual scoring (overview) of individuals such as parents of biparental populations or diversity
#              populations)

# 4) The score.markers (function score the alleles by finding the peaks provided in the panel (if
#              provided), otherwise returns all peaks present in the channel). Thisfinal function can be automatized
#              if several markers are located in the same channel by creating lists of panels taking advantage of R
#              capabilities and data structures.

# Set working directory for this session's output
setwd("GitHub/P3_Pipeline")

# Adding other pieces of code
source("manual_rescore.R")
source("replace_scores_df.R")
source("storing_inds_rev3.R")

# Designate data location
fsa_dir <- "fsa"

# Import files, show plots
fsa_data <- storing_inds_rev3(fsa_dir, channels = 5, rawPlot = FALSE, 
                              fourier = TRUE, saturated = TRUE, lets.pullup = FALSE)
# Code copied from Fragman manual


# Define vector of internal standard/ladder sizes
GS600LIZ <- c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 
              250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 
              460, 480, 500, 514, 520, 540, 560, 580, 600)

# Input and match ladder with the samples, with visual check:
ladder.info.attach(stored = fsa_data, 
                   ladder = GS600LIZ, 
                   ladd.init.thresh = 200,
                   prog = FALSE,
                   draw = FALSE)

# Score markers for all files, or specify a subset to score with `n.inds`
scores_29E6A <- score_markers_rev3(my.inds = fsa_data, 
                             channel = 1, 
                             channel.ladder = 5,
                             panel = "mic_29E6A", 
                             ladder = GS600LIZ, 
                             init.thresh = 100,
                             ploidy = length(mic_29E6A),
                             windowL = 2,
                             windowR= 0.5,
                             shift = 1.5, 
                             left.cond = c(0, 2.5),
                             right.cond = 0,
                             pref = 1, 
                             plotting = FALSE)
# Manually select peaks directly on the plotting pane
# Retrieve index locations of samples to be scored:
correct_inds <- which(names(my.inds) %in%
                        c("104.1a_V_4_Sample_20201013_110136.fsa",
                          "75.1a_V_1_Sample_20201013_110133.fsa",
                          "cwru_20.1_V_FA062920_2020-06-29_G12.fsa"))

# Cursor-click re-scoring a panel for specific samples
scores_rerun <- manual_rescore(my.inds = fsa_data, 
               channel = 1,
               n.inds = correct_inds,
               panel = mic_29E6A,
               ylim= c(0,4000),
               ladder = GS600LIZ)

# Integrate values from re-scored samples into prior results, with output as a data frame
# Script inputs:
# 1.	Original scores list object
# 2.	Second rerun scores list object (from auto or manual commands)
# 3.	Marker panel
# rerun scores are result of either automatic scoring with `score_markers_rev3()` or manual scoring with `scores_rerun`
scores_29E6A_rr_df <- replace_scores_df(scores_29E6A, scores_rerun, "mic_29E6A")
