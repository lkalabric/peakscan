# These commands allow use of the mouse cursor assign peaks to individual samples in the RStudio plotting panel, and to then integrate the new peak calls into the original batch scores.
# Setup: This tutorial assumes you already have your RStudio environment set up for the fragman workflow.
# Load scripts for manual peak scoring and values merging

# Set working directory for this session's output
setwd("GiHub/P3_Pipeline/")

source("manual_rescore.R")
source("replace_scores_df.R")
Example workflow for initial scores


# Designate data location
fsa_dir <- "fsa"

## Import files, show plots
fsa_data <- storing_inds_rev3(fsa_dir, channels = 5, rawPlot = FALSE, 
                              fourier = TRUE, saturated = TRUE, lets.pullup = FALSE)


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
Manually select peaks directly on the plotting pane
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
Integrate values from re-scored samples into prior results, with output as a data frame
Script inputs:
1.	Original scores list object
2.	Second rerun scores list object (from auto or manual commands)
3.	Marker panel
# rerun scores are result of either automatic scoring with `score_markers_rev3()` or manual scoring with `scores_rerun`
scores_29E6A_rr_df <- replace_scores_df(scores_29E6A, scores_rerun, "mic_29E6A")
