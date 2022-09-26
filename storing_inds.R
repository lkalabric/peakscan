# Jessica Blanton
# Last updated 1/7/22

## This script is an adaptation of "storing.inds" from the Fragman package.
## It primarily adjusts for changes in ABI's .fsa file formats.
## This script should run successfully for any format version <=3.
## Retains Fragman functions for fourier transformation, saturated peaks, and pullup correction.

## REVISION 3: 
# - Reinstated and updated channel parameter code - esp that dealt with ambiguous index locations.
# - Uses Dyechannel count specified in file directory, user specified parameter is now used as check only

## MAJOR CHANGES FROM ORIGINAL 
# - 1: Extract only from tags named as "DATA" in the raw file's directory instead of the full import. 
# - 2: Additional correct column selection for v3 formats.
# - 3: Modifications to the implementation of the "channel" parameter.

## MINOR CHANGES AND NOTES
# - 1: Allows for data directory to be specified with a relative path.
# - 2: Script ONLY imports .fsa files, not .txt (streamlines script).
# - 3: Renames channels with associated dye names

## 2) Plot titles and colors assumed a fixed channel order 1:5 of "FAM", "HEX", "NED", "ROX", "LIZ".
## In rare cases this may not be the actual order of dyes recorded, but this issue will certainly be captured downstream.
## Plot titles have simply been modified to reflect the near-certainty of the dye channel names (vs. absolute).

## 3) The "pullup" parameter introduces artifacts to many runs, avoid using.

storing_inds  <- function (folder, channels = NULL, fourier = TRUE, saturated = TRUE, 
                                lets.pullup = FALSE, plotting = FALSE, 
                                rawPlot = FALSE, llength = 3000, ulength = 80000) 
{
  message("This is a revision of the Fragman script storing.inds with the following changes:
  - Now imports and parse up to version 3 of the fsa format.  
  - Rarely used functionality for non-fsa files has been removed.
  For the original instructions and parameters, run '?storing.inds' \n")
  
  require(Fragman)
  
  # rm(list = ls())
  # folder <- "~/Dropbox/REB_Tulane/REB_Genotyping/R_scripts_JB/data_import_test/"
  # folder <- "~/Dropbox/REB_Tulane/REB_Genotyping/R_scripts_JB/test_intro_data"
  # fourier = TRUE
  # saturated = TRUE
  # lets.pullup = FALSE
  # plotting = FALSE
  # rawPlot = FALSE
  # llength = 3000
  # ulength = 80000
  # channels = 5
  # channels = 4
  # channels = NULL
  # folder= fsa_dir
  
  # coli  <- c("blue", "dark green", "yellow3",
  #                  "red", "orange", "purple")
  
  
  listp2 <- dir(folder, "*.fsa$")
  
  if (length(listp2) == 0) {
    stop(paste("We have not found files with extension .fsa. Please \nmake sure this is the right folder:", 
               folder, "\n"), call. = FALSE)
  }
  
  all.inds.mats <- list(NA)
  
  message("Reading FSA files")
  
  # count <- 0
  tot <- length(listp2)
  # pb <- txtProgressBar(style = 3)
  # setTxtProgressBar(pb, 0)
  
  for (i in 1:length(listp2)) {
    
    ## MINOR CHANGE: do not print .fsa "unimplemented legacy type" warnings
    fsaFile <- suppressWarnings(
      read.abif(paste0(folder,"/",listp2[i]))
    )
    
    ## MAJOR CHANGE:  limit search only to "DATA" tags instead of full import.
    ## Excludes occasional additional tags which fit selection length range.
    ## Also retreive dye number 
    
    # lens <- lapply(fsaFile$Data, length)
    # aaa <- table(unlist(lens))
    
    fsaFile.data <- fsaFile$Data[(grepl("DATA|Dye#|DyeN",names(fsaFile$Data)))]
    fsaFile.dir <- fsaFile$Directory[(grepl("DATA|Dye#|DyeN",fsaFile$Directory$name)), ]
    
    lens <- lapply(fsaFile.data, length)
    aaa <- table(unlist(lens))
    
    # This code selects tags containing data of length between the 
    # minimum and maximum number of indexes set for the runs, 
    # as defined by the parameters "llength" and "ulength." 
    # This is likely due to variable fragment locations in formats 
    # from different companies, beyond fsa versions from ABI
    
    # Get channel number from fsa directory
    channels_d <- fsaFile.data$`Dye#.1`
    
    if (is.numeric(channels) && (channels != channels_d)) {
      message(paste0("Your specified channel number of (", channels, ") is not the same as\n the number of dye channels recorded in the fsa file (", channels_d, ").\n Results will use ", channels_d, " channels."))
    }
    
    # Set channel parameter for each loop and accommodate double entries in v3 format
    if (fsaFile$Header$version==300 && !is.null(channels_d)) {
      channels_l <- 2*channels_d
    } else {
      channels_l <- channels_d
    }
    
    # Search for channels candidates, even if channel parameter is set
    cfound  <- as.vector(aaa[which(
      as.numeric(names(aaa)) > llength & 
        as.numeric(names(aaa)) < ulength)])
    
    #  if (is.null(channels)) {
    #   channels_l <- cfound
    # }
    
    
    # For files where location of indexes is not clear, as seen in some files with shorter runtimes
    if ((length(cfound) > 1) & !(channels_l %in% cfound)) {
      
      # if (length(channels_l) > 1) {
      
      cat(paste("\nYour data for file", listp2[1], "has multiple possible places where\nrun indexes could be stored and we don't know which is the correct one.\n"))
      prov <- aaa[which(as.numeric(names(aaa)) > llength & 
                          as.numeric(names(aaa)) < ulength)]
      
      if (fsaFile$Header$version==300){
        prov2 <- matrix(prov, nrow = 1)/2
      } else {
        prov2 <- matrix(prov, nrow = 1)
      }
      
      rownames(prov2) <- "number.of.channels.found"
      colnames(prov2) <- paste("Run_Length", names(prov), "indexes", sep = "_")
      prov2 <- rbind(prov2, prov2)
      rownames(prov2)[2] <- "number.to.type.if.selected"
      prov2[2, 1:ncol(prov2)] <- 1:ncol(prov2)
      
      cat("Please tell us which option has AT LEAST the number of expected channels\n\n")
      
      print(prov2)
      inut <- as.numeric(readline(prompt = "Enter one of the number.to.type: "))
      channels_l <- cfound[inut]
    }
    
    # Get length of runs (select length with entry numbers equal to number of channels)
    real.len <- as.numeric(
      names(aaa)[which(
        aaa == channels_l &
          as.numeric(names(aaa)) > llength &
          as.numeric(names(aaa)) < ulength)])
    
    # Identify which DATA tags are of correct length to be runs
    v0 <- as.vector(which(unlist(lens) == real.len))
    
    # Filter to only those which have identical datasize, and are multiple of dyechannel number
    chsize <- as.numeric(
      names(which((table(fsaFile.dir$datasize[v0]) %% channels_d) == 0)))
    
    v <- intersect(v0, which(fsaFile.dir$datasize==chsize))
    
    # Subset list to tags of correct length
    reads <- fsaFile.data[v]
    
    # Create dataframe of all selected tags
    prov0 <- as.data.frame(do.call(cbind, reads)) 
    
    # Perform additional column selection for format v3 final readouts
    # This selection is based on an assumed pattern to the order 
    # of the DATA tags, there may be a more deductive way to do this,
    # possibly by determining quantitative differences between the 
    # preliminary and final readouts. 
    if (fsaFile$Header$version==300) {
      
      # Get vector of indexes for columns to keep 
      keepcol <- c(
        # Sample channels_l
        (((channels_l-2)/2)+1):(channels_l-2),
        # Ladder channel
        channels_l)
      
      # Select columns to keep
      prov <- prov0[,keepcol]
      
    } else {
      # If not v3
      prov <- prov0
    }
    
    # Add column and row names
    
    dyenames <- fsaFile.data[grepl("DyeN",names(fsaFile.data))]%>%unlist
    
    if (length(dyenames) == ncol(prov)){
      colnames(prov) <-  paste0("ch_", 1:ncol(prov), "__", dyenames)
      rownames(prov) <- paste("index_", 1:nrow(prov), sep = "")
    } else {
      # Add basic column and row names
      colnames(prov) <- paste("channel_", 1:ncol(prov), sep = "")
      rownames(prov) <- paste("index_", 1:nrow(prov), sep = "")
    }
    
    # Add data frame to list of batch imports
    all.inds.mats[[i]] <- prov
    
    prov <- NULL # Reset variable for next loop of import
    
    # Name list entry with original file name
    names(all.inds.mats)[i] <- as.character(listp2[i])
    
  }
  
  
  if (fourier == TRUE) {
    message("Applying Fourier tranformation for smoothing...")
    all.inds.mats <- lapply(all.inds.mats, function(x) {
      apply(x, 2, transfft)
    })
  }
  
  if (saturated == TRUE) {
    message("Checking and correcting for saturated peaks...")
    all.inds.mats <- lapply(all.inds.mats, function(x) {
      apply(x, 2, saturate)
    })
  }
  
  if (lets.pullup == TRUE) {
    message("Applying pull up correction to the samples to decrease noise from channel to channel")
    if (plotting == TRUE) {
      all.inds.mats <- lapply(all.inds.mats, pullup, 
                              channel = channels_l, plotting = TRUE)
    }
    all.inds.mats <- lapply(all.inds.mats, pullup, 
                            channel = channels_l)
  }
  
  if (rawPlot == TRUE) {
    layout(matrix(1:2, 2, 1))
    coli <- c("blue", "dark green", "yellow3",
              "red", "orange", "purple")
    naname <- c("FAM?", "HEX?", "NED?", "ROX?", "LIZ?")
    # naname <- c("FAM", "HEX", "NED", "ROX", "LIZ")
    message("Plotting raw data")
    
    tot <- length(listp2)
    for (i in 1:dim(all.inds.mats[[1]])[2]) {
      plot(all.inds.mats[[1]][, i], col = transp(coli[i], 0.6),
           # plot(all.inds.mats[[1]][, i], col = transp("black", 0.6), 
           type = "l", ylab = "RFU", main = paste0(
             "Channel ",i," of ", dim(all.inds.mats[[1]])[2], " (",naname[i],")"),
           cex.main = 1, las = 2)
      rect(par("usr")[1], par("usr")[3], par("usr")[2],
           par("usr")[4], col = "white")
      if (length(all.inds.mats) > 1) {
        for (j in 1:length(all.inds.mats)) {
          lines(all.inds.mats[[j]][, i], 
                col = transp(coli[i],
                             # lines(all.inds.mats[[j]][, i], col = transp("black",
                             0.2), lwd = 0.6)
        }
      }
    }
  };
  
  layout(matrix(1, 1, 1))
  class(all.inds.mats) <- c("fsa_stored")
  message("\nOutput is a LIST where each element of the list is a DATAFRAME \n
          with the channels in columns for each FSA file.\n")
  return(all.inds.mats)
}
