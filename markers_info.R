# Shistosoma mansoni microsatellite information
# created by: Jessica Blanton
# last updated: 1/08/2022
# adapted by: Luciano Kalabric


# These lists imports and organize the expected sizes of the microsatellite sets
# for S. mansoni. These values were translated from "Microsatellite information 20220104 FC-REB2.xlsx", 
# assembled by W. Blank, G. Sabino, F. MacAllister, R.E. Blanton.

# Most recent modification reflects a shift of 1F8A fragments to be 1bp smaller than previously used

# Markers are multiplex for eletrophoresis by sets
# Set or batch 1
set <- "set1"
marker <- (marker= c("smms2", "smms13", "smms16"))
ch <- (ch= c(1, 4, 2))
dye <- (dye = c("FAM", "PET", "VIC"))
x_min <- (x_min= c(219, 171, 201)) # Allele with 1 repeat
x_max <- (x_max= c(295, 228, 258)) # Allele with 20 repeats
set1 <- data.frame(set, marker, ch, dye, x_min, x_max)

# Set or batch 2
set <- "set2"
marker <- (marker= c("smms3", "smms17", "smms18", "smms21"))
ch <- (ch= c(1, 2, 2, 4))
dye <- (dye = c("FAM", "VIC", "VIC", "PET"))
x_min <- (x_min= c(159, 271, 188, 160)) # Allele with 1 repeat
x_max <- (x_max= c(216, 328, 245, 217)) # Allele with 20 repeats
set2 <- data.frame(set, marker, ch, dye, x_min, x_max)

# Set or batch 4
set <- "set4"
marker <- (marker= c("sm13taga", "sm13-410", "sm1f8a", "smda23"))
ch <- (ch= c(1, 4, 2, 3))
dye <- (dye = c("FAM", "PET", "VIC", "NED"))
x_min <- (x_min= c(86, 179, 137, 175)) # Allele with 1 repeat
x_max <- (x_max= c(162, 236, 194, 251)) # Allele with 20 repeats
set4 <- data.frame(set, marker, ch, dye, x_min, x_max)

# Set or batch 5
set <- "set5"
marker <- (marker= c("sm29e6a", "sm13-478", "smu31768", "sm15j15a"))
ch <- (ch= c(1, 1, 3, 2))
dye <- (dye = c("FAM", "FAM", "NED", "VIC"))
x_min <- (x_min= c(150, 209, 179, 193)) # Allele with 1 repeat
x_max <- (x_max= c(207, 266, 236, 250)) # Allele with 20 repeats
set5 <- data.frame(set, marker, ch, dye, x_min, x_max)

# Set or batch 6
set <- "set6"
marker <- (marker= c("lg3_sc36b", "sc23b", "smd28"))
ch <- (ch= c(3, 3, 1))
dye <- (dye = c("NED", "NED", "FAM"))
x_min <- (x_min= c(223, 179, 228)) # Allele with 1 repeat
x_max <- (x_max= c(280, 236, 285)) # Allele with 20 repeats
set6 <- data.frame(set, marker, ch, dye, x_min, x_max)

# Set or batch 7
set <- "set7"
marker <- (marker= c("l46951", "r95529", "lg1_sc276", "lg5_sc475"))
ch <- (ch= c(1, 3, 1, 4))
dye <- (dye = c("FAM", "NED", "FAM", "PET"))
x_min <- (x_min= c(150, 216, 71, 269)) # Allele with 1 repeat
x_max <- (x_max= c(207, 273, 128, 326)) # Allele with 20 repeats
set7 <- data.frame(set, marker, ch, dye, x_min, x_max)

my_marker_sets <- rbind(set1, set2, set4, set5, set6, set7)

### MS painel

# channel 1; 6-FAM (blue) 
smms2 <- c(211, 215, 219, 223, 227, 231, 235, 239)
smms3 <- c(177, 180, 183, 186, 189, 192, 195, 198, 201, 204, 207)
sm13taga <- c(102, 106, 110, 114, 118, 122, 126, 130, 134, 138)
sm29e6a <- c(153, 156, 159, 162, 165, 168, 171, 174)
`sm13-478` <- c(224, 227, 230, 233, 236, 239, 242, 245, 248, 251, 254, 257)
smd28 <- c(225, 228, 231, 234, 237, 240)
l46951 <- c(156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189)
lg1_sc276 <- c(98, 101, 104, 107, 110)
ch1_list <- c("smms2", "smms3", "sm13taga", "sm29e6a", "sm13-478", "smd28", "l46951", "lg1_sc276")

# channel 2; VIC (green) 
smms16 <- c(210, 213, 216, 219, 222, 225, 228, 231, 234)
smms17 <- c(286, 289, 292, 295, 298, 301, 304, 307)
smms18 <- c(195, 198, 201, 204, 207, 210, 213, 216, 219, 222, 225, 228, 231)
sm1f8a <- c(149, 152, 155, 158, 161, 164, 167, 170)
sm15j15a <- c(208, 211, 214, 217, 220, 223, 226, 229, 232)
ch2_list <- c("smms16", "smms17", "smms18", "sm1f8a", "sm15j15a")

# channel 3; NED (yellow) 
smda23 <- c(191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235)
smu31768 <- c(185, 188, 191, 194, 197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227)
lg3_sc36b <- c(232, 235, 238, 241, 244, 247, 250, 253, 256, 259, 262, 265, 268)
sc23b <- c(191, 194, 197, 200, 203, 206, 209, 212, 215, 218)
r95529 <- c(219, 222, 225, 228, 231, 234, 237, 240, 243, 246, 249, 252)
ch3_list <-c("smda23", "smu31768", "lg3_sc36b", "sc23b", "r95529")

# channel 4; PET (red) 
smms13 <- c(183, 186, 189, 192, 195, 198, 201, 204)
smms21 <- c(172, 175, 178, 181, 184, 187)
`sm13-410` <- c(191, 194, 197, 200, 203)
lg5_sc475 <- c(281, 284, 287, 290, 293, 296, 299, 302, 305)
ch4_list <- c("smms13", "smms21", "sm13-410", "lg5_sc475")

smms_dye_list <- list(ch1_list, ch2_list, ch3_list, ch4_list)
names(smms_dye_list) <- c("FAM", "VIC", "NED", "PET")

