#############################################################
#
# Ref to the ARTICLE
# 
# Manuel Blank & Davide Bulgarelli
# mblank@dundee.ac.uk
# d.bulgarelli@dundee.ac.uk
# 
# Revison May 2017
# 
# script to reproduce ANCOM calculations presented in the manuscript
# 
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#Requirements for ANCOM package (installed from zip file)
#source("http://bioconductor.org/biocLite.R")
#biocLite("doParallel")
#biocLite("foreach")
#biocLite("ggplot2")
#biocLite("Rcpp")
#biocLite("shiny")
#biocLite("exactRankTests")
#biocLite("openxlsx")
#biocLite("DT")
#biocLite("coin")
#biocLite("grid")
#biocLite("futile.logger")


#Requirements for ANCOM package (installed from zip file)
library("doParallel")
library("foreach")
library("ggplot2")
library("Rcpp")
library("shiny")
library("exactRankTests")
library("openxlsx")
library("DT")
library("coin")
library("ancom.R")


#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set the working directory (mind that your working directory has a different path)
setwd("G:/Uni/Degree/Y5_MSci/Senga_rhm_manuscript/RHM_ancom")


#############################################################
#                     ANCOM analysis                        #
#############################################################


#############################################################
# Import OTU table in 'wide format', generated manually in Excel.
OTU_table_transposed <- read.delim("OTU_table_transposed.txt")

# ANCOM calculation
ancom_rhm_phyla <- ANCOM(OTU_table_transposed, multcorr = 2)

# Generate data frame with output of the calculation
ancom_result_dataframe <- as.data.frame(ancom_rhm_phyla$detected)

# Save result of ANCOM analysis as data frame
#write.table(ancom_result_dataframe, file="ws7.txt", sep="\t")

# End