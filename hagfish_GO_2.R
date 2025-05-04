# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("UniProt.ws")
# BiocManager::install("biomaRt")
library(UniProt.ws)
library(biomaRt)
library(ggplot2)
library(dplyr)

setwd("~/Desktop/Gene_Onto")


blast_cols <- c(
  "query_id",       # Query sequence identifier
  "subject_id",     # Subject (database) sequence identifier
  "perc_identity",  # Percentage identity
  "align_length",   # Alignment length
  "mismatches",     # Number of mismatches
  "gap_opens",      # Number of gap openings
  "q_start",        # Start position in query
  "q_end",          # End position in query
  "s_start",        # Start position in subject
  "s_end",          # End position in subject
  "evalue",         # E-value (significance)
  "bit_score"       # Bit score
)

blast_ionchannel <- read.table("pat_ah2p_vs_uniprot_ion_channel.txt",     
                               header = FALSE,           
                               sep = "\t",               
                               col.names = blast_cols,   
                               stringsAsFactors = FALSE
)


sens_egg <- read.csv("DE_ic_sens_egg.csv", check.names = FALSE)
sens_brainant <- read.csv("DE_ic_sens_brainant.csv", check.names = FALSE)
sens_brainpos <- read.csv("DE_ic_sens_brainpos.csv", check.names = FALSE)
sens_heart <- read.csv("DE_ic_sens_heart.csv", check.names = FALSE)

extract_spID <- function(df) {
  
  # in the 8th column, make 
  
}
