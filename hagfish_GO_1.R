# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("UniProt.ws")
# BiocManager::install("biomaRt")
library(UniProt.ws)
library(biomaRt)
library(ggplot2)
library(dplyr)

setwd("~/Desktop/count_output")

DE_sens_egg <- read.csv("DE_SensoryTentacles_vs_Egg.csv", check.names = FALSE)
DE_sens_brainant <- read.csv("DE_SensoryTentacles_vs_BrainAnt.csv", check.names = FALSE)
DE_sens_brainpos <- read.csv("DE_SensoryTentacles_vs_BrainPost.csv", check.names = FALSE)
DE_sens_heart <- read.csv("DE_SensoryTentacles_vs_Heart.csv", check.names = FALSE)

# extract all significantly expressed genes using p-value
DE_sens_egg <- DE_sens_egg[DE_sens_egg$PValue < 0.01, ]
DE_sens_brainant <- DE_sens_brainant[DE_sens_brainant$PValue < 0.01, ]
DE_sens_brainpos <- DE_sens_brainpos[DE_sens_brainpos$PValue < 0.01, ]
DE_sens_heart <- DE_sens_heart[DE_sens_heart$PValue < 0.01, ]



# read in protein blast results
setwd("~/Desktop")

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

blast_result <- read.table("blast_results.txt",     
                           header = FALSE,           # No header in the file
                           sep = "\t",               # Tab-delimited (default for BLAST -outfmt 6)
                           col.names = blast_cols,   # Assign column names
                           stringsAsFactors = FALSE  # Optional: avoid converting strings to factors
                          )


setwd("~/Desktop/Gene_Onto")

blast_ionchannel <- read.table("pat_ah2p_vs_uniprot_ion_channel.txt",     
                               header = FALSE,           
                               sep = "\t",               
                               col.names = blast_cols,   
                               stringsAsFactors = FALSE
                              )

# Extract gene names from each DE file
geneid_sens_egg <- DE_sens_egg[,1]
geneid_sens_brainant <- DE_sens_brainant[,1]
geneid_sens_brainpos <- DE_sens_brainpos[,1]
geneid_sens_heart <- DE_sens_heart[,1]

# Filter BLAST results for each comparison
#mapped_sens_egg <- blast_result[blast_result$query_id %in% geneid_sens_egg, ]
#mapped_sens_brainant <- blast_result[blast_result$query_id %in% geneid_sens_brainant, ]
#mapped_sens_brainpos <- blast_result[blast_result$query_id %in% geneid_sens_brainpos, ]
#mapped_sens_heart <- blast_result[blast_result$query_id %in% geneid_sens_heart, ]


# Filter ion channel BLAST results for each comparison
ic_mapped_sens_egg <- blast_ionchannel[blast_ionchannel$query_id %in% geneid_sens_egg, ]
ic_mapped_sens_brainant <- blast_ionchannel[blast_ionchannel$query_id %in% geneid_sens_brainant, ]
ic_mapped_sens_brainpos <- blast_ionchannel[blast_ionchannel$query_id %in% geneid_sens_brainpos, ]
ic_mapped_sens_heart <- blast_ionchannel[blast_ionchannel$query_id %in% geneid_sens_heart, ]


# filter out the good matches
good_blast <- function(df) {
  df %>% 
    filter(perc_identity > 50, evalue < 1e-10) %>%
    distinct(query_id) %>%  # Get unique query IDs
    pull(query_id)          # Extract as a vector
}

# Apply filtering to each tissue comparison (values stored in vector)
ic_siggenes_sens_egg <- good_blast(ic_mapped_sens_egg)
ic_siggenes_sens_brainant <- good_blast(ic_mapped_sens_brainant)
ic_siggenes_sens_brainpos <- good_blast(ic_mapped_sens_brainpos)
ic_siggenes_sens_heart <- good_blast(ic_mapped_sens_heart)

# Filter DE results to keep only ion channel genes that passed BLAST criteria
DE_ic_sens_egg <- DE_sens_egg[DE_sens_egg[,1] %in% ic_siggenes_sens_egg, ]
DE_ic_sens_brainant <- DE_sens_brainant[DE_sens_brainant[,1] %in% ic_siggenes_sens_brainant, ]
DE_ic_sens_brainpos <- DE_sens_brainpos[DE_sens_brainpos[,1] %in% ic_siggenes_sens_brainpos, ]
DE_ic_sens_heart <- DE_sens_heart[DE_sens_heart[,1] %in% ic_siggenes_sens_heart, ]

# name the first column
name_gene_id <- function(df) {
  names(df)[1] <- "gene_id"  # Rename the first column to "gene_id"
  return(df)
}

# Apply to all DE_ic_sens_{tissue} data frames
DE_ic_sens_egg <- name_gene_id(DE_ic_sens_egg)
DE_ic_sens_brainant <- name_gene_id(DE_ic_sens_brainant)
DE_ic_sens_brainpos <- name_gene_id(DE_ic_sens_brainpos)
DE_ic_sens_heart <- name_gene_id(DE_ic_sens_heart)


# for each match found in ic_mapped_sens_{tissue} with each entry stored 
# in the vector ic_siggenes_sens_{tissue}, make a corresponding list of 
# the subject id. Then merge the subject id to DE_ic_sens_{tissue}.

# Function to map query IDs to all their subject IDs
get_subject_ids <- function(blast_df, query_ids) {
  blast_df %>%
    filter(query_id %in% query_ids) %>%  # Keep only rows matching the query IDs
    group_by(query_id) %>%              # Group by query ID
    summarise(subject_ids = list(subject_id))  # Store all subject IDs as a list per query
}

# Apply the function to each tissue comparison
subject_map_sens_egg <- get_subject_ids(ic_mapped_sens_egg, ic_siggenes_sens_egg)
subject_map_sens_brainant <- get_subject_ids(ic_mapped_sens_brainant, ic_siggenes_sens_brainant)
subject_map_sens_brainpos <- get_subject_ids(ic_mapped_sens_brainpos, ic_siggenes_sens_brainpos)
subject_map_sens_heart <- get_subject_ids(ic_mapped_sens_heart, ic_siggenes_sens_heart)

# Function to merge subject IDs
merge_subject_ids <- function(de_df, subject_map_df) {
  subject_map_df %>%
    group_by(query_id) %>%
    summarise(subject_ids = paste(subject_ids, collapse = "; ")) %>%  # Combine multiple IDs
    left_join(de_df, ., by = c("gene_id" = "query_id"))              # Merge with DE data
}

# Apply to each tissue
DE_ic_sens_egg_merged <- merge_subject_ids(DE_ic_sens_egg, subject_map_sens_egg)
DE_ic_sens_brainant_merged <- merge_subject_ids(DE_ic_sens_brainant, subject_map_sens_brainant)
DE_ic_sens_brainpos_merged <- merge_subject_ids(DE_ic_sens_brainpos, subject_map_sens_brainpos)
DE_ic_sens_heart_merged <- merge_subject_ids(DE_ic_sens_heart, subject_map_sens_heart)

write.csv(DE_ic_sens_egg_merged, file = "DE_ic_sens_egg.csv")
write.csv(DE_ic_sens_brainant_merged, file = "DE_ic_sens_brainant.csv")
write.csv(DE_ic_sens_brainpos_merged, file = "DE_ic_sens_brainpos.csv")
write.csv(DE_ic_sens_heart_merged, file = "DE_ic_sens_heart.csv")


