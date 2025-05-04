#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

# Load edgeR
library(edgeR)

setwd("~/Desktop/count_output")

# Read count data
combined_matrix <- read.table("combined_counts.txt")

# Create DGEList object
dge <- DGEList(counts = combined_matrix, group = colnames(combined_matrix))

# Calculate normalization factors
dge <- calcNormFactors(dge)

# Define groups (make sure these match your column names in combined_matrix)
groups <- factor(c("Blood", "Brain_ant", "Brain_post", 
                   "Branchial_Pouch", "Egg", "Gut_ant", "Gut_post", 
                   "Heart", "Notochord", "Olfactory_bulb", 
                   "Sensory_tentacles", "Slime_Gland"))

# Create valid group names
valid_groups <- factor(make.names(groups))

# Create design matrix
design <- model.matrix(~0 + valid_groups)
colnames(design) <- levels(valid_groups)

# Set dispersion (adjust value based on your data)
dge$common.dispersion <- 0.1

# Fit the model
fit <- glmFit(dge, design)

# Define the specific comparisons we want with output filenames
comparisons <- list(
  "Sensory_tentacles_vs_Brain_post" = list(contrast = c("Sensory_tentacles", "Brain_post"), 
                                           filename = "DE_SensoryTentacles_vs_BrainPost.csv"),
  "Sensory_tentacles_vs_Brain_ant" = list(contrast = c("Sensory_tentacles", "Brain_ant"), 
                                          filename = "DE_SensoryTentacles_vs_BrainAnt.csv"),
  "Sensory_tentacles_vs_Heart" = list(contrast = c("Sensory_tentacles", "Heart"), 
                                      filename = "DE_SensoryTentacles_vs_Heart.csv"),
  "Sensory_tentacles_vs_Egg" = list(contrast = c("Sensory_tentacles", "Egg"), 
                                    filename = "DE_SensoryTentacles_vs_Egg.csv")
)

# Perform each comparison and save results
for (comp_name in names(comparisons)) {
  pair <- comparisons[[comp_name]]$contrast
  filename <- comparisons[[comp_name]]$filename
  
  # Create contrast
  contrast <- makeContrasts(
    contrasts = paste(pair[1], "-", pair[2], sep = ""),
    levels = design
  )
  
  # Perform likelihood ratio test
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Get all genes
  top_genes <- topTags(lrt, n = Inf)
  
  # Save results to separate file
  write.csv(top_genes$table, file = filename, row.names = TRUE)
  
  # Print progress
  print(paste("Saved results for", comp_name, "to", filename))
}

print("All comparisons completed and saved to separate files.")
