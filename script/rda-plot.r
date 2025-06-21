# Install and load the required package
library(phyloseq)
library(vegan)
library(tidyverse)
library(grid)

# Load ASV table
asv_table <- read.delim("<path-to-ASV-table>", row.names = 1)

# Load taxonomy
taxonomy <- read.delim("<path-to-taxonomy>", row.names = 1)

# Load metadata
metadata <- read.delim("<path-to-metadata>", row.names = 1)

# Create a phyloseq object to organize your ASV table, taxonomy, and metadata
ps <- phyloseq(
  otu_table(as.matrix(asv_table), taxa_are_rows = TRUE),
  sample_data(metadata),
  tax_table(as.matrix(taxonomy))
)

# Convert the phyloseq object into a format suitable for vegan
asv_matrix <- as(otu_table(ps), "matrix")
env_data <- as(sample_data(ps), "data.frame")

# Find the common sample names between the two datasets
common_samples <- intersect(colnames(asv_matrix), rownames(env_data))

# Subset ASV matrix and environmental data to only include the common samples
asv_matrix <- asv_matrix[, common_samples]
env_data <- env_data[common_samples, ]

# After subsetting, check again that the dimensions are correct and match
ncol(asv_matrix) == nrow(env_data) # Should return TRUE

# Transpose ASV matrix to match the sample-wise structure of env_data
asv_matrix_t <- t(asv_matrix)

# Perform RDA using the rda() function from the vegan package
rda_result <- rda(asv_matrix_t ~ <env-variable1> + <env-variable2> + <env-variable3> + ..., data = env_data)

# Assess the significance of environmental variables using ANOVA-like permutation tests
anova(rda_result, by = "term", permutations = 999)

# Plot the RDA