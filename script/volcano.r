# Load necessary packages
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(Maaslin2)

# Load the data for genes output from PICRUSt2
data <- read_tsv("<output_path_gene.tsv>")

# Prepare metadata for conditions
sample_conditions <- data.frame(
  sample = colnames(data)[2:ncol(data)], # assuming the first column is KO terms
  condition = c(rep("<sample-group>", <number-of-samples>), rep("<sample-group>", <number-of-samples>))
)

# Transpose the data so that samples are rows and features (KOs) are columns for MaAsLin2
data_transposed <- t(data[,-1])  # Remove KO terms and transpose
colnames(data_transposed) <- data$gene  # Assign gene names as column names

# Prepare the input list for MaAsLin2
input_data <- list(
  "data" = data_transposed,         # The omics data
  "metadata" = sample_conditions    # The metadata
)

# Ensure the metadata has rownames corresponding to the sample names
rownames(input_data$metadata) <- input_data$metadata$sample

# Run MaAsLin2 for differential abundance analysis
maaslin_results <- Maaslin2(
  input_data$data,
  input_data$metadata,
  fixed_effects = c("condition"),
  output = "MaAsLin2_results"
)

# Load the results from MaAsLin2
results_df <- read_tsv("MaAsLin2_results/all_results.tsv")

# Filter significant results and subset the data for both groups
# Generate volcano plot