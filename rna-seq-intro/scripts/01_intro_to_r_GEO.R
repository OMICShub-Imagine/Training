# 01_intro_to_r_GEO.R

# Show R version
print(R.version.string)

# Install (once)
install.packages("tidyverse") # Includes readr, ggplot2, dplyr

# Load (every session)
library(tidyverse)

#Setting up working directory (active folder)
setwd("Desktop/Training/rna-seq-intro/") #change with your own path

# Load data
counts <- read_csv("data/GSE60424_counts.csv")
metadata <- read_csv("data/GSE60424_metadata.csv")

# Explore data structure
dim(counts)        # Rows = genes, Columns = samples
head(counts)       # View first few genes
colnames(counts)   # Sample names
str(metadata)      # View metadata structure

# For simplicity, subset metadata to first 10 samples
metadata_sub <- metadata[1:10, ]

# Extract sample names matching counts columns
sample_titles <- metadata_sub$title
samples_in_counts <- intersect(sample_titles, colnames(counts))

cat("Samples selected:", samples_in_counts, "\n")

# Subset counts to genes + these samples
counts_sub <- counts[, c("genenames", samples_in_counts)]

dim(counts_sub)

# Summarize total counts per sample
library_sizes <- colSums(counts_sub[, -1])

barplot(
  library_sizes,
  col = "steelblue",
  las = 2,
  main = "Library Size per Sample",
  ylab = "Total Counts"
)

# Log2 transform and visualize expression
log_counts <- log2(as.matrix(counts_sub[, -1]) + 1)

boxplot(
  log_counts,
  las = 2,
  main = "Log2 Transformed Expression",
  ylab = "log2(counts + 1)"
)

sessionInfo()
