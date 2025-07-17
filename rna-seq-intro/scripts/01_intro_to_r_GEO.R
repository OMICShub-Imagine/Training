# 01_intro_to_r_GEO.R
# Load required libraries
library(tidyverse) # Includes readr, ggplot2, dplyr

# Load data
counts <- read_csv("data/GSE155432_counts.csv")
metadata <- read_csv("data/GSE155432_metadata.csv")

# Explore data structure
dim(counts)        # Rows = genes, Columns = samples
head(counts)       # View first few genes
colnames(counts)   # Sample names
str(metadata)      # View metadata structure

# Summarize total counts per sample
library_sizes <- colSums(counts[,-1]) # Remove gene column if needed

barplot(
  library_sizes,
  col = "steelblue",
  las = 2,
  main = "Library Size per Sample",
  ylab = "Total Counts"
)

# Log2 transform and visualize expression
log_counts <- log2(counts[,-1] + 1)
boxplot(
  as.data.frame(log_counts),
  las = 2,
  main = "Log2 Transformed Expression",
  ylab = "log2(counts + 1)"
)
