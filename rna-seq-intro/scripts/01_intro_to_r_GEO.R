# 01_intro_to_r_GEO.R
# Load required libraries
library(tidyverse)

# Load data
counts <- read_csv("data/GSE155432_counts.csv")
metadata <- read_csv("data/GSE155432_metadata.csv")

# Inspect structure
print(dim(counts))
print(head(counts))
print(metadata)

# Summarize total counts per sample
total_counts <- colSums(counts[,-1])
barplot(total_counts, main = "Library sizes", col = "steelblue", las = 2)

# Log2 transform
log_counts <- log2(counts[,-1] + 1)
boxplot(as.data.frame(log_counts), las = 2, main = "Log2 Expression Distributions")
