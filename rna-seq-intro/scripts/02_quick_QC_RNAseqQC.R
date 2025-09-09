# ===============================================
# Bulk RNA-seq QC overview
# ===============================================

# Load required packages

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(RNAseqQC) #if needed: install.packages("RNAseqQC")
library(PCAtools) # BiocManager::install("PCAtools")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tibble)
library(dplyr)
library(tidyr)
library(matrixStats)
library(cowplot)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load sample metadata
metadata <- read.table("NGS2025_09375_metadata.txt", header = TRUE, 
                       sep = "\t", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$Sample

metadata <- metadata %>%
  mutate(
    Condition = sub(" .*", "", Group),
    Handling = sub(".* ", "", Group)
  )

# Load count data
counts <- read.table("NGS2025_09375.count.genes.txt", header = TRUE, row.names = 1)
counts <- counts[, rownames(metadata)]

# Convert Ensembl IDs to gene symbols
ensembl_ids <- sub("\\..*$", "", rownames(counts))  # remove version numbers

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Create a data frame of counts with gene symbols
count_df <- as.data.frame(counts)
count_df$Symbol <- gene_symbols

# Optionally remove rows without symbols
#count_df <- count_df[!is.na(count_df$Symbol), ]

# Keep ENSG IDs as rownames
counts <- count_df[, setdiff(names(count_df), "Symbol")]

# RNAseqQC metrics
dds <- make_dds(counts = counts, metadata = metadata, ah_record = "AH89426")
dds <- filter_genes(dds, min_count = 5, min_rep = 2)

plot_total_counts(dds)
plot_library_complexity(dds)
plot_gene_detection(dds)
plot_biotypes(dds)

# Variance stabilization & replicate variability
vsd <- vst(dds)
mean_sd_plot(vsd)

ma_plots <- plot_sample_MAs(vsd, group = "Group")
cowplot::plot_grid(plotlist = ma_plots[1:8], ncol = 3)

norm_expr <- assay(vsd)

## Additional QC metrics

# Fraction of zero counts per sample
zero_frac <- colSums(counts == 0) / nrow(counts)
ggplot(data.frame(Sample=names(zero_frac), ZeroFrac=zero_frac), 
       aes(x=Sample, y=ZeroFrac, fill=ZeroFrac)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Fraction of zero counts per sample")

# Sample-sample correlation heatmap
sample_cor <- cor(norm_expr)
pheatmap(sample_cor, annotation_col=metadata[, c("Condition","Handling")],
         main="Sample-sample correlation heatmap")

# Replicate consistency (pairwise correlations)
rep_cor_vals <- sample_cor[lower.tri(sample_cor)]
summary(rep_cor_vals)

# Quick clustering
set.seed(1)
plot_sample_clustering(vsd, anno_vars = c("Condition", "Handling", "Rep"), distance = "euclidean")

# PCA analysis
p <- pca(norm_expr, metadata = metadata, removeVar = 0.1)

# Scree plot
screeplot(p, axisLab = TRUE, title = "Variance Explained by PCs")

# PCAtools biplot
biplot(p, 
       colby = "Handling", 
       shape = "Condition", 
       lab = metadata$Sample, 
       legendPosition = "right")

# Faceted PCA by Handling using ggplot2
scores <- as.data.frame(p$rotated)
scores$Sample <- rownames(scores)
scores$Condition <- metadata$Condition
scores$Handling <- metadata$Handling

ggplot(scores, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -0.5, size = 3) +
  facet_wrap(~Handling) +
  theme_bw() +
  labs(title = "PCA faceted by Handling", x = "PC1", y = "PC2")

# Sample-level clustering on top var genes
sample_dist <- dist(t(norm_expr))
sample_hclust <- hclust(sample_dist, method = "average")

plot(sample_hclust, labels = metadata$Sample,
     main = "Sample Clustering Dendrogram", cex = 0.8)

clusters <- cutree(sample_hclust, k = 4)
metadata$Cluster <- factor(clusters)

# Heatmap of top variable genes
top_var_genes <- head(order(matrixStats::rowVars(norm_expr), decreasing = TRUE), 500)
pheatmap(norm_expr[top_var_genes, ],
         annotation_col = metadata[, c("Condition","Handling")],
         show_rownames = FALSE,
         cluster_cols = TRUE,
         scale = "row")

# Detect potential outliers
z_scores <- scale(rowMeans(sample_cor))
outliers <- rownames(sample_cor)[abs(z_scores) > 2]
cat("Potential outliers:", outliers, "\n")
