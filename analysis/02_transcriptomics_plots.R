######################################################################
# 02_transcriptomics_plots_FINAL.R
######################################################################

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)

dir.create("figures", showWarnings = FALSE)
dir.create("figures/transcriptomics", showWarnings = FALSE)

# ---- make sure Step 0 objects exist ----
if (!exists("expr_small") | !exists("meta_small") | !exists("usp9x_df")) {
  stop("Run Step 0 first — expr_small / meta_small / usp9x_df missing.")
}

######################################################################
### 1) USP9X-ONLY SIMPLE PLOTS (your main focus)
######################################################################

# Boxplot + jitter
p_box <- ggplot(usp9x_df, aes(USP9X_group, USP9X_expr, fill=USP9X_group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  geom_jitter(width=0.15, alpha=0.6, size=2) +
  theme_minimal(base_size=14) +
  labs(title="USP9X Expression: High vs Low Tumors",
       x="", y="USP9X expression (log-scale input)")

ggsave("figures/transcriptomics/USP9X_boxplot.png",
       p_box, width=6, height=5, dpi=300)
print(p_box)

# Histogram (distribution)
p_hist <- ggplot(usp9x_df, aes(USP9X_expr)) +
  geom_histogram(bins=30, fill="steelblue", color="white") +
  theme_minimal(base_size=14) +
  labs(title="Distribution of USP9X Expression in Tumors",
       x="USP9X expression", y="Number of samples")

ggsave("figures/transcriptomics/USP9X_histogram.png",
       p_hist, width=6, height=5, dpi=300)
print(p_hist)

######################################################################
### 2) PCA (QC + shows real separation)
######################################################################

expr_small_filtered <- expr_small[rowVars(expr_small) > 0, ]
expr_t <- t(expr_small_filtered)

pca_res <- prcomp(expr_t, scale. = TRUE)

pc_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = meta_small$USP9X_group,
  Sample = rownames(pca_res$x)
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color=Group)) +
  geom_point(size=3, alpha=0.85) +
  theme_minimal(base_size=14) +
  labs(
    title="PCA: USP9X High vs Low Breast Tumors",
    x=paste0("PC1 (", round(pc_var[1],1), "%)"),
    y=paste0("PC2 (", round(pc_var[2],1), "%)")
  ) +
  scale_color_manual(values=c(High="#D73027", Low="#4575B4"))

ggsave("figures/transcriptomics/PCA_plot.png",
       p_pca, width=7, height=5, dpi=300)
print(p_pca)

######################################################################
### 3) HEATMAP: top 200 most variable genes  (FIX DUPLICATES)
######################################################################

# 1) remove duplicate samples in metadata
meta_small <- meta_small %>% distinct(short_barcode, .keep_all = TRUE)

# 2) remove duplicate columns in expression (just in case)
expr_small_filtered <- expr_small_filtered[, !duplicated(colnames(expr_small_filtered))]

# 3) re-order expression columns to match metadata order
expr_small_filtered <- expr_small_filtered[, meta_small$short_barcode]

# top variable genes
gene_var <- rowVars(expr_small_filtered)
top_genes <- names(sort(gene_var, decreasing=TRUE))[1:200]
heat_data <- expr_small_filtered[top_genes, ]

# annotation that matches heatmap columns exactly
ann <- data.frame(USP9X_group = meta_small$USP9X_group)
rownames(ann) <- meta_small$short_barcode

png("figures/transcriptomics/Heatmap_top200.png",
    width=2200, height=1900, res=220)

pheatmap(
  heat_data,
  show_rownames=FALSE,
  show_colnames=FALSE,
  scale="row",
  clustering_distance_rows="euclidean",
  clustering_distance_cols="euclidean",
  clustering_method="complete",
  annotation_col=ann,
  annotation_colors=list(
    USP9X_group=c(High="#D73027", Low="#4575B4")
  ),
  color=colorRampPalette(c("navy","white","firebrick3"))(120),
  main="Top 200 Most Variable Genes"
)

dev.off()


######################################################################
### 4) QC BONUS: Sample–sample correlation heatmap (simple & strong)
######################################################################

sample_cor <- cor(expr_small_filtered, method="pearson")

png("figures/transcriptomics/SampleCorrelation.png",
    width=2000, height=1800, res=220)

pheatmap( 
  sample_cor,
  annotation_col=ann,
  annotation_row=ann,
  show_rownames=FALSE,
  show_colnames=FALSE,
  main="Sample–Sample Correlation (QC)"
)

dev.off()

cat("\n✅ Step 2 done. All plots saved in figures/transcriptomics/\n")

