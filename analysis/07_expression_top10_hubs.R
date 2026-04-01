######################################################################
# 07_expression_top10_hubs.R
######################################################################

library(tidyverse)
library(pheatmap)
library(ggplot2)

# 1) Load expression matrix (genes x samples) and sample metadata
# >>>> change these paths/names to whatever you actually saved <<<<
expr_mat <- readr::read_tsv("results/vst_expr_matrix_USP9X.tsv")
sample_info <- readr::read_tsv("results/sample_info_USP9X.tsv")
# assume expr_mat has a column "SYMBOL" then one column per sample
expr_mat <- expr_mat %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

# make sure sample columns match sample IDs in metadata
expr_mat <- expr_mat[, sample_info$sample_id]

# 2) Load top 10 hubs from step 06
top10_hubs <- readr::read_tsv("results/string/STRING_top10_hubs.tsv")
hub_genes <- intersect(top10_hubs$gene, rownames(expr_mat))

cat("✅ Hub genes found in expression matrix:", hub_genes, "\n")

# 3) Subset and z-score by gene
hub_expr <- expr_mat[hub_genes, , drop = FALSE]
hub_expr_z <- t(scale(t(hub_expr)))  # z-score rows

# 4) Heatmap (USP9X High vs Low as annotation)
ann_col <- data.frame(
  USP9X_group = sample_info$USP9X_group   # e.g. "High"/"Low"
)
rownames(ann_col) <- sample_info$sample_id

pheatmap(
  hub_expr_z,
  annotation_col = ann_col,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Top 10 hub genes – z-scored expression"
)

ggsave("figures/string/exp_heatmap_top10_hubs.png",
       width = 7, height = 6, dpi = 300)

# 5) Example boxplot for the strongest hub
top_gene1 <- hub_genes[1]

plot_df <- data.frame(
  expr = hub_expr[top_gene1, ],
  sample = colnames(hub_expr),
  USP9X_group = sample_info$USP9X_group
)

ggplot(plot_df, aes(x = USP9X_group, y = expr, fill = USP9X_group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(
    title = paste0("Expression of hub gene ", top_gene1),
    x = "USP9X group", y = "Normalized expression"
  ) +
  theme_minimal()

ggsave(paste0("figures/string/boxplot_", top_gene1, "_by_USP9X_group.png"),
       width = 5, height = 4, dpi = 300)

