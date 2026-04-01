######################################################################
# 07_expression_top10_hubs_USP9X.R
# Plot expression of top 10 STRING hub genes (USP9X High vs Low)
######################################################################

library(tidyverse)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("figures", showWarnings = FALSE)
dir.create("figures/string", showWarnings = FALSE)

######################################################################
# 0) CHECK REQUIRED OBJECTS & FILES
######################################################################

# We rely on in-memory objects from 00_preprocess_STAR_UPDATED.R
if (!exists("expr_small") | !exists("meta_small")) {
  stop("Run 00_preprocess_STAR_UPDATED.R in this R session first so expr_small and meta_small exist.")
}

# Top 10 hubs from step 06 / 05
top10_file <- "results/string/STRING_top10_hubs.tsv"
if (!file.exists(top10_file)) {
  stop("File not found: ", top10_file,
       "\nRun 05_STRING_network_USP9X.R (and 06 if needed) first.")
}

deg_file <- "results/DEG_limma_USP9X_High_vs_Low.tsv"
if (!file.exists(deg_file)) {
  stop("DEG file not found: ", deg_file,
       "\nRun 03_DEG_limma.R first.")
}

######################################################################
# 1) LOAD TOP 10 HUBS + DEG TABLE
######################################################################

top10_hubs <- readr::read_tsv(top10_file, show_col_types = FALSE)
top10_symbols <- top10_hubs$gene

cat("✅ Top 10 hub symbols:\n")
print(top10_symbols)

deg_res <- readr::read_tsv(deg_file, show_col_types = FALSE) %>%
  dplyr::mutate(ensembl_clean = sub("\\..*$", "", ensembl_id))

######################################################################
# 2) MAP ENSEMBL -> SYMBOL (same logic as in step 05)
######################################################################

gene_map <- bitr(
  deg_res$ensembl_clean,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)

deg_mapped <- deg_res %>%
  dplyr::inner_join(gene_map, by = c("ensembl_clean" = "ENSEMBL")) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# keep mapping from SYMBOL -> Ensembl rownames in expr_small
hub_map <- deg_mapped %>%
  dplyr::filter(SYMBOL %in% top10_symbols) %>%
  dplyr::mutate(SYMBOL = factor(SYMBOL, levels = top10_symbols)) %>%
  dplyr::arrange(SYMBOL)

cat("\n✅ Hub mapping SYMBOL -> Ensembl:\n")
print(hub_map[, c("SYMBOL", "ensembl_id")])

ens_ids <- hub_map$ensembl_id

######################################################################
# 3) SUBSET expr_small TO HUB GENES
######################################################################

# expr_small rownames are Ensembl IDs (with version), same as ensembl_id
common_ens <- intersect(ens_ids, rownames(expr_small))
cat("\nCommon Ensembl IDs between hubs and expr_small:", length(common_ens), "\n")

if (length(common_ens) < 2) {
  stop("Found fewer than 2 hub genes in expr_small. Check ID mapping.")
}

# align metadata with expression (as you did in 02/03)
meta_small <- meta_small %>% dplyr::distinct(short_barcode, .keep_all = TRUE)
expr_small <- expr_small[, !duplicated(colnames(expr_small))]
expr_small <- expr_small[, meta_small$short_barcode]

hub_expr <- expr_small[common_ens, meta_small$short_barcode, drop = FALSE]

# rename rows to SYMBOLs for prettier plots
symbol_order <- hub_map$SYMBOL[match(rownames(hub_expr), hub_map$ensembl_id)]
rownames(hub_expr) <- as.character(symbol_order)

cat("\n✅ Hub expression matrix dimensions (genes x samples): ",
    nrow(hub_expr), "x", ncol(hub_expr), "\n")

######################################################################
# 4) Z-SCORE ROWS & HEATMAP
######################################################################

hub_expr_z <- t(scale(t(hub_expr)))  # z-score each gene (row)

if (!"USP9X_group" %in% colnames(meta_small)) {
  stop("meta_small must contain USP9X_group (High/Low).")
}

ann <- data.frame(
  USP9X_group = meta_small$USP9X_group
)
rownames(ann) <- meta_small$short_barcode

png("figures/string/Heatmap_top10_hubs_expr.png",
    width = 2000, height = 1600, res = 220)

pheatmap(
  hub_expr_z,
  annotation_col           = ann,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames            = TRUE,
  show_colnames            = FALSE,
  main = "Top 10 hub genes – z-scored expression (USP9X High vs Low)"
)

dev.off()

######################################################################
# 5) BOXLOT FOR 1 HUB GENE (FIRST IN LIST)
######################################################################

first_gene <- rownames(hub_expr)[1]

plot_df <- data.frame(
  expr          = as.numeric(hub_expr[first_gene, ]),
  sample        = meta_small$short_barcode,
  USP9X_group   = meta_small$USP9X_group
)

p_box <- ggplot(plot_df, aes(x = USP9X_group, y = expr, fill = USP9X_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Expression of hub gene ", first_gene),
    x     = "",
    y     = "Log-scale expression (expr_small)"
  )

ggsave(paste0("figures/string/Boxplot_", first_gene, "_hub_expr.png"),
       p_box, width = 6, height = 5, dpi = 300)

cat("\n✅ Step 07 complete. Heatmap + boxplot in figures/string/\n")

