######################################################################
# 08_USP9X_hub_module_correlation.R
# Correlation between USP9X expression and top 10 hub module
######################################################################

library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("figures", showWarnings = FALSE)
dir.create("figures/hubs", showWarnings = FALSE)

######################################################################
# 0) CHECK NEEDED OBJECTS / FILES
######################################################################

# Need expression + metadata from Step 00 in *this* session
if (!exists("expr_small") | !exists("meta_small")) {
  stop("Run 00_preprocess_STAR_UPDATED.R in this R session first so expr_small and meta_small exist.")
}

# Load USP9X-only data (you saved this in Step 00)
if (!exists("usp9x_df")) {
  if (!file.exists("data/USP9X_only_expression.rds")) {
    stop("data/USP9X_only_expression.rds not found. Run Step 00 first.")
  }
  usp9x_df <- readRDS("data/USP9X_only_expression.rds")
}

top10_file <- "results/hubs/STRING_top10_hubs.tsv"
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
# 1) LOAD TOP 10 HUBS + DEG TABLE AND MAP TO ENSEMBL
######################################################################

top10_hubs <- readr::read_tsv(top10_file, show_col_types = FALSE)
top10_symbols <- top10_hubs$gene

deg_res <- readr::read_tsv(deg_file, show_col_types = FALSE) %>%
  dplyr::mutate(ensembl_clean = sub("\\..*$", "", ensembl_id))

gene_map <- bitr(
  deg_res$ensembl_clean,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)

deg_mapped <- deg_res %>%
  dplyr::inner_join(gene_map, by = c("ensembl_clean" = "ENSEMBL")) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

hub_map <- deg_mapped %>%
  dplyr::filter(SYMBOL %in% top10_symbols)

ens_ids <- hub_map$ensembl_id

######################################################################
# 2) GET HUB EXPRESSION AND MODULE SCORE
######################################################################

# align metadata and expression as before
meta_small <- meta_small %>% dplyr::distinct(short_barcode, .keep_all = TRUE)
expr_small <- expr_small[, !duplicated(colnames(expr_small))]
expr_small <- expr_small[, meta_small$short_barcode]

common_ens <- intersect(ens_ids, rownames(expr_small))
if (length(common_ens) < 2) {
  stop("Fewer than 2 hub genes found in expr_small. Check mapping.")
}

hub_expr <- expr_small[common_ens, meta_small$short_barcode, drop = FALSE]

# rename rows to SYMBOLs
symbol_order <- hub_map$SYMBOL[match(rownames(hub_expr), hub_map$ensembl_id)]
rownames(hub_expr) <- as.character(symbol_order)

# z-score per gene and compute module score (mean of z-scores)
hub_expr_z <- t(scale(t(hub_expr)))
hub_module_score <- colMeans(hub_expr_z, na.rm = TRUE)

hub_module_df <- data.frame(
  short_barcode      = colnames(hub_expr),
  hub_module_score   = hub_module_score
)

######################################################################
# 3) MERGE WITH USP9X EXPRESSION
######################################################################

merged_df <- hub_module_df %>%
  dplyr::inner_join(
    usp9x_df %>% dplyr::select(short_barcode, USP9X_expr, USP9X_group),
    by = "short_barcode"
  )

cat("✅ Merged samples:", nrow(merged_df), "\n")

######################################################################
# 4) CORRELATION + SCATTERPLOT
######################################################################

cor_res <- cor.test(
  merged_df$USP9X_expr,
  merged_df$hub_module_score,
  method = "spearman"
)

cat("\n✅ Spearman correlation between USP9X and hub-module score:\n")
print(cor_res)

rho_txt <- round(cor_res$estimate, 3)
p_txt   <- signif(cor_res$p.value, 3)

p_scatter <- ggplot(merged_df,
                    aes(x = USP9X_expr, y = hub_module_score,
                        color = USP9X_group)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "USP9X expression vs hub-module activity",
    x = "USP9X expression (log-scale)",
    y = "Hub module score (mean z-score of top 10 hubs)",
    subtitle = paste0("Spearman rho = ", rho_txt, ", p = ", p_txt)
  )

ggsave("figures/hubs/USP9X_vs_hub_module_scatter.png",
       p_scatter, width = 7, height = 5, dpi = 300)

cat("\n✅ Step 08 complete. Scatterplot saved in figures/hubs/\n")

