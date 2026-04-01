######################################################################
# 09_hub_vs_USP9X_individual_genes.R
# Show expression of each top-10 hub gene vs USP9X
######################################################################

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("figures/hubs", showWarnings = FALSE)

######################################################################
# 0) CHECK OBJECTS / FILES
######################################################################

if (!exists("expr_small") | !exists("meta_small")) {
  stop("Run 00_preprocess_STAR_UPDATED.R in this R session first so expr_small/meta_small exist.")
}

if (!exists("usp9x_df")) {
  if (!file.exists("data/USP9X_only_expression.rds")) {
    stop("USP9X_only_expression.rds not found. Run Step 00 first.")
  }
  usp9x_df <- readRDS("data/USP9X_only_expression.rds")
}

top10_file <- "results/string/STRING_top10_hubs.tsv"
deg_file   <- "results/DEG_limma_USP9X_High_vs_Low.tsv"

if (!file.exists(top10_file)) stop("Missing ", top10_file)
if (!file.exists(deg_file))   stop("Missing ", deg_file)

######################################################################
# 1) LOAD TOP 10 HUBS + MAP SYMBOL -> ENSEMBL
######################################################################

top10_hubs    <- readr::read_tsv(top10_file, show_col_types = FALSE)
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
# 2) GET HUB EXPRESSION MATRIX (10 genes x samples)
######################################################################

meta_small <- meta_small %>% dplyr::distinct(short_barcode, .keep_all = TRUE)
expr_small <- expr_small[, !duplicated(colnames(expr_small))]
expr_small <- expr_small[, meta_small$short_barcode]

common_ens <- intersect(ens_ids, rownames(expr_small))
if (length(common_ens) < 2) stop("Fewer than 2 hubs found in expr_small.")

hub_expr <- expr_small[common_ens, meta_small$short_barcode, drop = FALSE]

# rename rows as SYMBOLs in the right order
symbol_order <- hub_map$SYMBOL[match(rownames(hub_expr), hub_map$ensembl_id)]
rownames(hub_expr) <- as.character(symbol_order)

######################################################################
# 3) LONG DATAFRAME: one row = (sample, gene)
######################################################################

hub_long <- hub_expr %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    cols      = -gene,
    names_to  = "short_barcode",
    values_to = "hub_expr"
  ) %>%
  dplyr::inner_join(
    usp9x_df %>% dplyr::select(short_barcode, USP9X_expr, USP9X_group),
    by = "short_barcode"
  )

cat("✅ Long hub-expression table rows:", nrow(hub_long), "\n")

######################################################################
# 4) FACETED SCATTER: EACH HUB vs USP9X
######################################################################

p_scatter_all <- ggplot(
  hub_long,
  aes(x = USP9X_expr, y = hub_expr, color = USP9X_group)
) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Each hub gene vs USP9X expression",
    x        = "USP9X expression (log-scale)",
    y        = "Hub gene expression (log-scale)",
    color    = "USP9X group"
  )

ggsave("figures/hubs/Scatter_each_hub_vs_USP9X.png",
       p_scatter_all, width = 10, height = 7, dpi = 300)

######################################################################
# 5) FACETED BOXPLOTS: EACH HUB BY USP9X GROUP
######################################################################

p_box_all <- ggplot(
  hub_long,
  aes(x = USP9X_group, y = hub_expr, fill = USP9X_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 0.8) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Expression of top 10 hubs in USP9X High vs Low tumors",
    x     = "",
    y     = "Hub gene expression (log-scale)"
  )

ggsave("figures/hubs/Boxplot_each_hub_by_group.png",
       p_box_all, width = 10, height = 7, dpi = 300)

cat("\n✅ Step 09 complete. Faceted scatter + boxplots in figures/hubs/\n")

