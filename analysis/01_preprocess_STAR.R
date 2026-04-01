######################################################################
# 00_preprocess_STAR_UPDATED.R
# Goal: Prep High vs Low USP9X tumors + USP9X-only mini dataset
######################################################################

library(tidyverse)
library(matrixStats)

dir.create("data", showWarnings = FALSE)

### --------------------------------------------------------------
### 1. LOAD EXPRESSION MATRIX
### --------------------------------------------------------------
counts <- read_tsv("data/TCGA-BRCA.star_counts.tsv.gz") %>%
  rename(gene_id = 1)

### --------------------------------------------------------------
### 2. LOAD PHENOTYPE / CLINICAL DATA
### --------------------------------------------------------------
pheno <- read_tsv("data/TCGA-BRCA.clinical.tsv") %>%
  mutate(short_barcode = substr(sample, 1, 15))

### --------------------------------------------------------------
### 3. CLEAN COUNT MATRIX COLUMN NAMES
### --------------------------------------------------------------
expr_cols  <- colnames(counts)[-1]
expr_short <- substr(expr_cols, 1, 15)
colnames(counts) <- c("gene_id", expr_short)

### --------------------------------------------------------------
### 4. FILTER FOR PRIMARY TUMOR SAMPLES
### --------------------------------------------------------------
pheno_filtered <- pheno %>%
  filter(`sample_type.samples` == "Primary Tumor")

samples_keep <- intersect(pheno_filtered$short_barcode, expr_short)

expr <- counts %>% select(gene_id, all_of(samples_keep))

### --------------------------------------------------------------
### 5. CONVERT TO MATRIX
### --------------------------------------------------------------
expr_mat <- expr %>% column_to_rownames("gene_id") %>% as.matrix()

# keep both raw + log (log is only for plotting/QC)
expr_mat_raw <- expr_mat
expr_mat_log <- expr_mat_raw   # already log2 / normalized

### --------------------------------------------------------------
### 6. GET USP9X EXPRESSION (ROBUST)
### --------------------------------------------------------------
# correct USP9X Ensembl gene ID
usp9x_id <- "ENSG00000124486"

usp_row <- rownames(expr_mat_raw)[startsWith(rownames(expr_mat_raw), usp9x_id)]

if (length(usp_row) == 0) {
  stop("USP9X (ENSG00000124486) not found in matrix. Check rownames.")
}

usp_expr <- as.numeric(expr_mat_raw[usp_row[1], ])
names(usp_expr) <- colnames(expr_mat_raw)


### --------------------------------------------------------------
### 7. DEFINE TOP 10% HIGH vs BOTTOM 10% LOW
### --------------------------------------------------------------
q10 <- quantile(usp_expr, 0.10, na.rm = TRUE)
q90 <- quantile(usp_expr, 0.90, na.rm = TRUE)

usp_group <- ifelse(
  usp_expr >= q90, "High",
  ifelse(usp_expr <= q10, "Low", "Mid")
)

pheno_filtered <- pheno_filtered %>%
  filter(short_barcode %in% samples_keep) %>%
  mutate(USP9X_group = usp_group[short_barcode])

### --------------------------------------------------------------
### 8. KEEP ONLY HIGH + LOW GROUPS
### --------------------------------------------------------------
keep_samples <- pheno_filtered$short_barcode[pheno_filtered$USP9X_group != "Mid"]

expr_small <- expr_mat_log[, keep_samples, drop = FALSE]  # log matrix for plotting
meta_small <- pheno_filtered %>% filter(short_barcode %in% keep_samples)

### --------------------------------------------------------------
### 9. MAKE USP9X-ONLY MINI DATAFRAME + SAVE
### --------------------------------------------------------------
usp9x_df <- data.frame(
  short_barcode = names(usp_expr),
  USP9X_expr = as.numeric(usp_expr),
  USP9X_group = usp_group
) %>%
  filter(USP9X_group != "Mid")

write.csv(usp9x_df, "data/USP9X_only_expression.csv", row.names = FALSE)
saveRDS(usp9x_df, "data/USP9X_only_expression.rds")

### --------------------------------------------------------------
### 10. PRINT CHECKS (so you know it worked)
### --------------------------------------------------------------
cat("\n✅ Preprocessing complete!\n")
cat("USP9X row used:", usp_row[1], "\n")
cat("Samples retained (High + Low):", length(keep_samples), "\n\n")

cat("Group counts:\n")
print(table(usp9x_df$USP9X_group))

cat("\nHead of USP9X-only data:\n")
print(head(usp9x_df, 10))

cat("\nUSP9X expression summary:\n")
print(summary(usp9x_df$USP9X_expr))

