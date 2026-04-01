######################################################################
# 03_DEG_limma.R
# Differential Expression: USP9X High vs Low (using log-scale input)
######################################################################

library(tidyverse)
library(limma)
library(matrixStats)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures/transcriptomics", showWarnings = FALSE)

# ---- make sure Step 0 objects exist ----
if (!exists("expr_small") | !exists("meta_small")) {
  stop("Run Step 0 first — expr_small / meta_small missing.")
}

######################################################################
### 1) CLEAN / ALIGN
######################################################################

# remove duplicate samples if any + align order
meta_small <- meta_small %>% distinct(short_barcode, .keep_all = TRUE)
expr_small <- expr_small[, !duplicated(colnames(expr_small))]
expr_small <- expr_small[, meta_small$short_barcode]

group <- factor(meta_small$USP9X_group, levels=c("Low","High"))

######################################################################
### 2) LIMMA DEG
######################################################################

design <- model.matrix(~ group)
fit <- lmFit(expr_small, design)
fit <- eBayes(fit)

deg_res <- topTable(fit, coef="groupHigh", number=Inf, sort.by="P")

# add Ensembl without version for readability
deg_res$ensembl_id <- rownames(deg_res)
deg_res$ensembl_clean <- sub("\\..*$", "", deg_res$ensembl_id)

# save results
write_tsv(deg_res, "results/DEG_limma_USP9X_High_vs_Low.tsv")

cat("\n✅ LIMMA DEG done!\n")
cat("Top genes:\n")
print(head(deg_res[, c("ensembl_clean","logFC","P.Value","adj.P.Val")], 10))

######################################################################
### 3) VOLCANO (real FDR-based)
######################################################################

volcano_df <- deg_res %>%
  mutate(
    neglog10FDR = -log10(adj.P.Val),
    Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No")
  )

p_volcano <- ggplot(volcano_df, aes(logFC, neglog10FDR, color=Significant)) +
  geom_point(alpha=0.5, size=1.2) +
  theme_minimal(base_size=14) +
  labs(
    title="Volcano Plot: USP9X High vs Low (limma)",
    x="log2FC (High vs Low)",
    y="-log10(FDR)"
  ) +
  scale_color_manual(values=c("No"="gray70", "Yes"="red")) +
  coord_cartesian(ylim=c(0, 20)) +        # 👈 makes it readable
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed")

ggsave("figures/transcriptomics/Volcano_limma.png",
       p_volcano, width=7, height=5, dpi=300)
print(p_volcano)

cat("\n✅ Step 3 complete. Results in results/ and volcano saved.\n")

