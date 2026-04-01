######################################################################
# 04_enrichment_USP9X.R
# Enrichment after limma DEG (USP9X High vs Low)
# Input: results/DEG_limma_USP9X_High_vs_Low.tsv
######################################################################

library(tidyverse)          # dplyr + ggplot2, avoids select() conflict
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("results/enrichment", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("figures/enrichment", showWarnings = FALSE)

######################################################################
### 0) LOAD DEG RESULTS
######################################################################

deg_res <- read_tsv("results/DEG_limma_USP9X_High_vs_Low.tsv")

######################################################################
### 1) CLEAN IDs + SIGNIFICANT DEG LIST
######################################################################

deg_res <- deg_res %>%
  dplyr::mutate(ensembl_clean = sub("\\..*$", "", ensembl_id))

deg_sig <- deg_res %>%
  dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 1)

cat("\n✅ Significant DEGs:", nrow(deg_sig), "\n")
write_tsv(deg_sig, "results/DEG_limma_USP9X_sig.tsv")

######################################################################
### 2) MAP ENSEMBL -> ENTREZ
######################################################################

gene_map <- bitr(
  deg_sig$ensembl_clean,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

deg_sig_mapped <- deg_sig %>%
  dplyr::inner_join(gene_map, by = c("ensembl_clean" = "ENSEMBL"))

entrez_sig <- unique(deg_sig_mapped$ENTREZID)

cat("✅ Mapped DEGs to Entrez:", length(entrez_sig), "\n")

######################################################################
### 3) GO Biological Process ORA
######################################################################

ego_bp <- enrichGO(
  gene          = entrez_sig,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

write_tsv(as.data.frame(ego_bp), "results/enrichment/GO_BP_ORA.tsv")

png("figures/enrichment/GO_BP_dotplot.png", width=1800, height=1200, res=200)
print(dotplot(ego_bp, showCategory=15) + ggtitle("GO BP Enrichment (ORA)"))
dev.off()

######################################################################
### 4) KEGG ORA
######################################################################

ekegg <- enrichKEGG(
  gene          = entrez_sig,
  organism      = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

write_tsv(as.data.frame(ekegg), "results/enrichment/KEGG_ORA.tsv")

png("figures/enrichment/KEGG_dotplot.png", width=1800, height=1200, res=200)
print(dotplot(ekegg, showCategory=15) + ggtitle("KEGG Enrichment (ORA)"))
dev.off()

######################################################################
### 5) GSEA (GO BP) using ranked logFC
######################################################################

# Map ALL genes (not only significant) for ranking
gene_map_all <- bitr(
  deg_res$ensembl_clean,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

gene_rank <- deg_res %>%
  dplyr::select(ensembl_clean, logFC) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(gene_map_all, by = c("ensembl_clean" = "ENSEMBL")) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarise(logFC = mean(logFC), .groups="drop")

rank_vec <- gene_rank$logFC
names(rank_vec) <- gene_rank$ENTREZID
rank_vec <- sort(rank_vec, decreasing = TRUE)

gsea_bp <- gseGO(
  geneList      = rank_vec,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  verbose       = FALSE
)

write_tsv(as.data.frame(gsea_bp), "results/enrichment/GO_BP_GSEA.tsv")

png("figures/enrichment/GO_BP_GSEA_ridgeplot.png", width=1800, height=1200, res=200)
print(ridgeplot(gsea_bp, showCategory=15) + ggtitle("GO BP (GSEA)"))
dev.off()

######################################################################
### 6) PRINT TOP TERMS
######################################################################

cat("\nTop GO BP ORA terms:\n")
print(head(as.data.frame(ego_bp)[, c("ID","Description","p.adjust")], 10))

cat("\nTop KEGG ORA pathways:\n")
print(head(as.data.frame(ekegg)[, c("ID","Description","p.adjust")], 10))

cat("\nTop GO BP GSEA terms:\n")
print(head(as.data.frame(gsea_bp)[, c("ID","Description","p.adjust","NES")], 10))

cat("\n✅ Step 4 enrichment complete!\n")
cat("Tables saved in results/enrichment/\n")
cat("Plots saved in figures/enrichment/\n")

