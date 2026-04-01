######################################################################
# 05_STRING_network_USP9X.R
# STRING / PPI network after limma DEGs (USP9X High vs Low)
######################################################################

# ---- Libraries ----
library(tidyverse)
library(STRINGdb)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(ggraph)
library(ggplot2)

# ---- Output folders ----
dir.create("results", showWarnings = FALSE)
dir.create("results/string", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("figures/string", showWarnings = FALSE)

######################################################################
# 1) LOAD DEG RESULTS
######################################################################

deg_file_sig  <- "results/DEG_limma_USP9X_sig.tsv"
deg_file_full <- "results/DEG_limma_USP9X_High_vs_Low.tsv"

if (file.exists(deg_file_sig)) {
  deg_res <- read_tsv(deg_file_sig, show_col_types = FALSE)
  message("Loaded significant DEGs file.")
} else {
  deg_res <- read_tsv(deg_file_full, show_col_types = FALSE)
  message("Loaded full DEG file, filtering now.")
  
  deg_res <- deg_res %>%
    mutate(ensembl_clean = sub("\\..*$", "", ensembl_id)) %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1)
}

# ensure we have clean Ensembl IDs
if (!"ensembl_clean" %in% colnames(deg_res)) {
  deg_res <- deg_res %>%
    mutate(ensembl_clean = sub("\\..*$", "", ensembl_id))
}

cat("\n✅ Significant DEGs loaded:", nrow(deg_res), "\n")

######################################################################
# 2) MAP ENSEMBL -> SYMBOL
######################################################################

gene_map <- bitr(
  deg_res$ensembl_clean,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)

deg_mapped <- deg_res %>%
  inner_join(gene_map, by = c("ensembl_clean"="ENSEMBL")) %>%
  distinct(SYMBOL, .keep_all = TRUE)

cat("✅ Mapped to SYMBOL:", nrow(deg_mapped), "\n")

######################################################################
# 3) PICK TOP GENES (avoid huge network)
######################################################################

# fix weird column type issue (Rle) and go to plain data.frame
deg_mapped2 <- deg_mapped %>%
  mutate(
    adj.P.Val = as.numeric(adj.P.Val),
    logFC     = as.numeric(logFC)
  ) %>%
  as.data.frame()

TOP_N <- 300              # change to 200 if too dense
TOP_N <- min(TOP_N, nrow(deg_mapped2))

# order by FDR then |logFC| using base R
ord <- order(deg_mapped2$adj.P.Val, -abs(deg_mapped2$logFC))
top_genes <- deg_mapped2[ord[1:TOP_N], , drop = FALSE]

cat("✅ Using top genes for STRING:", nrow(top_genes), "\n")

######################################################################
# 4) QUERY STRING
######################################################################

string_db <- STRINGdb$new(
  version="11.5",
  species=9606,
  score_threshold=400
)

string_input <- data.frame(gene = top_genes$SYMBOL)

mapped_string <- string_db$map(
  string_input,
  "gene",
  removeUnmappedRows = TRUE
)

cat("✅ STRING mapped genes:", nrow(mapped_string), "\n")

ppi <- string_db$get_interactions(mapped_string$STRING_id)
cat("✅ Interactions pulled:", nrow(ppi), "\n")

######################################################################
# 5) BUILD NETWORK + HUB GENES  (patched: no dplyr::select)
######################################################################

# keep only from/to columns for edges
ppi_edges <- ppi[, c("from", "to")]

# build igraph object
g <- graph_from_data_frame(ppi_edges, directed = FALSE)

# node table: STRING id + gene symbol
nodes_df <- mapped_string[, c("STRING_id", "gene")]
nodes_df <- unique(nodes_df)

# add labels to graph
V(g)$label <- nodes_df$gene[match(V(g)$name, nodes_df$STRING_id)]

# hub stats
hub_table <- data.frame(
  STRING_id   = V(g)$name,
  gene        = V(g)$label,
  degree      = degree(g),
  betweenness = betweenness(g)
)
hub_table <- hub_table[order(-hub_table$degree), ]

write_tsv(hub_table, "results/string/STRING_hub_genes.tsv")

cat("\nTop hub genes:\n")
print(head(hub_table, 10))

# ---- Top 50 hub subnetwork (cleaner for slides) ----
top_hubs <- hub_table %>%
  dplyr::filter(degree > 0) %>%
  dplyr::slice(1:50)

sub_nodes <- top_hubs$STRING_id

g_sub <- induced_subgraph(g, vids = V(g)[name %in% sub_nodes])

png("figures/string/STRING_network_top50_hubs.png",
    width = 2000, height = 2000, res = 200)

plot(
  g_sub,
  vertex.label     = V(g_sub)$label,
  vertex.size      = 6,
  vertex.label.cex = 0.7,
  edge.color       = "grey80"
)

dev.off()

# ---- Top 10 hub subnetwork (very focused view) ----
top_hubs10 <- hub_table %>%
  dplyr::filter(degree > 0) %>%
  dplyr::slice(1:10)

sub_nodes10 <- top_hubs10$STRING_id

g_sub10 <- induced_subgraph(g, vids = V(g)[name %in% sub_nodes10])

png("figures/string/STRING_network_top10_hubs.png",
    width = 2000, height = 2000, res = 200)

plot(
  g_sub10,
  vertex.label     = V(g_sub10)$label,
  vertex.size      = 8,      # slightly bigger so they show nicely
  vertex.label.cex = 0.9,
  edge.color       = "grey60"
)

dev.off()



######################################################################
# 6) PLOT FULL NETWORK (all STRING nodes)
######################################################################

# store node degree as an attribute so ggraph can see it
V(g)$degree_attr <- degree(g)

set.seed(42)

p_net <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(size = degree_attr), color = "pink") +
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle(paste0("STRING Network: USP9X High vs Low (Top ", TOP_N, " DEGs)"))

ggsave("figures/string/STRING_network_topDEGs.png",
       p_net, width = 9, height = 7, dpi = 300)

print(p_net)

######################################################################
# 7) ANNOTATE & EXPORT TOP HUB GENES
######################################################################

# take the same hub_table we already made
# merge in logFC and adj.P.Val from DE results
hub_annot <- hub_table %>%
  dplyr::left_join(
    deg_mapped2 %>% dplyr::select(SYMBOL, logFC, adj.P.Val),
    by = c("gene" = "SYMBOL")
  )

# save full hub table (all nodes in STRING network)
write_tsv(hub_annot, "results/string/STRING_hub_genes_annotated.tsv")

# TOP 50 and TOP 10 versions, nicely annotated for write-up
top50_annot <- hub_annot %>% dplyr::slice(1:50)
top10_annot <- hub_annot %>% dplyr::slice(1:10)

write_tsv(top50_annot, "results/string/STRING_top50_hubs_annotated.tsv")
write_tsv(top10_annot, "results/string/STRING_top10_hubs_annotated.tsv")

# (optional) quick barplot of top 10 by degree
top10_annot$gene <- factor(top10_annot$gene,
                           levels = top10_annot$gene[order(top10_annot$degree)])

p_bar <- ggplot(top10_annot, aes(x = gene, y = degree)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 10 hub proteins in USP9X network",
    x = "Gene", y = "Degree (number of interactions)"
  ) +
  theme_minimal()

ggsave("figures/string/STRING_top10_hubs_barplot.png",
       p_bar, width = 6, height = 4, dpi = 300)


