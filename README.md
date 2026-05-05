# USP9X-Driven Coordination of Proteostasis and Stress Response Pathways in Breast Cancer

> **An integrative transcriptomic and network-omics analysis of USP9X dosage effects in TCGA-BRCA**

---

## Overview

USP9X is a deubiquitinase that stabilizes substrates involved in signaling, chromatin regulation, and cell survival. While individual USP9X targets have been studied in cell lines, the **network-level consequences of USP9X expression variation in primary breast tumors remain poorly characterized**.

This project uses bulk RNA-seq and clinical data from the **TCGA-BRCA cohort** to ask:
- Do USP9X-High and USP9X-Low tumors have distinct transcriptomic states?
- Which downstream gene programs and pathways are coordinately regulated with USP9X?
- Do USP9X-associated differentially expressed genes converge on a coherent protein interaction hub module?

---

## Key Findings

- USP9X-High and USP9X-Low breast tumors **separate clearly on PCA** (PC1: 15%, PC2: 9.6%), indicating broad transcriptomic shifts linked to USP9X dosage
- **Thousands of genes** are differentially expressed (FDR < 0.05, |log2FC| ≥ 1), with strong skew toward upregulation in USP9X-High tumors
- STRING network analysis identified a **10-gene hub module** (EP300, TRIP12, BPTF, NIPBL, DDX3X, SETD2, XPO1, UBXN7, KDM6A, HUWE1) enriched for chromatin regulation and ubiquitin-proteasome pathways
- Hub-module activity correlates strongly with continuous USP9X expression across tumors (**Spearman ρ = 0.89, p = 5.7e-78**)
- GO/KEGG enrichment reveals consistent overrepresentation of **cell cycle, protein polyubiquitination, chromosome segregation, and autophagy** pathways

---

## Repository Structure

```
USP9X-multiomics-brca/
├── analysis/                        # R analysis scripts (run in order)
│   ├── 01_preprocess_STAR.R         # Data loading, QC, USP9X group stratification
│   ├── 02_transcriptomics_plots.R   # PCA, correlation heatmap, expression plots
│   ├── 03_DEG_limma.R               # Differential expression (limma-voom)
│   ├── 04_enrichment_USP9X.R        # GO and KEGG enrichment (clusterProfiler)
│   ├── 05_STRING_network.R          # PPI network construction and hub detection
│   ├── 06_enrichment_top10_hubs.R   # Enrichment restricted to top 10 hub genes
│   ├── 07_expression_top10_hubs.R   # Hub gene expression boxplots and heatmap
│   ├── 08_USP9X_hub_module_correlation.R  # Hub-module score vs USP9X expression
│   └── 09_hub_vs_USP9X_individual_genes.R # Per-hub scatterplots vs USP9X
├── data/
│   ├── TCGA-BRCA.clinical.tsv.gz    # Clinical annotations (TCGA-BRCA)
│   ├── USP9X_only_expression.csv    # Extracted USP9X expression per sample
│   └── USP9X_only_expression.rds    # RDS version for R workflows
│   # NOTE: TCGA-BRCA.star_counts.tsv.gz (132MB) excluded from repo
│   # Download from: https://portal.gdc.cancer.gov/
├── figures/
│   ├── transcriptomics/             # PCA, correlation heatmap, USP9X distribution
│   ├── enrichment/                  # GO BP dotplot, KEGG dotplot, GSEA ridgeplot
│   ├── string/                      # PPI network plots, hub barplot, hub heatmap
│   └── hubs/                        # Hub expression boxplots, scatter plots
└── results/
    ├── DEG_limma_USP9X_High_vs_Low.tsv    # Full DEG results table
    ├── DEG_limma_USP9X_sig.tsv            # Significant DEGs only
    ├── enrichment/                         # GO BP, KEGG ORA results
    ├── pathways/                           # gProfiler results (up/down)
    └── string/                             # Hub gene tables and annotations
```

---

## Methods

| Step | Method | Tool |
|------|--------|------|
| Cohort stratification | Top/bottom 10% USP9X expression | R base |
| Quality control | Sample-sample Pearson correlation, PCA | R / ggplot2 |
| Differential expression | limma-voom, FDR < 0.05, \|log2FC\| ≥ 1 | limma |
| Pathway enrichment | ORA and GSEA | clusterProfiler |
| Network construction | STRING v11, high-confidence edges | STRINGdb |
| Hub detection | Degree centrality, top 10/50 hubs | igraph |
| Hub-module score | Mean z-score across top 10 hubs | R base |

---

## Data Access

The large STAR count matrix is excluded from this repository due to GitHub file size limits.

Download `TCGA-BRCA.star_counts.tsv.gz` from the [GDC Data Portal](https://portal.gdc.cancer.gov/):
- Project: TCGA-BRCA
- Data type: Gene Expression Quantification
- Workflow: STAR - Counts

Place the downloaded file at `data/TCGA-BRCA.star_counts.tsv.gz` before running scripts.

---

## How to Run

```r
# Install required packages (if needed)
install.packages(c("ggplot2", "dplyr", "pheatmap"))
BiocManager::install(c("limma", "edgeR", "clusterProfiler", "STRINGdb"))

# Run scripts in order
source("analysis/01_preprocess_STAR.R")
source("analysis/02_transcriptomics_plots.R")
source("analysis/03_DEG_limma.R")
# ... continue sequentially
```

---

## Selected Figures

| Figure | Description |
|--------|-------------|
| `figures/transcriptomics/PCA_plot.png` | USP9X-High vs Low tumor separation |
| `figures/transcriptomics/Volcano_limma.png` | Differential expression summary |
| `figures/string/STRING_network_top10_hubs.png` | Top 10 hub gene network |
| `figures/string/STRING_top10_hubs_barplot.png` | Hub degree ranking |
| `figures/hubs/USP9X_vs_hub_module_scatter.png` | Hub-module score vs USP9X (ρ = 0.89) |
| `figures/enrichment/GO_BP_dotplot.png` | GO Biological Process enrichment |
| `figures/enrichment/KEGG_dotplot.png` | KEGG pathway enrichment |

---

## Biological Context

USP9X-High tumors in TCGA-BRCA activate a coordinated hub module enriched for:
- **Chromatin regulation**: EP300, SETD2, KDM6A, BPTF, NIPBL
- **Ubiquitin-proteasome pathway**: TRIP12, HUWE1, UBXN7
- **Nucleocytoplasmic transport**: XPO1, DDX3X

This is consistent with USP9X's known roles in stabilizing YAP1, MCL1, and β-catenin in solid tumors, and suggests that USP9X dosage is embedded in a broader regulatory circuit rather than acting through a single substrate.

---

## Planned Extensions

- **Epigenomic layer**: H3K4me3/H3K27ac ChIP-seq and ATAC-seq at USP9X and hub loci (ENCODE/UCSC)
- **Proteomic validation**: TCGA RPPA data for USP9X hub partners (EP300, TRIP12)
- **Clinical outcomes**: Association of hub-module score with subtype, grade, and survival
- **External validation**: Reproducing hub-module correlation in METABRIC cohort

---

## Author

**Unnati Moradiya**  
MS Bioinformatics, Northeastern University  
[LinkedIn](https://www.linkedin.com/in/unnatimoradiya) | [GitHub](https://github.com/unnnattiiii)

---

## References

- Law et al. (2014) voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*
- Ritchie et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*
- Szklarczyk et al. (2021) STRING v11: protein–protein association networks. *Nucleic Acids Research*
- Yu et al. (2012) clusterProfiler: an R Package for comparing biological themes among gene clusters. *OMICS*
- Li et al. (2018) USP9X regulates centrosome duplication and promotes breast carcinogenesis. *Nature Communications*
- The Cancer Genome Atlas Network (2012) Comprehensive molecular portraits of human breast tumours. *Nature*
