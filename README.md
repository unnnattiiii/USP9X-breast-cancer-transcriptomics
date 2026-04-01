# Integrative Multi-Omics Analysis of USP9X in Breast Cancer

## Overview
USP9X is a deubiquitinase enzyme implicated in cancer progression through regulation of key oncogenic and tumor suppressive pathways. This project performs an integrative transcriptomic analysis of USP9X expression in breast cancer using TCGA-BRCA RNA-seq data, identifying co-expressed genes, enriched pathways, and hub genes associated with USP9X activity.

**Biological Question:** What transcriptional programs are associated with high vs low USP9X expression in breast cancer, and what are the key hub genes and pathways driving this signature?

---

## Pipeline
```
TCGA-BRCA RNA-seq counts → Preprocessing & stratification → 
Differential expression (limma) → Pathway enrichment (GO/KEGG/GSEA) → 
STRING PPI network → Hub gene identification → Correlation analysis
```

---

## Repository Structure
```
USP9X-multiomics-brca/
├── analysis/
│   ├── 01_preprocess_STAR.R
│   ├── 02_transcriptomics_plots.R
│   ├── 03_DEG_limma.R
│   ├── 04_enrichment_USP9X.R
│   ├── 05_STRING_network.R
│   ├── 06_enrichment_top10_hubs_USP9X.R
│   ├── 07_expression_top10_hubs.R
│   ├── 08_USP9X_hub_module_correlation.R
│   └── 09_hub_vs_USP9X_individual_genes.R
├── data/
│   ├── TCGA-BRCA.clinical.tsv        # Clinical metadata
│   ├── USP9X_only_expression.csv     # USP9X expression values
│   └── USP9X_only_expression.rds     # USP9X expression (R object)
├── results/
│   ├── DEG_limma_USP9X_High_vs_Low.tsv
│   ├── DEG_limma_USP9X_sig.tsv
│   ├── enrichment/
│   ├── pathways/
│   └── string/
├── figures/
│   ├── enrichment/
│   ├── hubs/
│   ├── string/
│   └── transcriptomics/
├── .gitignore
└── README.md

Note: TCGA-BRCA.star_counts.tsv.gz (132MB) excluded due to GitHub file size limits.
Download from GDC portal: https://portal.gdc.cancer.gov/
Project: TCGA-BRCA | Data Type: Gene Expression Quantification | Workflow: STAR - Counts
```

---

## Methods

**Data Source:** TCGA-BRCA RNA-seq gene expression counts (STAR aligner) downloaded from GDC portal. Clinical metadata including PAM50 subtypes and survival data included.

**Preprocessing:** Raw STAR counts loaded and normalized. Samples stratified into USP9X-High and USP9X-Low groups based on median expression cutoff.

**Differential Expression:** limma-voom used for differential expression analysis between USP9X-High and USP9X-Low groups. Results filtered at FDR < 0.05.

**Pathway Enrichment:** GO (Biological Process), KEGG, and GSEA performed using clusterProfiler on significant DEGs.

**Network Analysis:** STRING PPI network constructed for hub gene identification. Top 10 hub genes identified by degree centrality.

**Correlation Analysis:** Expression correlation between USP9X and hub genes examined across the full cohort.

---

## Key Results

- Identified significant DEGs between USP9X-High and USP9X-Low breast cancer samples
- Enriched pathways characterized the biological programs associated with USP9X activity
- STRING network analysis identified top hub genes co-expressed with USP9X
- Correlation analysis confirmed USP9X-hub gene relationships at individual gene level

---

## Tools and Packages

| Tool | Version | Purpose |
|------|---------|---------|
| R | 4.x | Analysis environment |
| limma | Bioconductor | Differential expression |
| clusterProfiler | Bioconductor | Pathway enrichment |
| STRINGdb | Bioconductor | PPI network analysis |
| ggplot2 | CRAN | Visualization |
| DESeq2 | Bioconductor | Normalization support |

---

## Data Access

Raw TCGA-BRCA counts are not included in this repository due to file size constraints.

**To reproduce this analysis:**
1. Download TCGA-BRCA STAR counts from GDC portal: https://portal.gdc.cancer.gov/
2. Filter for: Project = TCGA-BRCA, Data Type = Gene Expression Quantification, Workflow Type = STAR - Counts
3. Place downloaded file in `data/` directory as `TCGA-BRCA.star_counts.tsv.gz`
4. Run scripts in order: 01 → 02 → 03 → ... → 09

---

## Author
Unnati Moradiya
MS Bioinformatics, Northeastern University
Northeastern University | Sept 2025 – Dec 2025
