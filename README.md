# RNAseq Wizard

A Shiny app for interactive bulk RNA-seq analysis: differential expression, QC, heatmaps, boxplots, GSEA, GSVA, and ORA.

## Quick start

Install the required packages (see [Requirements](#requirements) below), then run directly from GitHub — no download needed:

```r
shiny::runGitHub("RNAseqWizard", "nourabdelfattah")
```

## Requirements

R >= 4.2 with the following packages:

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyWidgets", "DT",
  "ggplot2", "ggrepel", "ggpubr", "patchwork",
  "dplyr", "stringr", "pheatmap", "RColorBrewer",
  "openxlsx", "rmarkdown"
))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "limma", "fgsea", "GSVA",
  "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db",
  "biomaRt"
))

install.packages("msigdbr")
```

## Input files

| File | Format | Notes |
|------|--------|-------|
| Counts | CSV | Genes x Samples raw counts; first column = gene IDs |
| Metadata | CSV | Samples x Variables; first column = sample IDs matching count column names |

## Running the app

```r
shiny::runApp("app.R")
```

Or open `app.R` in RStudio and click **Run App**.

## Workflow

1. **Data** — upload counts + metadata, filter low-expressed genes
2. **QC** — library size bar chart, density plots, PCA
3. **DE** — two-group or all-pairs edgeR differential expression; volcano, MA, heatmap
4. **Boxplots** — CPM boxplots for genes of interest with pairwise stats
5. **Pathways** — GSEA (fgsea), GSVA, or ORA; supports MSigDB collections or custom GMT files

## Species

Mouse and human are supported. When running mouse data against human MSigDB collections, enable **Mouse -> Human (biomaRt)** conversion in the Pathways tab.
