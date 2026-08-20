# ============================================================
# FIGURE LEGENDS & METHODS COMPILER
# ============================================================

legend_libsize <- function(color_by = NULL) {
  col_str <- if (!is.null(color_by)) paste0(" coloured by ", color_by) else ""
  paste0(
    "Library Size. Bar plot showing total mapped read counts (millions) per sample",
    col_str, ". Samples are ordered by library size. ",
    "All samples should have comparable library sizes; large discrepancies may indicate ",
    "quality issues or batch effects."
  )
}

legend_density <- function(color_by = NULL) {
  col_str <- if (!is.null(color_by)) paste0(", coloured by ", color_by) else ""
  paste0(
    "Expression Density. Kernel density estimate of log₂ CPM values across all retained genes ",
    "per sample", col_str, ". Curves should have similar shapes and peaks; ",
    "strongly shifted or bimodal curves may indicate outlier samples."
  )
}

legend_pca <- function(color_by = NULL, shape_by = NULL, n_genes = 2500) {
  dims_str <- paste0("top ", format(n_genes, big.mark = ","), " most variable genes")
  col_str  <- if (!is.null(color_by) && color_by != "None")
    paste0(" Samples are coloured by ", color_by, ".") else ""
  shp_str  <- if (!is.null(shape_by) && shape_by != "None")
    paste0(" Shape encodes ", shape_by, ".") else ""
  paste0(
    "Principal Component Analysis (PCA). Scatter plot of PC1 vs PC2 computed from the ",
    dims_str, " (log₂ CPM, row-scaled).", col_str, shp_str,
    " Percentages in axis labels indicate variance explained by each component."
  )
}

legend_volcano <- function(group_a, group_b, fdr = 0.05, lfc = 1) {
  paste0(
    "Volcano Plot: ", group_a, " vs ", group_b, ". ",
    "Each point represents a gene. The x-axis shows log₂ fold-change (", group_a,
    " relative to ", group_b, ") and the y-axis shows −log₂₀(FDR). ",
    "Red points are significantly up-regulated and blue points significantly down-regulated ",
    "(FDR < ", fdr, ", |log₂FC| > ", lfc, "). ",
    "Gene labels indicate the top hits by FDR."
  )
}

legend_ma <- function(group_a, group_b, fdr = 0.05, lfc = 1) {
  paste0(
    "MA Plot: ", group_a, " vs ", group_b, ". ",
    "The x-axis shows average log₂ CPM across all samples; the y-axis shows ",
    "log₂ fold-change. Dashed horizontal lines mark the |log₂FC| = ", lfc,
    " threshold. Coloured points meet FDR < ", fdr, "."
  )
}

legend_omni <- function(groups, fdr = 0.05) {
  grp_str <- paste(groups, collapse = ", ")
  paste0(
    "Omnibus Differential Expression. MA-style plot of average log₂ CPM vs ",
    "−log₂₀(FDR) from a quasi-likelihood F-test across all groups (",
    grp_str, "). Red points are significant at FDR < ", fdr, ". ",
    "Gene labels indicate the top hits. Use pairwise comparisons to determine direction."
  )
}

legend_heatmap <- function(n_genes, ann_cols = NULL, scale = "row") {
  ann_str <- if (!is.null(ann_cols) && length(ann_cols) > 0)
    paste0(" Column annotations show: ", paste(ann_cols, collapse = ", "), ".") else ""
  scale_str <- switch(scale,
    row    = "Colour scale reflects row z-scores (genes scaled across samples).",
    column = "Colour scale reflects column z-scores (samples scaled across genes).",
    none   = "Colour scale reflects unscaled log₂ CPM values.",
    "")
  paste0(
    "Heatmap of ", n_genes, " genes. ",
    "Rows are genes; columns are samples.", ann_str, " ",
    scale_str,
    " Blue–white–red gradient; rows clustered by Euclidean distance / complete linkage."
  )
}

legend_boxplot <- function(genes, x_by, test = "t.test", posthoc = "bonferroni",
                            comp_style = "pairwise") {
  gene_str  <- if (length(genes) == 1) paste0("gene ", genes[1])
               else paste0(length(genes), " genes")
  test_str  <- switch(test,
    t.test       = "Student’s t-test",
    wilcox.test  = "Wilcoxon rank-sum test",
    anova        = "one-way ANOVA",
    kruskal.test = "Kruskal–Wallis test",
    test)
  ph_str    <- if (posthoc != "none")
    paste0(" Post-hoc p-values corrected by ", posthoc, " method.") else ""
  comp_str  <- if (comp_style == "pairwise") "all pairwise comparisons" else "reference group comparisons"
  paste0(
    "Boxplot of ", gene_str, " expression (log₂ CPM) grouped by ", x_by, ". ",
    "Individual data points are overlaid as dots. ",
    "Statistical comparisons shown for ", comp_str, " using ", test_str, ".", ph_str
  )
}

legend_gsea <- function(collection = NULL, species = NULL, pval = 0.05, n_top = 15) {
  src_str <- if (!is.null(collection)) paste0("MSigDB ", collection) else "custom GMT file"
  sp_str  <- if (!is.null(species)) paste0(" (", species, ")") else ""
  paste0(
    "GSEA lollipop plot. Gene Set Enrichment Analysis was performed using gene sets from ",
    src_str, sp_str, " with fgsea (min.sz=5, max.sz=600). ",
    "Each point represents a significantly enriched gene set (adjusted p < ", pval, "). ",
    "Normalised Enrichment Score (NES) > 0 indicates up-regulation; NES < 0 down-regulation. ",
    "Top ", n_top, " pathways per direction are shown."
  )
}

legend_gsva <- function(collection = NULL, species = NULL, fdr = 0.05) {
  src_str <- if (!is.null(collection)) paste0("MSigDB ", collection) else "custom GMT file"
  sp_str  <- if (!is.null(species)) paste0(" (", species, ")") else ""
  paste0(
    "GSVA heatmap. Gene Set Variation Analysis scores were computed with GSVA ",
    "using gene sets from ", src_str, sp_str, ". ",
    "Differential pathway activity between groups was tested with limma; ",
    "heatmap shows significant pathways (FDR < ", fdr, "), scaled by row. ",
    "Blue–white–red gradient."
  )
}

# ============================================================
# METHODS COMPILER
# ============================================================

compile_methods <- function(params) {
  sections <- list()

  # Data loading
  sections$loading <- paste0(
    "Raw count data were imported and intersected with sample metadata. ",
    if (isTRUE(params$filter_by_type) && !is.null(params$gene_types))
      paste0("Genes were restricted to the following biotypes: ",
             paste(params$gene_types, collapse = ", "), ". ")
    else "",
    "Low-expression genes were removed, retaining those with >",
    params$min_cpm %||% 0.5, " CPM in at least ", params$min_samples %||% 2,
    " samples. ",
    "Libraries were normalised using TMM (trimmed mean of M-values) as implemented in edgeR."
  )

  # PCA
  if (!is.null(params$pca_ngenes)) {
    sections$pca <- paste0(
      "For visualisation, principal component analysis (PCA) was performed on the ",
      format(params$pca_ngenes, big.mark = ","),
      " most variable genes (log₂ CPM, row-scaled)."
    )
  }

  # DE
  if (!is.null(params$de_mode)) {
    de_str <- switch(params$de_mode,
      "2group"  = paste0("between ", params$de_a, " and ", params$de_b),
      "allpairs"= "across all pairwise group combinations",
      "omnibus" = "across all groups via omnibus F-test",
      "")
    sections$de <- paste0(
      "Differential expression analysis was performed ", de_str,
      " using the edgeR glmQLFTest pipeline. Dispersions were estimated with ",
      "estimateDisp() (robust=TRUE). Genes were considered differentially expressed at ",
      "FDR < ", params$de_fdr %||% 0.05,
      if (!is.null(params$de_lfc) && params$de_mode != "omnibus")
        paste0(" and |log₂FC| > ", params$de_lfc)
      else "",
      "."
    )
  }

  # Pathway
  if (!is.null(params$path_method)) {
    gene_src <- if (isTRUE(params$path_use_msigdb) && !is.null(params$path_collection))
      paste0("MSigDB ", params$path_collection, " (msigdbr, ", params$path_species, ")")
    else "a custom GMT file"

    sections$pathway <- switch(params$path_method,
      gsea = paste0(
        "Gene Set Enrichment Analysis (GSEA) was performed with fgsea using gene sets from ",
        gene_src, ". Genes were pre-ranked by log₂ fold-change. ",
        "Pathways with adjusted p < ", params$path_pval %||% 0.05, " are reported."
      ),
      gsva = paste0(
        "Gene Set Variation Analysis (GSVA) scores were computed using gene sets from ",
        gene_src, ". Differential pathway activity was tested with limma::eBayes(). ",
        "Significant pathways (FDR < ", params$path_pval %||% 0.05, ") are shown."
      ),
      ora = paste0(
        "Over-representation analysis (ORA) was performed with clusterProfiler using ",
        "GO (", params$ora_ont %||% "BP", ") and KEGG databases for ",
        params$path_species %||% "mouse", ". Significance threshold: p < ",
        params$path_pval %||% 0.05, "."
      ),
      ""
    )
  }

  # Stats for boxplots
  if (!is.null(params$bp_test)) {
    test_str <- switch(params$bp_test,
      t.test       = "Student’s t-test",
      wilcox.test  = "Wilcoxon rank-sum test",
      anova        = "one-way ANOVA",
      kruskal.test = "Kruskal–Wallis test",
      params$bp_test)
    ph_str <- if (!is.null(params$bp_posthoc) && params$bp_posthoc != "none")
      paste0(" Post-hoc p-values were corrected using the ",
             params$bp_posthoc, " method.") else ""
    sections$boxplot <- paste0(
      "Gene expression comparisons in boxplots used ", test_str, ".", ph_str,
      " Normality was assessed per group using the Shapiro–Wilk test."
    )
  }

  paste(unlist(sections), collapse = "\n\n")
}

# ============================================================
# PPTX BUILDER
# ============================================================

build_bulk_pptx <- function(figures_list, methods_text = NULL,
                              params = list(), out_path = "BulkSeq_Report.pptx") {
  if (!requireNamespace("officer",  quietly = TRUE))
    stop("Install officer: install.packages('officer')")
  if (!requireNamespace("ggplot2",  quietly = TRUE))
    stop("Install ggplot2.")

  prs <- officer::read_pptx()

  # Title slide
  prs <- officer::add_slide(prs, layout = "Title Slide", master = "Office Theme")
  prs <- officer::ph_with(prs, value = "BulkSeq Wizard — Analysis Report",
                           location = officer::ph_location_type("ctrTitle"))
  prs <- officer::ph_with(prs, value = format(Sys.time(), "%B %d, %Y"),
                           location = officer::ph_location_type("subTitle"))

  tmp_dir <- tempfile("bulkseq_pptx_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  for (fig in figures_list) {
    name    <- fig$name    %||% "Figure"
    legend  <- fig$legend  %||% ""
    plot_obj <- fig$plot

    if (is.null(plot_obj)) next

    img_path <- file.path(tmp_dir, paste0(gsub("[^A-Za-z0-9_]", "_", name), ".png"))
    tryCatch({
      if (inherits(plot_obj, "pheatmap")) {
        png(img_path, width = 2400, height = 2000, res = 200)
        grid::grid.newpage(); grid::grid.draw(plot_obj$gtable)
        dev.off()
      } else {
        ggplot2::ggsave(img_path, plot_obj, width = 9, height = 6, dpi = 200)
      }
    }, error = function(e) {
      message("PPTX: could not render '", name, "': ", conditionMessage(e))
      return(NULL)
    })

    if (!file.exists(img_path)) next

    prs <- officer::add_slide(prs, layout = "Title and Content", master = "Office Theme")
    prs <- officer::ph_with(prs, value = name,
                             location = officer::ph_location_type("title"))
    prs <- officer::ph_with(prs, value = officer::external_img(img_path, width = 8, height = 5),
                             location = officer::ph_location(left=0.5, top=1.2, width=8, height=5))

    if (nchar(trimws(legend)) > 0) {
      prs <- officer::ph_with(prs,
               value = officer::fpar(officer::ftext(legend,
                         prop = officer::fp_text(font.size = 8))),
               location = officer::ph_location(left=0.5, top=6.35, width=8.5, height=0.9))
    }
  }

  # Methods slide(s)
  if (!is.null(methods_text) && nchar(trimws(methods_text)) > 0) {
    prs <- officer::add_slide(prs, layout = "Title and Content", master = "Office Theme")
    prs <- officer::ph_with(prs, value = "Methods",
                             location = officer::ph_location_type("title"))
    prs <- officer::ph_with(prs,
             value = officer::fpar(officer::ftext(methods_text,
                       prop = officer::fp_text(font.size = 9))),
             location = officer::ph_location_type("body"))
  }

  print(prs, target = out_path)
  invisible(out_path)
}
