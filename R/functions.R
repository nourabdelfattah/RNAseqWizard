suppressPackageStartupMessages({
  library(edgeR); library(limma)
  library(ggplot2); library(ggpubr); library(ggrepel)
  library(pheatmap); library(reshape2)
  library(dplyr); library(stringr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================
# THEME + COLOR PALETTE — defined inline with fallbacks
# ============================================================

if (!exists("theme_NourMin", mode = "function")) {
  theme_NourMin <- function(base_size = 11, ...) {
    ggplot2::theme_light(base_size = base_size) +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(fill = NA, colour = "grey90", linewidth = 1),
        strip.background = ggplot2::element_rect(fill = NA, colour = NA),
        strip.text       = ggplot2::element_text(colour = "black", face = "bold"),
        axis.text        = ggplot2::element_text(colour = "black"),
        plot.title       = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(0.9)),
        legend.background = ggplot2::element_rect(fill = NA),
        plot.margin      = ggplot2::unit(c(0.3, 0.3, 0.3, 0.3), "lines")
      )
  }
}

if (!exists("Nour18")) {
  Nour18 <- c(
    petrol="#003f5c", lightpetrol="#008796", darkblue="#2f4b7c", blue="#227ed4",
    darkpurple="#665191", purple="#9e71c7", darkpink="#a05195", pink="#e082c1",
    fuchsia="#d45087", salmon="#f57689", darkred="#bd3e3e", lightred="#f95d6a",
    pumpkin="#e07b18", darkorange="#ff7c43", lightorange="#ffa600",
    yellow="#f5d256", green="#42a665", grass="#8db032"
  )
  Nour_pal <- function(palette = "main", reverse = FALSE, ...) {
    pals <- list(
      main = unname(Nour18[c("petrol","darkblue","darkpurple","darkpink",
                              "fuchsia","lightred","darkorange","lightorange")]),
      cool = unname(Nour18[c("petrol","lightpetrol","darkblue","blue",
                              "darkpurple","purple","green","grass")]),
      hot  = unname(Nour18[c("darkred","lightred","pumpkin","darkorange",
                              "lightorange","yellow","fuchsia","salmon")]),
      all  = unname(Nour18)
    )
    pal <- pals[[palette]] %||% pals$main
    if (reverse) pal <- rev(pal)
    grDevices::colorRampPalette(pal, ...)
  }
  scale_fill_Nour <- function(palette = "main", discrete = TRUE, ...) {
    pal <- Nour_pal(palette)
    if (discrete) ggplot2::discrete_scale("fill",   paste0("Nour_", palette), palette = pal, ...)
    else          ggplot2::scale_fill_gradientn(colours = pal(256), ...)
  }
  scale_color_Nour <- function(palette = "main", discrete = TRUE, ...) {
    pal <- Nour_pal(palette)
    if (discrete) ggplot2::discrete_scale("colour", paste0("Nour_", palette), palette = pal, ...)
    else          ggplot2::scale_color_gradientn(colours = pal(256), ...)
  }
}

# ============================================================
# METADATA COLOR MAP
# ============================================================

build_meta_colors <- function(meta, palette = "all") {
  pal_fn <- Nour_pal(palette)
  result <- list()
  for (col in colnames(meta)) {
    if (col %in% c("lib.size", "norm.factors")) next
    vals <- sort(unique(as.character(meta[[col]])))
    if (length(vals) >= 2 && length(vals) <= 20) {
      result[[col]] <- setNames(pal_fn(length(vals)), vals)
    }
  }
  result
}

# Rebuild colors after relabeling — carries over user-edited colors for groups
# that still exist; assigns fresh palette colors for any new group values.
rebuild_meta_colors <- function(meta, existing = NULL, palette = "all") {
  new_cols <- build_meta_colors(meta, palette)
  if (is.null(existing)) return(new_cols)
  for (col in names(new_cols)) {
    if (!col %in% names(existing)) next
    for (grp in names(new_cols[[col]])) {
      if (grp %in% names(existing[[col]]))
        new_cols[[col]][[grp]] <- existing[[col]][[grp]]
    }
  }
  new_cols
}

# ============================================================
# DATA LOADING
# ============================================================

.detect_gene_col <- function(cnts) {
  cands <- c("gene_name","geneName","gene_id","Gene","genes","symbol","SYMBOL","Name")
  found <- cands[cands %in% colnames(cnts)]
  if (length(found)) found[1] else NULL
}

load_and_filter <- function(counts_input, meta_input,
                             gene_types      = c("protein_coding","lincRNA"),
                             filter_by_type  = TRUE,
                             min_cpm         = 0.5,
                             min_samples     = 2) {
  cnts <- if (is.character(counts_input))
    read.csv(counts_input, check.names = FALSE, stringsAsFactors = FALSE)
  else counts_input

  meta <- if (is.character(meta_input))
    read.csv(meta_input, as.is = TRUE, row.names = 1, stringsAsFactors = FALSE)
  else meta_input

  gene_col   <- .detect_gene_col(cnts)
  type_col   <- if ("gene_type" %in% colnames(cnts)) "gene_type" else NULL
  count_cols <- intersect(colnames(cnts), rownames(meta))

  if (length(count_cols) == 0)
    stop("No matching sample names between count columns and metadata rownames.\n",
         "Count cols: ",  paste(head(colnames(cnts), 6), collapse=", "), "\n",
         "Meta rows:  ",  paste(head(rownames(meta), 6), collapse=", "))

  if (filter_by_type && !is.null(type_col))
    cnts <- cnts[cnts[[type_col]] %in% gene_types, ]

  mat   <- cnts[, count_cols, drop = FALSE]
  genes <- if (!is.null(gene_col)) cnts[[gene_col]] else as.character(seq_len(nrow(cnts)))

  valid <- !is.na(genes) & nchar(trimws(genes)) > 0
  mat   <- mat[valid, ]; genes <- genes[valid]

  dge  <- DGEList(counts  = mat,
                  genes   = data.frame(genes = genes, stringsAsFactors = FALSE),
                  samples = meta[count_cols, , drop = FALSE])
  keep <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)
  dge
}

make_log2cpm <- function(dge) {
  gene_names <- dge$genes$genes
  idx        <- !duplicated(gene_names)
  mat        <- cpm(dge[idx, ], log = TRUE)
  rownames(mat) <- gene_names[idx]
  mat
}

# ============================================================
# QC PLOTS
# ============================================================

plot_library_size <- function(dge, color_by = NULL, palette = "main", color_map = NULL,
                               sort_by_size = FALSE) {
  df <- data.frame(sample   = rownames(dge$samples),
                   lib_size = dge$samples$lib.size / 1e6,
                   dge$samples, stringsAsFactors = FALSE)
  df$sample <- factor(df$sample,
                       levels = if (sort_by_size) df$sample[order(df$lib_size)]
                                else rownames(dge$samples))
  p <- ggplot(df, aes(x = sample, y = lib_size))
  if (!is.null(color_by) && color_by %in% colnames(df))
    p <- p + aes(fill = .data[[color_by]])
  else
    p <- p + aes(fill = sample)
  fill_scale <- if (!is.null(color_map) && !is.null(color_by))
    ggplot2::scale_fill_manual(values = color_map) else scale_fill_Nour(palette = palette)
  p + geom_col() + coord_flip() +
    fill_scale +
    labs(x = NULL, y = "Library size (M reads)") +
    theme_NourMin(base_size = 11) +
    theme(legend.position = "none")
}

plot_density <- function(dge, color_by = NULL, palette = "main", color_map = NULL) {
  mat <- cpm(dge, log = TRUE)
  df  <- reshape2::melt(mat, varnames = c("gene","sample"), value.name = "log2cpm")
  df$sample <- as.character(df$sample)
  df <- merge(df, cbind(sample = rownames(dge$samples), dge$samples, stringsAsFactors = FALSE),
              by = "sample")
  p <- ggplot(df, aes(x = log2cpm, group = sample))
  if (!is.null(color_by) && color_by %in% colnames(df))
    p <- p + aes(color = .data[[color_by]]) +
      (if (!is.null(color_map)) ggplot2::scale_color_manual(values = color_map)
       else scale_color_Nour(palette = palette))
  p + geom_density() +
    labs(x = "Log2 CPM", y = "Density") +
    theme_NourMin(base_size = 11)
}

# ============================================================
# PCA
# ============================================================

plot_pca <- function(log2cpm, meta, color_by = NULL, shape_by = NULL,
                     n_genes = 2500, label_samples = FALSE, palette = "main",
                     color_map = NULL) {
  vars <- apply(log2cpm, 1, var, na.rm = TRUE)
  top  <- names(sort(vars, decreasing = TRUE))[1:min(n_genes, sum(!is.na(vars)))]
  mat  <- t(scale(t(log2cpm[top, ]))); mat[is.na(mat)] <- 0

  pca <- prcomp(t(mat), scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                   sample = rownames(pca$x), stringsAsFactors = FALSE)
  df <- merge(df, cbind(sample = rownames(meta), meta, stringsAsFactors = FALSE),
              by = "sample", all.x = TRUE)

  p <- ggplot(df, aes(PC1, PC2))
  if (!is.null(color_by) && color_by %in% colnames(df)) {
    p <- p + aes(color = .data[[color_by]]) +
      (if (!is.null(color_map)) ggplot2::scale_color_manual(values = color_map)
       else scale_color_Nour(palette = palette))
  }
  if (!is.null(shape_by) && shape_by %in% colnames(df))
    p <- p + aes(shape = .data[[shape_by]])

  p <- p + geom_point(size = 3) +
    labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
    theme_NourMin(base_size = 12)

  if (label_samples)
    p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE)
  p
}

# ============================================================
# DIFFERENTIAL EXPRESSION
# ============================================================

run_DE <- function(dge, group_col, group_a, group_b, fdr = 0.05, lfc = 1) {
  keep   <- dge$samples[[group_col]] %in% c(group_a, group_b)
  sub    <- dge[, keep, keep.lib.sizes = FALSE]
  group  <- factor(sub$samples[[group_col]], levels = c(group_a, group_b))
  design <- model.matrix(~0 + group)
  ca     <- make.names(group_a); cb <- make.names(group_b)
  colnames(design) <- c(ca, cb)

  sub  <- estimateDisp(sub, design, robust = TRUE)
  fit  <- glmQLFit(sub, design, robust = TRUE)
  con  <- makeContrasts(contrasts = paste0(ca, "-", cb), levels = design)
  res  <- glmQLFTest(fit, contrast = con)

  tab  <- topTags(res, n = nrow(sub))$table
  tab  <- tab[!duplicated(tab$genes), ]
  rownames(tab) <- tab$genes

  list(
    all  = tab,
    sig  = tab[!is.na(tab$FDR) & tab$FDR < fdr  & abs(tab$logFC) > lfc, ],
    up   = tab[!is.na(tab$FDR) & tab$FDR < fdr  & tab$logFC >  lfc, ],
    down = tab[!is.na(tab$FDR) & tab$FDR < fdr  & tab$logFC < -lfc, ]
  )
}

run_DE_allpairs <- function(dge, group_col, groups = NULL, fdr = 0.05, lfc = 1) {
  if (is.null(groups))
    groups <- sort(unique(as.character(dge$samples[[group_col]])))
  pairs <- combn(groups, 2, simplify = FALSE)
  results <- lapply(pairs, function(p) {
    tryCatch(run_DE(dge, group_col, p[1], p[2], fdr = fdr, lfc = lfc),
             error = function(e) { message("Pair ", p[1], "/", p[2], ": ", conditionMessage(e)); NULL })
  })
  names(results) <- sapply(pairs, function(p) paste0(p[1], " vs ", p[2]))
  Filter(Negate(is.null), results)
}

run_DE_omnibus <- function(dge, group_col, groups = NULL, fdr = 0.05) {
  if (!is.null(groups))
    dge <- dge[, dge$samples[[group_col]] %in% groups, keep.lib.sizes = FALSE]

  group  <- factor(dge$samples[[group_col]])
  design <- model.matrix(~group)
  dge    <- estimateDisp(dge, design, robust = TRUE)
  fit    <- glmQLFit(dge, design, robust = TRUE)
  res    <- glmQLFTest(fit, coef = 2:ncol(design))

  tab    <- topTags(res, n = nrow(dge))$table
  tab    <- tab[!duplicated(tab$genes), ]
  rownames(tab) <- tab$genes

  list(all = tab, sig = tab[!is.na(tab$FDR) & tab$FDR < fdr, ])
}

# ============================================================
# VOLCANO / MA / OMNIBUS PLOTS
# ============================================================

plot_volcano <- function(de, fdr = 0.05, lfc = 1, n_label = 10) {
  df     <- de$all
  df$sig <- dplyr::case_when(
    !is.na(df$FDR) & df$FDR < fdr & df$logFC >  lfc ~ "Up",
    !is.na(df$FDR) & df$FDR < fdr & df$logFC < -lfc ~ "Down",
    TRUE ~ "NS"
  )
  df$nlp <- pmin(-log10(df$FDR + 1e-300), 300)
  top    <- rownames(head(df[df$sig != "NS", ][order(df[df$sig != "NS","FDR"]), ], n_label))
  df$label <- ifelse(rownames(df) %in% top, rownames(df), "")

  ggplot(df, aes(logFC, nlp, color = sig, label = label)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c(Up = unname(Nour18["lightred"]),
                                  Down = unname(Nour18["darkblue"]), NS = "grey70")) +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-lfc, lfc),   linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 15, show.legend = FALSE) +
    labs(x = "log2 Fold Change", y = "-log10(FDR)", color = "") +
    theme_NourMin(base_size = 12)
}

plot_ma <- function(de, fdr = 0.05, lfc = 1) {
  df     <- de$all
  df$sig <- dplyr::case_when(
    !is.na(df$FDR) & df$FDR < fdr & df$logFC >  lfc ~ "Up",
    !is.na(df$FDR) & df$FDR < fdr & df$logFC < -lfc ~ "Down",
    TRUE ~ "NS"
  )
  ggplot(df, aes(logCPM, logFC, color = sig)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c(Up = unname(Nour18["lightred"]),
                                  Down = unname(Nour18["darkblue"]), NS = "grey70")) +
    geom_hline(yintercept = c(-lfc, 0, lfc),
               linetype = c("dashed","solid","dashed"), color = "grey50") +
    labs(x = "Average log2 CPM", y = "log2 Fold Change", color = "") +
    theme_NourMin(base_size = 12)
}

plot_omni <- function(de_omni, fdr = 0.05) {
  df      <- de_omni$all
  df$sig  <- !is.na(df$FDR) & df$FDR < fdr
  df$nlp  <- pmin(-log10(df$FDR + 1e-300), 300)
  top     <- rownames(head(df[df$sig, ][order(df[df$sig,"FDR"]), ], 10))
  df$label <- ifelse(rownames(df) %in% top, rownames(df), "")

  ggplot(df, aes(logCPM, nlp, color = sig, label = label)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("TRUE"  = unname(Nour18["lightred"]),
                                  "FALSE" = "grey70"),
                       labels = c("TRUE" = "Significant", "FALSE" = "NS")) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 15, show.legend = FALSE) +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "grey50") +
    labs(x = "Average log2 CPM", y = "-log10(FDR)", color = "") +
    theme_NourMin(base_size = 12)
}

# ============================================================
# HEATMAP
# ============================================================

plot_heatmap <- function(log2cpm, genes, meta, ann_cols = NULL,
                          cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
                          ann_colors = NULL) {
  genes <- intersect(genes, rownames(log2cpm))
  if (length(genes) == 0) stop("No matching genes in expression matrix.")

  mat       <- log2cpm[genes, , drop = FALSE]; mode(mat) <- "numeric"
  cell_h    <- max(6, min(14, floor(400 / length(genes))))
  hm_colors <- colorRampPalette(
    c(unname(Nour18["darkblue"]), "white", unname(Nour18["lightred"])))(200)

  ann <- if (!is.null(ann_cols) && length(ann_cols) > 0)
    meta[colnames(mat), ann_cols, drop = FALSE] else NULL

  pheatmap::pheatmap(mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                     scale = scale, annotation_col = ann,
                     annotation_colors = ann_colors,
                     border_color = "black", fontsize = 9, cellwidth = 20,
                     cellheight = cell_h, color = hm_colors, silent = TRUE)
}

# ============================================================
# BOXPLOTS
# ============================================================

parse_genes <- function(text) {
  g <- unlist(strsplit(text, "[,;\n\r\t ]+"))
  unique(trimws(g[nchar(trimws(g)) > 0]))
}

plot_boxplots <- function(log2cpm, meta, genes, x_by,
                           color_by    = x_by,
                           facet_by    = NULL,
                           test        = "t.test",
                           posthoc     = "bonferroni",
                           comp_style  = "pairwise",
                           ref_group   = NULL,
                           p_display   = "p.format",
                           palette     = "main",
                           color_map   = NULL) {

  genes <- intersect(genes, rownames(log2cpm))
  if (length(genes) == 0) stop("None of the specified genes found.")

  mat <- log2cpm[genes, , drop = FALSE]
  df  <- reshape2::melt(t(mat), varnames = c("sample","gene"), value.name = "log2cpm")
  df$sample <- as.character(df$sample); df$gene <- as.character(df$gene)
  df <- merge(df, cbind(sample = rownames(meta), meta, stringsAsFactors = FALSE),
              by = "sample", all.x = TRUE)

  groups   <- sort(unique(as.character(df[[x_by]])))
  n_groups <- length(groups)

  cmps <- if (n_groups == 2) {
    list(groups)
  } else if (comp_style == "pairwise") {
    combn(groups, 2, simplify = FALSE)
  } else {
    rg <- ref_group %||% groups[1]
    lapply(setdiff(groups, rg), function(g) c(rg, g))
  }

  pairwise_m <- if (test %in% c("t.test","anova"))    "t.test"      else "wilcox.test"
  global_m   <- if (test %in% c("t.test","anova"))    "anova"       else "kruskal.test"

  p <- ggplot(df, aes(x = .data[[x_by]], y = log2cpm, fill = .data[[color_by]])) +
    geom_boxplot(width = 0.5, lwd = 0.3, colour = "black") +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.8,
                 position = position_dodge(0.8)) +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    scale_y_continuous(expand = expansion(mult = 0.28, add = 0)) +
    (if (!is.null(color_map)) ggplot2::scale_fill_manual(values = color_map)
     else scale_fill_Nour(palette = palette)) +
    labs(y = "Log2 CPM", x = NULL) +
    theme_NourMin(base_size = 10) +
    theme(axis.text.x  = element_text(angle = 25, hjust = 1),
          plot.margin   = ggplot2::margin(5, 10, 5, 12, "mm"))

  has_facet  <- !is.null(facet_by) && facet_by != "None" && facet_by %in% colnames(df)
  multi_gene <- length(genes) > 1
  if      (has_facet && multi_gene) p <- p + facet_grid(reformulate(facet_by, "gene"), scales = "free_y")
  else if (multi_gene)              p <- p + facet_wrap(~gene, scales = "free_y")
  else if (has_facet)               p <- p + facet_wrap(reformulate(facet_by), scales = "free_y")

  ph <- if (posthoc == "none" || n_groups == 2) "none" else posthoc

  if (n_groups >= 2 && length(cmps) > 0) {
    # Faceted plots: stat_compare_means only draws stats for the first facet column.
    # Use rstatix per-panel stats whenever a facet_by column is active.
    if (has_facet && requireNamespace("rstatix", quietly = TRUE)) {
      facet_vars   <- c(if (multi_gene) "gene", facet_by)
      formula_rsx  <- as.formula(paste("log2cpm ~", x_by))

      stat_df <- tryCatch({
        grp_df  <- dplyr::group_by(df, dplyr::across(dplyr::all_of(facet_vars)))
        test_fn <- if (pairwise_m == "t.test") rstatix::pairwise_t_test
                   else rstatix::pairwise_wilcox_test
        stat_r  <- test_fn(grp_df, formula_rsx, comparisons = cmps, p.adjust.method = ph)
        # facet.by was added to rstatix::add_xy_position in a later version;
        # fall back gracefully if the installed version doesn't support it.
        tryCatch(
          rstatix::add_xy_position(stat_r, x = x_by, data = df, formula = formula_rsx,
                                    scales = "free", fun = "max", step.increase = 0.12,
                                    facet.by = facet_vars),
          error = function(e2)
            rstatix::add_xy_position(stat_r, x = x_by, data = df, formula = formula_rsx,
                                      scales = "free", fun = "max", step.increase = 0.12)
        )
      }, error = function(e) {
        message("Per-panel stats error: ", conditionMessage(e)); NULL
      })

      if (!is.null(stat_df) && nrow(stat_df) > 0) {
        lbl <- if (p_display == "p.signif") "p.adj.signif" else "p.adj"
        p   <- p + ggpubr::stat_pvalue_manual(stat_df, label = lbl,
                                               tip.length = 0.01, size = 2.5,
                                               hide.ns = TRUE)
      }
    } else {
      # Non-faceted (or rstatix unavailable): stat_compare_means works correctly
      if (n_groups > 2)
        p <- p + ggpubr::stat_compare_means(method = global_m, label.y.npc = 0.97, size = 2.5)
      p <- p + ggpubr::stat_compare_means(comparisons = cmps, method = pairwise_m,
                                           p.adjust.method = ph, label = p_display,
                                           size = 2.5, vjust = -0.3)
    }
  }
  p
}

# ============================================================
# NORMALITY CHECK
# ============================================================

check_normality <- function(log2cpm, gene, meta, x_by) {
  gene <- intersect(gene, rownames(log2cpm))[1]
  if (is.na(gene)) return(NULL)
  vals   <- as.numeric(log2cpm[gene, ])
  df     <- data.frame(val = vals, group = as.character(meta[[x_by]]))
  groups <- sort(unique(df$group))
  lapply(setNames(groups, groups), function(g) {
    x <- df$val[df$group == g]; n <- length(x)
    if (n < 3) return(list(n=n, W=NA, p=NA, normal=NA, note="n < 3 — skipped"))
    sw <- tryCatch(shapiro.test(x), error = function(e) NULL)
    if (is.null(sw)) return(list(n=n, W=NA, p=NA, normal=NA, note="error"))
    list(n=n, W=round(sw$statistic[[1]],4), p=round(sw$p.value,4),
         normal=sw$p.value > 0.05, note="")
  })
}

recommend_test <- function(norm_results, n_groups) {
  all_ok <- all(sapply(norm_results, function(r) is.na(r$normal) || isTRUE(r$normal)))
  if (n_groups == 2) { if (all_ok) "t.test" else "wilcox.test" }
  else               { if (all_ok) "anova"  else "kruskal.test" }
}

# ============================================================
# SPECIES + MOUSE → HUMAN CONVERSION
# ============================================================

detect_species <- function(genes, n = 200) {
  g <- head(genes[!grepl("^[0-9]|Rik$|^Gm[0-9]", genes)], n)
  if (mean(g == toupper(g), na.rm = TRUE) > 0.8) "human" else "mouse"
}

convert_mouse_human <- function(mouse_genes) {
  tryCatch({
    archive    <- "https://dec2021.archive.ensembl.org"
    human_mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = archive)
    mouse_mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = archive)
    map <- biomaRt::getLDS(attributes = "mgi_symbol", filters = "mgi_symbol",
                           values = mouse_genes, mart = mouse_mart,
                           attributesL = "hgnc_symbol", martL = human_mart, uniqueRows = TRUE)
    map <- map[!duplicated(map[,1]) & nchar(map[,2]) > 0, ]
    setNames(map[,2], map[,1])
  }, error = function(e) { warning("biomaRt: ", conditionMessage(e)); NULL })
}

apply_human_symbols <- function(log2cpm, conversion) {
  common <- intersect(rownames(log2cpm), names(conversion))
  if (!length(common)) stop("No overlap with conversion map.")
  mat <- log2cpm[common, , drop = FALSE]
  rownames(mat) <- conversion[common]
  mat[!duplicated(rownames(mat)), ]
}

# ============================================================
# MSIGDB GENE SETS
# ============================================================

MSIGDB_COLLECTIONS <- list(
  "Hallmark (H)"                 = list(category = "H",  subcategory = NULL),
  "KEGG (C2)"                    = list(category = "C2", subcategory = "CP:KEGG_LEGACY",
                                         subcategory_alt = "CP:KEGG"),
  "Reactome (C2)"                = list(category = "C2", subcategory = "CP:REACTOME"),
  "WikiPathways (C2)"            = list(category = "C2", subcategory = "CP:WIKIPATHWAYS"),
  "GO: Biological Process (C5)"  = list(category = "C5", subcategory = "GO:BP"),
  "GO: Molecular Function (C5)"  = list(category = "C5", subcategory = "GO:MF"),
  "Immunologic signatures (C7)"  = list(category = "C7", subcategory = NULL),
  "Oncogenic signatures (C6)"    = list(category = "C6", subcategory = NULL)
)

load_msigdb_sets <- function(species = "Homo sapiens", collection = "Hallmark (H)") {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Install msigdbr: BiocManager::install('msigdbr')")
  spec <- MSIGDB_COLLECTIONS[[collection]]
  if (is.null(spec)) stop("Unknown collection: ", collection)

  # Always use species= (not db_species=); the species= arg returns gene_symbol
  # in the target species' own symbol format (e.g. "Abcc4" for mouse) which is
  # what the ranked DE list uses. db_species= selects a different gene-set
  # database and is intentionally avoided here.
  .fetch <- function(subcat) {
    args <- list(species = species, category = spec$category)
    if (!is.null(subcat)) args$subcategory <- subcat
    do.call(msigdbr::msigdbr, args)
  }

  mdf <- tryCatch(.fetch(spec$subcategory), error = function(e) {
    if (!is.null(spec$subcategory_alt))
      tryCatch(.fetch(spec$subcategory_alt),
               error = function(e2) stop(conditionMessage(e)))
    else stop(conditionMessage(e))
  })

  if (nrow(mdf) == 0)
    stop("No gene sets returned for '", collection, "' / species='", species,
         "'. Try a different collection or species.")
  split(mdf$gene_symbol, mdf$gs_name)
}

# ============================================================
# GSEA
# ============================================================

run_gsea <- function(de_all, gene_sets, pval = 0.05, n_top = 15) {
  gene_list <- sort(setNames(de_all$logFC, rownames(de_all)), decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list)) & !is.na(gene_list)]

  pathways <- if (is.character(gene_sets)) fgsea::gmtPathways(gene_sets) else gene_sets

  # Overlap diagnostic — stop early with an informative message
  all_pw_genes <- unique(unlist(pathways))
  n_overlap    <- length(intersect(names(gene_list), all_pw_genes))
  if (n_overlap < 10)
    stop("Only ", n_overlap, " of your genes match the pathway database (",
         length(all_pw_genes), " unique symbols in ", length(pathways), " sets). ",
         "Likely cause: mouse symbols vs human MSigDB. ",
         "Enable 'Mouse -> Human (biomaRt)' conversion or select 'Mus musculus' as species.")

  # Use fgseaMultilevel when available (fgsea >= 1.14), fall back to classic fgsea
  fgsea_res <- tryCatch(
    fgsea::fgseaMultilevel(pathways = pathways, stats = gene_list,
                            minSize = 5, maxSize = 600, eps = 0),
    error = function(e)
      fgsea::fgsea(pathways = pathways, stats = gene_list,
                   minSize = 5, maxSize = 600, nperm = 1000)
  )

  # fgseaMultilevel returns a data.table with a list column (leadingEdge).
  # Extract non-list columns first, then flatten leadingEdge to a string.
  non_list_cols <- names(fgsea_res)[!sapply(fgsea_res, is.list)]
  res <- as.data.frame(fgsea_res)[, non_list_cols, drop = FALSE]
  res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = ";")
  res <- res[order(-res$NES), , drop = FALSE]

  if (nrow(res) == 0)
    return(list(results = res, plot = NULL,
                n_overlap = n_overlap, n_sets = length(pathways)))

  res$pathway <- stringr::str_replace_all(
    res$pathway, "^(GO|HALLMARK|REACTOME|KEGG_LEGACY|KEGG|WP)_", "")
  res$pathway <- stringr::str_replace_all(res$pathway, "_", " ")
  res$dir     <- ifelse(res$NES > 0, "Up-regulated", "Down-regulated")
  res$sig     <- !is.na(res$padj) & res$padj < pval

  # top n_top up + top n_top down by absolute NES
  plot_df <- rbind(
    head(res[res$NES > 0, ], n_top),
    head(res[res$NES < 0, ], n_top)
  )

  p <- ggplot(plot_df, aes(reorder(pathway, NES), NES)) +
    geom_segment(aes(xend = pathway, y = 0, yend = NES), color = "grey60") +
    geom_point(aes(fill = dir, alpha = sig), size = 3, shape = 21, stroke = 1) +
    scale_fill_manual(values = c("Up-regulated"   = unname(Nour18["lightred"]),
                                 "Down-regulated" = unname(Nour18["darkblue"]))) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35),
                       labels = c("TRUE" = paste0("padj < ", pval), "FALSE" = "NS"),
                       name = "") +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    coord_flip() +
    labs(x = NULL, y = "Normalized Enrichment Score", fill = "") +
    theme_NourMin(base_size = 10) +
    theme(axis.text.y = element_text(face = "bold"))

  list(results = res, plot = p, n_overlap = n_overlap, n_sets = length(pathways))
}

# ============================================================
# ORA
# ============================================================

run_ora <- function(de, species = "mouse", ont = "BP", pval = 0.05, qval = 0.2) {
  require(clusterProfiler)
  orgdb      <- if (species == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
  org_prefix <- if (species == "mouse") "mmu" else "hsa"

  to_entrez <- function(syms)
    tryCatch(clusterProfiler::bitr(syms, "SYMBOL", "ENTREZID", OrgDb = orgdb)$ENTREZID,
             error = function(e) character(0))

  bg <- to_entrez(rownames(de$all))
  up <- to_entrez(rownames(de$up))
  dn <- to_entrez(rownames(de$down))
  out <- list()

  if (length(up) >= 5) {
    out$up_go   <- clusterProfiler::enrichGO(up, OrgDb=orgdb, ont=ont, universe=bg,
                                              pvalueCutoff=pval, qvalueCutoff=qval, readable=TRUE)
    out$up_kegg <- clusterProfiler::enrichKEGG(up, organism=org_prefix, universe=bg, pvalueCutoff=pval)
  }
  if (length(dn) >= 5) {
    out$dn_go   <- clusterProfiler::enrichGO(dn, OrgDb=orgdb, ont=ont, universe=bg,
                                              pvalueCutoff=pval, qvalueCutoff=qval, readable=TRUE)
    out$dn_kegg <- clusterProfiler::enrichKEGG(dn, organism=org_prefix, universe=bg, pvalueCutoff=pval)
  }
  out
}

# ============================================================
# GSVA
# ============================================================

run_gsva <- function(log2cpm, gene_sets, meta, group_col, group_a, group_b,
                      fdr = 0.05, lfc = log2(1.5)) {
  require(GSVA); require(limma)

  if (group_a == group_b)
    stop("group_a and group_b must be different.")

  pathways <- if (is.character(gene_sets)) fgsea::gmtPathways(gene_sets) else gene_sets

  scores <- GSVA::gsva(as.matrix(log2cpm), pathways, min.sz = 10, max.sz = 500,
                        verbose = FALSE)

  keep   <- meta[[group_col]] %in% c(group_a, group_b)
  valid  <- intersect(rownames(meta)[keep], colnames(scores))
  if (length(valid) < 2)
    stop("Fewer than 2 samples found for the selected groups in the score matrix.")

  scores <- scores[, valid, drop = FALSE]
  meta_s <- meta[valid, , drop = FALSE]
  group  <- factor(meta_s[[group_col]], levels = c(group_a, group_b))

  design <- model.matrix(~0 + group)
  ca     <- make.names(group_a); cb <- make.names(group_b)
  colnames(design) <- c(ca, cb)

  fit  <- limma::lmFit(scores, design)
  cont <- limma::makeContrasts(contrasts = paste0(ca, "-", cb), levels = design)
  fit  <- limma::contrasts.fit(fit, cont)
  fit  <- limma::eBayes(fit)

  list(scores = scores,
       all    = limma::topTable(fit, number = Inf),
       de     = limma::topTable(fit, number = Inf, p.value = fdr, lfc = lfc),
       design = design)
}

# ============================================================
# DT HELPER — collapse list columns
# ============================================================

.dt_safe <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  for (col in colnames(df))
    if (is.list(df[[col]]))
      df[[col]] <- vapply(df[[col]], function(x) paste(x, collapse=";"), character(1))
  df
}
