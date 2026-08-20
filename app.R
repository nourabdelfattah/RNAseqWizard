suppressPackageStartupMessages({
  library(shiny)
  library(shinycssloaders)
  library(DT)
  library(ggplot2)
  library(grid)
  library(dplyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

source("R/functions.R")
source("R/legends.R")

# ============================================================
# HELPERS
# ============================================================

spinner_col <- "#2c7bb6"

legend_box <- function(id) {
  tagList(
    tags$hr(),
    tags$small(tags$b("Figure legend")),
    textAreaInput(id, NULL, rows = 4, resize = "vertical",
                  placeholder = "Auto-generated legend — edit as needed.")
  )
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', sans-serif; font-size: 13px; }
    .well { background-color: #f8f9fa; border: 1px solid #dee2e6; }
    h5 { font-weight: 600; margin-top: 10px; }
    .btn-primary { background-color: #2c7bb6; border-color: #2c7bb6; }
    .warn-box { background:#fff3cd; border:1px solid #ffc107; border-radius:4px;
                padding:8px 12px; margin-bottom:10px; font-size:12px; }
    .norm-result { font-size:11px; background:#f8f9fa; border:1px solid #dee2e6;
                   padding:6px; border-radius:4px; white-space:pre-wrap; }
  "))),

  navbarPage(
    title = "BulkSeq Wizard",
    id    = "main_tabs",

    # ----------------------------------------------------------------
    # 1. LOAD DATA
    # ----------------------------------------------------------------
    tabPanel("Load Data", value = "tab_load",
      sidebarLayout(
        sidebarPanel(width = 3,
          h5("Input Files"),
          fileInput("counts_file", "Count Matrix (CSV)", accept = ".csv"),
          fileInput("meta_file",   "Sample Metadata (CSV)", accept = ".csv"),
          tags$p(tags$em("Metadata: rows = sample IDs, columns = variables.")),
          tags$hr(),
          h5("Filtering"),
          checkboxInput("filter_by_type", "Filter by gene_type column (if present)", TRUE),
          conditionalPanel("input.filter_by_type",
            checkboxGroupInput("gene_types", NULL,
              choices  = c("protein_coding","lincRNA","pseudogene","antisense"),
              selected = c("protein_coding","lincRNA"))
          ),
          numericInput("min_cpm",     "Min CPM threshold",     value = 0.5, min = 0, step = 0.1),
          numericInput("min_samples", "Min samples above CPM", value = 2,   min = 1),
          tags$hr(),
          actionButton("load_btn", "Load & Filter Data", class = "btn-primary",
                       style = "width:100%"),
          tags$br(), tags$br(),
          verbatimTextOutput("filter_summary")
        ),
        mainPanel(width = 9,
          # Sample-mismatch warning (hidden until data loads)
          uiOutput("sample_warning"),
          fluidRow(
            column(6, h5("Count Matrix (first 5 rows)"), DTOutput("counts_preview")),
            column(6, h5("Metadata (all rows)"),         DTOutput("meta_preview"))
          ),
          tags$hr(),

          # ---- Sample Labels & Column Order ----
          h5("Sample Labels & Column Order"),
          fluidRow(
            column(3,
              radioButtons("relabel_mode", "Label source",
                choices = c("Original IDs"        = "none",
                            "From metadata column" = "col",
                            "Concat two columns"   = "concat"),
                selected = "none")
            ),
            column(4, uiOutput("relabel_opts_ui")),
            column(3,
              uiOutput("relabel_sort_ui"),
              radioButtons("relabel_sortdir", NULL,
                           choices = c("Ascending" = "asc", "Descending" = "desc"),
                           selected = "asc", inline = TRUE)
            ),
            column(2,
              tags$br(),
              actionButton("apply_relabel", "Apply",
                           class = "btn-primary", style = "width:100%; margin-bottom:6px"),
              actionButton("reset_relabel", "Reset", style = "width:100%")
            )
          ),
          div(class = "norm-result", verbatimTextOutput("relabel_preview")),
          tags$hr(),

          # ---- Metadata Color Assignments ----
          h5("Metadata Color Assignments"),
          tags$p(tags$em(
            "Colors are auto-assigned from the full Nour palette and used consistently
             across all plots. Edit and click Apply to override.")),
          uiOutput("color_editor_ui"),
          conditionalPanel("output.has_colors",
            actionButton("apply_colors", "Apply Colors",
                         class = "btn-primary", style = "margin-top:6px")
          )
        )
      )
    ),

    # ----------------------------------------------------------------
    # 2. QC
    # ----------------------------------------------------------------
    tabPanel("QC", value = "tab_qc",
      sidebarLayout(
        sidebarPanel(width = 3,
          selectInput("qc_color",  "Color/fill by", choices = NULL),
          checkboxInput("qc_sort_libsize", "Sort bar chart by library size", value = FALSE),
          numericInput("qc_height", "Plot height (px)", value = 420, min = 200, step = 20),
          tags$hr(),
          downloadButton("dl_libsize", "Download Library Size PDF", style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_density", "Download Density PDF",      style = "width:100%"),
          legend_box("legend_libsize"),
          legend_box("legend_density")
        ),
        mainPanel(width = 9,
          fluidRow(
            column(6, h5("Library Size"),
                   plotOutput("lib_size_plot") %>% withSpinner(color = spinner_col)),
            column(6, h5("Expression Density (log₂ CPM)"),
                   plotOutput("density_plot")  %>% withSpinner(color = spinner_col))
          )
        )
      )
    ),

    # ----------------------------------------------------------------
    # 3. PCA
    # ----------------------------------------------------------------
    tabPanel("PCA", value = "tab_pca",
      sidebarLayout(
        sidebarPanel(width = 3,
          selectInput("pca_color",  "Color by",  choices = NULL),
          selectInput("pca_shape",  "Shape by",  choices = NULL),
          numericInput("pca_ngenes","Top variable genes", value = 2500, min = 100, step = 100),
          checkboxInput("pca_label","Label samples", value = FALSE),
          actionButton("run_pca", "Run PCA", class = "btn-primary", style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_pca", "Download PDF", style = "width:100%"),
          legend_box("legend_pca")
        ),
        mainPanel(width = 9,
          plotOutput("pca_plot", height = "520px") %>% withSpinner(color = spinner_col)
        )
      )
    ),

    # ----------------------------------------------------------------
    # 4. DE ANALYSIS
    # ----------------------------------------------------------------
    tabPanel("DE Analysis", value = "tab_de",
      sidebarLayout(
        sidebarPanel(width = 3,
          h5("Comparison Mode"),
          radioButtons("de_mode", NULL,
            choices = c("Two groups"     = "2group",
                        "All pairs"      = "allpairs",
                        "Omnibus F-test" = "omnibus"),
            selected = "2group"),
          tags$hr(),
          selectInput("de_col", "Group column", choices = NULL),
          conditionalPanel("input.de_mode == '2group'",
            selectInput("de_a", "Group A (numerator)",   choices = NULL),
            selectInput("de_b", "Group B (denominator)", choices = NULL)
          ),
          conditionalPanel("input.de_mode == 'allpairs' || input.de_mode == 'omnibus'",
            checkboxGroupInput("de_groups_sel", "Groups to include", choices = NULL)
          ),
          tags$hr(),
          numericInput("de_fdr", "FDR cutoff",      value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput("de_lfc", "|log₂FC| cutoff", value = 1,    min = 0, step = 0.25),
          tags$hr(),
          actionButton("run_de", "Run DE Analysis", class = "btn-primary", style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_de_table", "Download DE Table", style = "width:100%"),
          tags$br(), tags$br(),
          verbatimTextOutput("de_summary"),
          conditionalPanel("input.de_mode == 'allpairs'",
            tags$hr(),
            selectInput("de_pair_view", "View pair", choices = NULL)
          ),
          legend_box("legend_volcano"),
          legend_box("legend_ma")
        ),
        mainPanel(width = 9,
          tabsetPanel(id = "de_inner",
            tabPanel("Table",
              DTOutput("de_table") %>% withSpinner(color = spinner_col)),
            tabPanel("Volcano",
              fluidRow(column(2, downloadButton("dl_volcano","PDF")), column(10)),
              tags$br(),
              plotOutput("volcano_plot", height = "520px") %>% withSpinner(color = spinner_col)),
            tabPanel("MA Plot",
              fluidRow(column(2, downloadButton("dl_ma","PDF")), column(10)),
              tags$br(),
              plotOutput("ma_plot", height = "520px") %>% withSpinner(color = spinner_col)),
            tabPanel("Omnibus Plot",
              fluidRow(column(2, downloadButton("dl_omni","PDF")), column(10)),
              tags$br(),
              plotOutput("omni_plot", height = "520px") %>% withSpinner(color = spinner_col))
          )
        )
      )
    ),

    # ----------------------------------------------------------------
    # 5. HEATMAP
    # ----------------------------------------------------------------
    tabPanel("Heatmap", value = "tab_hm",
      sidebarLayout(
        sidebarPanel(width = 3,
          radioButtons("hm_source", "Gene source",
            choices = c("Custom list" = "custom", "DE results" = "de"), selected = "custom"),
          conditionalPanel("input.hm_source == 'custom'",
            textAreaInput("hm_genes", "Genes (one per line or comma-separated)", rows = 8)
          ),
          conditionalPanel("input.hm_source == 'de'",
            checkboxInput("hm_de_sig", "Only FDR-significant genes", TRUE),
            numericInput("hm_n_genes", "Max genes per comparison", value = 50, min = 2, step = 5),
            uiOutput("hm_pairs_ui"),
            radioButtons("hm_de_gene_rule", "Gene rule (multi-pair)",
                         choices = c("Union" = "union", "Intersection" = "intersection"),
                         selected = "union", inline = TRUE)
          ),
          tags$hr(),
          checkboxGroupInput("hm_ann_cols", "Annotation columns", choices = NULL),
          checkboxInput("hm_cluster_rows", "Cluster rows",    TRUE),
          checkboxInput("hm_cluster_cols", "Cluster columns", FALSE),
          selectInput("hm_scale", "Scale by",
                      choices = c("row","column","none"), selected = "row"),
          tags$hr(),
          actionButton("run_hm", "Generate Heatmap", class = "btn-primary", style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_heatmap", "Download PDF", style = "width:100%"),
          legend_box("legend_heatmap")
        ),
        mainPanel(width = 9,
          div(style = "overflow-y: auto; max-height: 80vh;",
            plotOutput("heatmap_plot", height = "auto") %>% withSpinner(color = spinner_col)
          )
        )
      )
    ),

    # ----------------------------------------------------------------
    # 6. BOXPLOTS
    # ----------------------------------------------------------------
    tabPanel("Boxplots", value = "tab_bp",
      sidebarLayout(
        sidebarPanel(width = 3,
          textAreaInput("bp_genes", "Genes (one per line or comma-separated)", rows = 5),
          selectInput("bp_x",     "X-axis (grouping)", choices = NULL),
          selectInput("bp_facet", "Facet by",          choices = NULL),
          selectInput("bp_color", "Color/fill by",     choices = NULL),
          tags$hr(),
          h5("Statistics"),
          actionButton("bp_norm_check", "Check normality", style = "width:100%"),
          tags$br(), tags$br(),
          div(class = "norm-result", verbatimTextOutput("bp_norm_result")),
          tags$br(),
          selectInput("bp_test", "Statistical test",
                      choices = c("t.test","wilcox.test","anova","kruskal.test"),
                      selected = "t.test"),
          selectInput("bp_posthoc", "Post-hoc correction",
                      choices = c("bonferroni","BH","none"), selected = "bonferroni"),
          radioButtons("bp_comp_style", "Comparison style",
                       choices = c("All pairwise" = "pairwise", "vs reference" = "reference"),
                       selected = "pairwise", inline = TRUE),
          conditionalPanel("input.bp_comp_style == 'reference'",
            selectInput("bp_ref_group", "Reference group", choices = NULL)
          ),
          selectInput("bp_p_display", "P-value display",
                      choices = c("p.format","p.signif"), selected = "p.format"),
          tags$hr(),
          fluidRow(
            column(6, numericInput("bp_width",     "Width (in)",  value = 7, min = 2, step = 0.5)),
            column(6, numericInput("bp_height_in", "Height (in)", value = 5, min = 2, step = 0.5))
          ),
          actionButton("run_bp", "Generate", class = "btn-primary", style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_boxplot", "Download PDF", style = "width:100%"),
          legend_box("legend_boxplot")
        ),
        mainPanel(width = 9,
          plotOutput("boxplot_plot", height = "520px") %>% withSpinner(color = spinner_col)
        )
      )
    ),

    # ----------------------------------------------------------------
    # 7. PATHWAYS
    # ----------------------------------------------------------------
    tabPanel("Pathways", value = "tab_path",
      sidebarLayout(
        sidebarPanel(width = 3,
          radioButtons("path_method", "Method",
            choices = c("GSEA (fgsea)" = "gsea", "GSVA" = "gsva",
                        "ORA (clusterProfiler)" = "ora")),
          tags$hr(),
          radioButtons("path_geneset_src", "Gene Sets",
            choices = c("MSigDB (msigdbr)" = "msigdb", "Upload GMT file" = "gmt"),
            selected = "msigdb"),
          conditionalPanel("input.path_geneset_src == 'msigdb'",
            selectInput("path_species_msig", "Species",
                        choices = c("Homo sapiens","Mus musculus"), selected = "Mus musculus"),
            selectInput("path_collection", "Collection",
                        choices = names(MSIGDB_COLLECTIONS), selected = "Hallmark (H)")
          ),
          conditionalPanel("input.path_geneset_src == 'gmt'",
            fileInput("gmt_file", "GMT file", accept = c(".gmt",".txt"))
          ),
          tags$hr(),
          radioButtons("path_species", "Symbol conversion",
            choices = c("Mouse" = "mouse", "Human" = "human",
                        "Mouse → Human (biomaRt)" = "mouse2human"),
            selected = "mouse"),
          conditionalPanel("input.path_method == 'ora'",
            selectInput("ora_ont", "GO Ontology",
                        choices = c("BP","MF","CC","ALL"), selected = "BP")
          ),
          # Pair selector for GSEA/ORA when allpairs DE was run
          uiOutput("path_pair_ui"),
          conditionalPanel("input.path_method == 'gsva'",
            tags$hr(),
            selectInput("gsva_de_col", "Group column", choices = NULL),
            selectInput("gsva_de_a",   "Group A",      choices = NULL),
            selectInput("gsva_de_b",   "Group B",      choices = NULL)
          ),
          numericInput("path_pval", "p/padj cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
          conditionalPanel("input.path_method == 'gsea'",
            numericInput("gsea_n_top", "Top pathways per direction", value = 15, min = 3, step = 5)
          ),
          conditionalPanel("input.path_method == 'gsva'",
            numericInput("gsva_n_top", "Top/bottom pathways in heatmap", value = 25, min = 3, step = 5),
            tags$hr(),
            h5("Heatmap size"),
            numericInput("gsva_hm_height", "Display height (px)", value = 640, min = 200, step = 50),
            fluidRow(
              column(6, numericInput("gsva_dl_w", "PDF width (in)",  value = 10, min = 3, step = 1)),
              column(6, numericInput("gsva_dl_h", "PDF height (in)", value = 8,  min = 3, step = 1))
            )
          ),
          tags$hr(),
          actionButton("run_path", "Run Analysis", class = "btn-primary", style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_path_results", "Download Results",  style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_path_plot",    "Download Plot PDF", style = "width:100%"),
          legend_box("legend_path")
        ),
        mainPanel(width = 9,
          tabsetPanel(id = "path_tabs",
            tabPanel("Plot",
              plotOutput("path_plot", height = "640px") %>% withSpinner(color = spinner_col)),
            tabPanel("Table",
              DTOutput("path_table") %>% withSpinner(color = spinner_col))
          )
        )
      )
    ),

    # ----------------------------------------------------------------
    # 8. METHODS
    # ----------------------------------------------------------------
    tabPanel("Methods", value = "tab_methods",
      fluidPage(fluidRow(column(10,
        h5("Auto-compiled Methods"),
        tags$p(tags$em("Click Refresh after running analyses. Edit before copying.")),
        actionButton("refresh_methods", "Refresh", class = "btn-primary"),
        tags$br(), tags$br(),
        textAreaInput("methods_text", NULL, rows = 18, width = "100%",
                      placeholder = "Run analyses to populate."),
        tags$hr(),
        actionButton("copy_methods", "Copy to clipboard",
                     onclick = "navigator.clipboard.writeText(document.getElementById('methods_text').value)")
      )))
    ),

    # ----------------------------------------------------------------
    # 9. EXPORT (PPTX)
    # ----------------------------------------------------------------
    tabPanel("Export", value = "tab_export",
      sidebarLayout(
        sidebarPanel(width = 3,
          h5("Figures to include"),
          checkboxInput("exp_libsize",  "Library Size",          TRUE),
          checkboxInput("exp_density",  "Expression Density",    TRUE),
          checkboxInput("exp_pca",      "PCA",                   TRUE),
          checkboxInput("exp_volcano",  "Volcano Plot",          TRUE),
          checkboxInput("exp_ma",       "MA Plot",               TRUE),
          checkboxInput("exp_omni",     "Omnibus Plot",          FALSE),
          checkboxInput("exp_heatmap",  "Heatmap",               TRUE),
          checkboxInput("exp_boxplot",  "Boxplots",              TRUE),
          checkboxInput("exp_pathway",  "Pathway Plot",          TRUE),
          checkboxInput("exp_methods",  "Include Methods slide", TRUE),
          tags$hr(),
          textInput("export_filename", "Output filename", value = "BulkSeq_Report.pptx"),
          actionButton("build_pptx", "Build PPTX", class = "btn-primary", style = "width:100%"),
          tags$br(), tags$br(),
          downloadButton("dl_pptx", "Download PPTX", style = "width:100%")
        ),
        mainPanel(width = 9,
          h5("Session figures"),
          verbatimTextOutput("export_summary")
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  rv <- reactiveValues(
    dge              = NULL,
    log2cpm          = NULL,
    meta             = NULL,
    meta_colors      = NULL,
    sample_mismatch  = list(meta_only = character(0), counts_only = character(0)),
    de               = NULL,
    de_mode          = NULL,
    pca_plot         = NULL,
    volcano_plot     = NULL,
    ma_plot          = NULL,
    omni_plot        = NULL,
    hm_plot          = NULL,
    hm_n_genes_shown = 0L,
    bp_plot          = NULL,
    path_plot        = NULL,
    path_res         = NULL,
    pptx_path        = NULL,
    params           = list()
  )

  none_vec <- c("None" = "None")

  # User-defined metadata columns only (no edgeR internals)
  user_mc <- reactive({
    req(rv$meta)
    setdiff(colnames(rv$meta), c("lib.size", "norm.factors"))
  })

  # Color map for a single column
  cmap <- function(col) {
    if (is.null(rv$meta_colors) || is.null(col) || !col %in% names(rv$meta_colors)) NULL
    else rv$meta_colors[[col]]
  }

  # Current DE result (respects pair selector for allpairs)
  current_de <- reactive({
    req(rv$de)
    if (!is.null(rv$de_mode) && rv$de_mode == "allpairs") {
      req(input$de_pair_view, input$de_pair_view %in% names(rv$de))
      rv$de[[input$de_pair_view]]
    } else rv$de
  })

  # ============================================================
  # LOAD DATA
  # ============================================================

  observeEvent(input$load_btn, {
    req(input$counts_file, input$meta_file)
    withProgress(message = "Loading & filtering...", {
      tryCatch({
        # Detect sample mismatches before filtering
        meta_raw <- read.csv(input$meta_file$datapath, as.is = TRUE, row.names = 1)
        hdr      <- read.csv(input$counts_file$datapath, check.names = FALSE, nrows = 0)
        gc_det   <- .detect_gene_col(hdr)
        cnt_cols <- setdiff(colnames(hdr), gc_det %||% character(0))
        rv$sample_mismatch <- list(
          meta_only   = setdiff(rownames(meta_raw), cnt_cols),
          counts_only = setdiff(cnt_cols, rownames(meta_raw))
        )

        dge <- load_and_filter(
          input$counts_file$datapath,
          input$meta_file$datapath,
          gene_types     = input$gene_types,
          filter_by_type = input$filter_by_type,
          min_cpm        = input$min_cpm,
          min_samples    = input$min_samples
        )
        rv$dge         <- dge
        rv$log2cpm     <- make_log2cpm(dge)
        rv$meta        <- dge$samples
        rv$meta_colors <- build_meta_colors(rv$meta)
        rv$de          <- NULL

        rv$params$min_cpm        <- input$min_cpm
        rv$params$min_samples    <- input$min_samples
        rv$params$filter_by_type <- input$filter_by_type
        rv$params$gene_types     <- input$gene_types

        mc   <- user_mc()
        grps <- sort(unique(as.character(dge$samples[[mc[1]]])))

        for (id in c("qc_color","pca_color","de_col","bp_x","bp_color"))
          updateSelectInput(session, id, choices = mc, selected = mc[1])
        updateSelectInput(session, "pca_shape",
                          choices = c(none_vec, mc), selected = "None")
        updateSelectInput(session, "bp_facet",
                          choices = c(none_vec, mc), selected = "None")
        updateCheckboxGroupInput(session, "hm_ann_cols",
                                 choices = mc, selected = mc[1:min(2,length(mc))])
        updateCheckboxGroupInput(session, "de_groups_sel",
                                 choices = grps, selected = grps)
        updateSelectInput(session, "bp_ref_group", choices = grps, selected = grps[1])
        # GSVA group selectors
        updateSelectInput(session, "gsva_de_col", choices = mc, selected = mc[1])
        updateSelectInput(session, "gsva_de_a", choices = grps,
                          selected = grps[min(1,length(grps))])
        updateSelectInput(session, "gsva_de_b", choices = grps,
                          selected = grps[min(2,length(grps))])

        showNotification(
          paste0("Loaded: ", nrow(dge), " genes × ", ncol(dge), " samples"),
          type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Load error:", conditionMessage(e)),
                         type = "error", duration = 10)
      })
    })
  })

  # Sample mismatch warning banner
  output$sample_warning <- renderUI({
    mm <- rv$sample_mismatch
    msgs <- c(
      if (length(mm$meta_only) > 0)
        paste0(length(mm$meta_only), " sample(s) in metadata but NOT in counts: ",
               paste(mm$meta_only, collapse = ", ")),
      if (length(mm$counts_only) > 0)
        paste0(length(mm$counts_only), " sample(s) in counts but NOT in metadata: ",
               paste(mm$counts_only, collapse = ", "))
    )
    if (!length(msgs)) return(NULL)
    div(class = "warn-box",
        tags$b("⚠ Sample mismatch — these samples were excluded:"),
        tags$ul(lapply(msgs, tags$li)))
  })

  # DE group dropdowns
  observeEvent(input$de_col, {
    req(rv$meta, input$de_col %in% colnames(rv$meta))
    grps <- sort(unique(as.character(rv$meta[[input$de_col]])))
    updateSelectInput(session, "de_a", choices = grps, selected = grps[min(1,length(grps))])
    updateSelectInput(session, "de_b", choices = grps, selected = grps[min(2,length(grps))])
    updateCheckboxGroupInput(session, "de_groups_sel", choices = grps, selected = grps)
    updateSelectInput(session, "bp_ref_group", choices = grps, selected = grps[1])
  })

  # GSVA group dropdowns
  observeEvent(input$gsva_de_col, {
    req(rv$meta, input$gsva_de_col %in% colnames(rv$meta))
    grps <- sort(unique(as.character(rv$meta[[input$gsva_de_col]])))
    updateSelectInput(session, "gsva_de_a", choices = grps, selected = grps[min(1,length(grps))])
    updateSelectInput(session, "gsva_de_b", choices = grps, selected = grps[min(2,length(grps))])
  })

  # ---- Previews ----
  output$counts_preview <- renderDT({
    req(input$counts_file)
    .dt_safe(head(read.csv(input$counts_file$datapath, check.names = FALSE), 5))
  }, options = list(scrollX = TRUE, dom = "t"), rownames = FALSE)

  output$meta_preview <- renderDT({
    req(input$meta_file)
    d <- tryCatch(
      read.csv(input$meta_file$datapath, as.is = TRUE, check.names = FALSE),
      error = function(e) data.frame(Error = conditionMessage(e))
    )
    .dt_safe(d)
  }, options = list(scrollX = TRUE, dom = "t", pageLength = -1))

  output$filter_summary <- renderText({
    req(rv$dge)
    paste0("Retained: ", nrow(rv$dge), " genes\n",
           "Samples:  ", ncol(rv$dge), "\n",
           "Norm:     TMM")
  })

  # ============================================================
  # SAMPLE RELABELING
  # ============================================================

  # Dynamic label-option inputs (rendered from rv$meta so choices always current)
  output$relabel_opts_ui <- renderUI({
    req(rv$meta)
    mc <- user_mc()
    switch(input$relabel_mode %||% "none",
      col    = selectInput("relabel_col", "Column", choices = mc, selected = mc[1]),
      concat = tagList(
        selectInput("relabel_col1", "Column 1", choices = mc, selected = mc[1]),
        selectInput("relabel_col2", "Column 2", choices = mc,
                    selected = mc[min(2, length(mc))]),
        textInput("relabel_sep", "Separator", value = "_")
      ),
      NULL
    )
  })

  output$relabel_sort_ui <- renderUI({
    req(rv$meta)
    mc <- user_mc()
    selectInput("relabel_sortby", "Sort samples by",
                choices = c("None" = "none", mc), selected = "none")
  })

  output$relabel_preview <- renderText({
    req(rv$meta)
    old <- rownames(rv$meta)
    new <- switch(input$relabel_mode %||% "none",
      none   = return("Labels unchanged."),
      col    = { col <- input$relabel_col; req(col, col %in% colnames(rv$meta))
                 as.character(rv$meta[[col]]) },
      concat = { c1 <- input$relabel_col1; c2 <- input$relabel_col2
                 sep <- input$relabel_sep %||% "_"
                 req(c1, c2, c1 %in% colnames(rv$meta), c2 %in% colnames(rv$meta))
                 paste0(rv$meta[[c1]], sep, rv$meta[[c2]]) }
    )
    n <- min(length(old), 8)
    paste0("Preview (first ", n, "):\n",
           paste0("  ", old[1:n], "  →  ", new[1:n], collapse = "\n"))
  })

  observeEvent(input$apply_relabel, {
    req(rv$log2cpm, rv$meta)
    old <- rownames(rv$meta)

    new <- switch(input$relabel_mode %||% "none",
      none   = old,
      col    = { col <- input$relabel_col; req(col %in% colnames(rv$meta))
                 as.character(rv$meta[[col]]) },
      concat = { c1 <- input$relabel_col1; c2 <- input$relabel_col2
                 sep <- input$relabel_sep %||% "_"
                 req(c1 %in% colnames(rv$meta), c2 %in% colnames(rv$meta))
                 paste0(rv$meta[[c1]], sep, rv$meta[[c2]]) }
    )

    # Make unique if duplicated
    if (anyDuplicated(new)) {
      showNotification("Duplicate labels — appending suffix.", type = "warning")
      cntr <- setNames(integer(length(unique(new))), unique(new))
      for (i in seq_along(new)) {
        cntr[[new[i]]] <- cntr[[new[i]]] + 1L
        if (cntr[[new[i]]] > 1) new[i] <- paste0(new[i], "_", cntr[[new[i]]])
      }
    }

    # Sort order
    sortby <- input$relabel_sortby %||% "none"
    ord <- if (sortby != "none" && sortby %in% colnames(rv$meta))
      order(rv$meta[[sortby]], decreasing = input$relabel_sortdir == "desc")
    else seq_len(nrow(rv$meta))

    rv$meta                   <- rv$meta[ord, , drop = FALSE]
    rownames(rv$meta)         <- new[ord]
    rv$log2cpm                <- rv$log2cpm[, ord, drop = FALSE]
    colnames(rv$log2cpm)      <- new[ord]
    rv$dge$samples            <- rv$dge$samples[ord, , drop = FALSE]
    rownames(rv$dge$samples)  <- new[ord]
    # Keep counts column names in sync so cpm() → plot_density merge works
    rv$dge$counts             <- rv$dge$counts[, ord, drop = FALSE]
    colnames(rv$dge$counts)   <- new[ord]

    rv$meta_colors <- rebuild_meta_colors(rv$meta, rv$meta_colors)
    showNotification(paste0("Applied to ", length(ord), " samples."),
                     type = "message", duration = 3)
  })

  observeEvent(input$reset_relabel, {
    req(rv$dge)
    rv$log2cpm   <- make_log2cpm(rv$dge)
    rv$meta      <- rv$dge$samples
    rv$meta_colors <- rebuild_meta_colors(rv$meta, rv$meta_colors)
    showNotification("Labels reset to original IDs.", type = "message", duration = 3)
  })

  # ============================================================
  # COLOR EDITOR
  # ============================================================

  output$has_colors <- reactive({ !is.null(rv$meta_colors) && length(rv$meta_colors) > 0 })
  outputOptions(output, "has_colors", suspendWhenHidden = FALSE)

  output$color_editor_ui <- renderUI({
    req(rv$meta_colors)
    has_picker <- requireNamespace("colourpicker", quietly = TRUE)
    if (!has_picker)
      return(tags$p(tags$em(
        "Install colourpicker for an interactive color editor: ",
        tags$code("install.packages('colourpicker')"),
        " then restart the app.")))

    col_blocks <- lapply(names(rv$meta_colors), function(col) {
      pickers <- lapply(names(rv$meta_colors[[col]]), function(grp) {
        inp_id <- paste0("clr_", make.names(col), "_X_", make.names(grp))
        column(3,
          colourpicker::colourInput(
            inputId    = inp_id,
            label      = grp,
            value      = rv$meta_colors[[col]][[grp]],
            showColour = "both",
            palette    = "square"
          )
        )
      })
      tagList(tags$b(col, ":"), fluidRow(pickers), tags$br())
    })
    tagList(col_blocks)
  })

  observeEvent(input$apply_colors, {
    req(rv$meta_colors)
    if (!requireNamespace("colourpicker", quietly = TRUE)) return()
    new_cols <- rv$meta_colors
    for (col in names(new_cols)) {
      for (grp in names(new_cols[[col]])) {
        inp_id <- paste0("clr_", make.names(col), "_X_", make.names(grp))
        val <- input[[inp_id]]
        if (!is.null(val) && nzchar(trimws(val)))
          new_cols[[col]][[grp]] <- val
      }
    }
    rv$meta_colors <- new_cols
    showNotification("Colors updated.", type = "message", duration = 2)
  })

  # ============================================================
  # QC
  # ============================================================

  qc_h     <- reactive(input$qc_height)
  col_by_qc <- reactive({
    req(rv$meta)
    if (!is.null(input$qc_color) && input$qc_color %in% colnames(rv$meta))
      input$qc_color else NULL
  })

  output$lib_size_plot <- renderPlot({
    req(rv$dge)
    plot_library_size(rv$dge, color_by = col_by_qc(), color_map = cmap(col_by_qc()),
                      sort_by_size = isTRUE(input$qc_sort_libsize))
  }, height = qc_h)

  output$density_plot <- renderPlot({
    req(rv$dge)
    plot_density(rv$dge, color_by = col_by_qc(), color_map = cmap(col_by_qc()))
  }, height = qc_h)

  observeEvent(rv$meta, {
    updateTextAreaInput(session, "legend_libsize", value = legend_libsize(col_by_qc()))
    updateTextAreaInput(session, "legend_density", value = legend_density(col_by_qc()))
  })

  output$dl_libsize <- downloadHandler(
    filename = function() "library_size.pdf",
    content  = function(f) {
      req(rv$dge)
      ggplot2::ggsave(f, plot_library_size(rv$dge, color_by = col_by_qc(),
                                            color_map = cmap(col_by_qc()),
                                            sort_by_size = isTRUE(input$qc_sort_libsize)),
                      width = 8, height = 5, device = "pdf")
    }
  )
  output$dl_density <- downloadHandler(
    filename = function() "expression_density.pdf",
    content  = function(f) {
      req(rv$dge)
      ggplot2::ggsave(f, plot_density(rv$dge, color_by = col_by_qc(),
                                       color_map = cmap(col_by_qc())),
                      width = 8, height = 5, device = "pdf")
    }
  )

  # ============================================================
  # PCA
  # ============================================================

  observeEvent(input$run_pca, {
    req(rv$log2cpm, rv$meta)
    withProgress(message = "Running PCA...", {
      shape_by <- if (!is.null(input$pca_shape) && input$pca_shape != "None")
        input$pca_shape else NULL
      rv$pca_plot <- plot_pca(
        rv$log2cpm, rv$meta,
        color_by      = input$pca_color,
        shape_by      = shape_by,
        n_genes       = input$pca_ngenes,
        label_samples = input$pca_label,
        color_map     = cmap(input$pca_color)
      )
      rv$params$pca_ngenes <- input$pca_ngenes
      updateTextAreaInput(session, "legend_pca",
                          value = legend_pca(input$pca_color, shape_by, input$pca_ngenes))
    })
  })

  output$pca_plot <- renderPlot({ req(rv$pca_plot); rv$pca_plot }, height = 520)

  output$dl_pca <- downloadHandler(
    filename = function() "pca.pdf",
    content  = function(f) {
      req(rv$pca_plot)
      ggplot2::ggsave(f, rv$pca_plot, width = 7, height = 6, device = "pdf")
    }
  )

  # ============================================================
  # DE ANALYSIS
  # ============================================================

  observeEvent(input$run_de, {
    req(rv$dge, input$de_col)
    withProgress(message = "Running edgeR glmQLFTest...", {
      tryCatch({
        mode <- input$de_mode; rv$de_mode <- mode
        if (mode == "2group") {
          req(input$de_a, input$de_b)
          rv$de <- run_DE(rv$dge, input$de_col, input$de_a, input$de_b,
                          fdr = input$de_fdr, lfc = input$de_lfc)
          rv$params$de_a <- input$de_a; rv$params$de_b <- input$de_b
        } else if (mode == "allpairs") {
          grps <- if (length(input$de_groups_sel) >= 2) input$de_groups_sel else NULL
          rv$de <- run_DE_allpairs(rv$dge, input$de_col, groups = grps,
                                    fdr = input$de_fdr, lfc = input$de_lfc)
          updateSelectInput(session, "de_pair_view",
                            choices = names(rv$de), selected = names(rv$de)[1])
        } else {
          grps <- if (length(input$de_groups_sel) >= 2) input$de_groups_sel else NULL
          rv$de <- run_DE_omnibus(rv$dge, input$de_col, groups = grps, fdr = input$de_fdr)
        }
        rv$params$de_mode <- mode
        rv$params$de_fdr  <- input$de_fdr
        rv$params$de_lfc  <- input$de_lfc
        rv$params$de_col  <- input$de_col
        showNotification("DE complete.", type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("DE error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })
  })

  output$de_summary <- renderText({
    req(rv$de, rv$de_mode)
    if (rv$de_mode == "2group")
      paste0("Up:    ", nrow(rv$de$up), "\nDown:  ", nrow(rv$de$down),
             "\nTotal: ", nrow(rv$de$sig))
    else if (rv$de_mode == "allpairs") paste0(length(rv$de), " pairs computed")
    else paste0("Omnibus sig: ", nrow(rv$de$sig), " genes")
  })

  output$de_table <- renderDT({
    req(rv$de, rv$de_mode)
    d <- if (rv$de_mode == "2group") rv$de$all[order(rv$de$all$FDR), ]
    else if (rv$de_mode == "allpairs") {
      req(input$de_pair_view, input$de_pair_view %in% names(rv$de))
      rv$de[[input$de_pair_view]]$all[order(rv$de[[input$de_pair_view]]$all$FDR), ]
    } else rv$de$all[order(rv$de$all$FDR), ]
    .dt_safe(d[, setdiff(colnames(d), "genes")])
  }, options = list(pageLength = 20, scrollX = TRUE), filter = "top")

  output$volcano_plot <- renderPlot({
    req(rv$de_mode != "omnibus"); d <- current_de(); req(d)
    p <- plot_volcano(d, fdr = input$de_fdr, lfc = input$de_lfc)
    rv$volcano_plot <- p; p
  }, height = 520)

  output$ma_plot <- renderPlot({
    req(rv$de_mode != "omnibus"); d <- current_de(); req(d)
    p <- plot_ma(d, fdr = input$de_fdr, lfc = input$de_lfc)
    rv$ma_plot <- p; p
  }, height = 520)

  output$omni_plot <- renderPlot({
    req(rv$de_mode == "omnibus", rv$de)
    p <- plot_omni(rv$de, fdr = input$de_fdr)
    rv$omni_plot <- p; p
  }, height = 520)

  observe({
    req(rv$de, rv$de_mode)
    if (rv$de_mode == "2group") {
      updateTextAreaInput(session, "legend_volcano",
        value = legend_volcano(input$de_a, input$de_b, input$de_fdr, input$de_lfc))
      updateTextAreaInput(session, "legend_ma",
        value = legend_ma(input$de_a, input$de_b, input$de_fdr, input$de_lfc))
    } else if (rv$de_mode == "allpairs") {
      parts <- strsplit(input$de_pair_view %||% "", " vs ")[[1]]
      if (length(parts) == 2) {
        updateTextAreaInput(session, "legend_volcano",
          value = legend_volcano(parts[1], parts[2], input$de_fdr, input$de_lfc))
        updateTextAreaInput(session, "legend_ma",
          value = legend_ma(parts[1], parts[2], input$de_fdr, input$de_lfc))
      }
    }
  })

  output$dl_de_table <- downloadHandler(
    filename = function() paste0("DE_", input$de_mode, ".csv"),
    content  = function(f) {
      req(rv$de)
      d <- if (rv$de_mode == "2group") rv$de$all
      else if (rv$de_mode == "allpairs") {
        dfs <- lapply(names(rv$de), function(nm) cbind(pair = nm, rv$de[[nm]]$all))
        do.call(rbind, dfs)
      } else rv$de$all
      write.csv(d[order(d$FDR), ], f)
    }
  )
  output$dl_volcano <- downloadHandler("volcano.pdf", function(f) {
    req(rv$volcano_plot); ggplot2::ggsave(f, rv$volcano_plot, width=7, height=6, device="pdf") })
  output$dl_ma <- downloadHandler("ma_plot.pdf", function(f) {
    req(rv$ma_plot); ggplot2::ggsave(f, rv$ma_plot, width=7, height=6, device="pdf") })
  output$dl_omni <- downloadHandler("omnibus_plot.pdf", function(f) {
    req(rv$omni_plot); ggplot2::ggsave(f, rv$omni_plot, width=7, height=6, device="pdf") })

  # ============================================================
  # HEATMAP
  # ============================================================

  output$hm_pairs_ui <- renderUI({
    req(rv$de_mode == "allpairs", rv$de)
    checkboxGroupInput("hm_de_pairs", "Comparisons to include",
                       choices = names(rv$de), selected = names(rv$de))
  })

  # Pathway tab: pair selector shown only in allpairs DE mode (for GSEA/ORA)
  output$path_pair_ui <- renderUI({
    req(rv$de_mode == "allpairs", rv$de)
    selectInput("path_de_pair_sel", "DE comparison (for GSEA/ORA)",
                choices = names(rv$de), selected = names(rv$de)[1])
  })

  observeEvent(input$run_hm, {
    req(rv$log2cpm, rv$meta)
    genes <- if (input$hm_source == "custom") {
      parse_genes(input$hm_genes)
    } else {
      req(rv$de, rv$de_mode)
      if (rv$de_mode == "allpairs") {
        pairs <- if (!is.null(input$hm_de_pairs) && length(input$hm_de_pairs) > 0)
          input$hm_de_pairs else names(rv$de)
        gene_sets <- Filter(length, lapply(pairs, function(p) {
          d   <- rv$de[[p]]
          src <- if (isTRUE(input$hm_de_sig) && nrow(d$sig) > 0) d$sig else d$all
          if (nrow(src) == 0) return(character(0))
          rownames(head(src[order(src$FDR), ], input$hm_n_genes))
        }))
        if (!length(gene_sets)) { showNotification("No DE genes found.", type="warning"); return() }
        if (input$hm_de_gene_rule == "intersection" && length(gene_sets) > 1)
          Reduce(intersect, gene_sets) else unique(unlist(gene_sets))
      } else {
        src <- if (isTRUE(input$hm_de_sig)) rv$de$sig else rv$de$all
        if (nrow(src) == 0) { showNotification("No significant genes.", type="warning"); return() }
        rownames(head(src[order(src$FDR), ], input$hm_n_genes))
      }
    }
    if (length(genes) < 2) { showNotification("Need ≥ 2 genes.", type="warning"); return() }

    ann_cols   <- if (length(input$hm_ann_cols) > 0) input$hm_ann_cols else NULL
    ann_colors <- if (!is.null(ann_cols) && !is.null(rv$meta_colors))
      rv$meta_colors[intersect(ann_cols, names(rv$meta_colors))] else NULL

    withProgress(message = "Building heatmap...", {
      tryCatch({
        n_g <- length(intersect(genes, rownames(rv$log2cpm)))
        rv$hm_n_genes_shown <- n_g
        rv$hm_plot <- plot_heatmap(
          rv$log2cpm, genes, rv$meta,
          ann_cols     = ann_cols,
          cluster_rows = input$hm_cluster_rows,
          cluster_cols = input$hm_cluster_cols,
          scale        = input$hm_scale,
          ann_colors   = ann_colors
        )
        updateTextAreaInput(session, "legend_heatmap",
                            value = legend_heatmap(n_g, ann_cols, input$hm_scale))
      }, error = function(e) {
        showNotification(paste("Heatmap error:", conditionMessage(e)),
                         type="error", duration=8)
      })
    })
  })

  hm_height <- reactive({
    n <- rv$hm_n_genes_shown
    if (n == 0L) return(640L)
    ch <- max(6L, min(14L, floor(400L / n)))
    as.integer(min(max(640L, n * ch + 260L), 4000L))
  })

  output$heatmap_plot <- renderPlot({
    req(rv$hm_plot)
    grid::grid.newpage()
    grid::grid.draw(rv$hm_plot$gtable)
  }, height = hm_height)

  output$dl_heatmap <- downloadHandler("heatmap.pdf", function(f) {
    req(rv$hm_plot); pdf(f, width=10, height=8)
    grid::grid.newpage(); grid::grid.draw(rv$hm_plot$gtable); dev.off()
  })

  # ============================================================
  # BOXPLOTS
  # ============================================================

  output$bp_norm_result <- renderText({
    req(rv$log2cpm, rv$meta, input$bp_x)
    genes <- parse_genes(input$bp_genes)
    if (!length(genes)) return("Enter at least one gene.")
    res <- check_normality(rv$log2cpm, genes[1], rv$meta, input$bp_x)
    if (is.null(res)) return("Gene not found.")
    n_groups <- length(unique(as.character(rv$meta[[input$bp_x]])))
    rec      <- recommend_test(res, n_groups)
    lines <- sapply(names(res), function(g) {
      r  <- res[[g]]
      nm <- if (is.na(r$normal)) "—" else if (isTRUE(r$normal)) "normal" else "non-normal"
      sprintf("  %s: n=%d, W=%.3f, p=%.4f (%s)", g, r$n, r$W %||% NA, r$p %||% NA, nm)
    })
    paste0("Shapiro-Wilk (", genes[1], "):\n",
           paste(lines, collapse="\n"), "\n\nRecommended: ", rec)
  }) %>% bindEvent(input$bp_norm_check)

  observeEvent(input$bp_norm_check, {
    req(rv$log2cpm, rv$meta, input$bp_x)
    genes <- parse_genes(input$bp_genes); if (!length(genes)) return()
    res <- check_normality(rv$log2cpm, genes[1], rv$meta, input$bp_x)
    if (is.null(res)) return()
    n_groups <- length(unique(as.character(rv$meta[[input$bp_x]])))
    updateSelectInput(session, "bp_test", selected = recommend_test(res, n_groups))
  })

  observeEvent(input$run_bp, {
    req(rv$log2cpm, rv$meta, nchar(trimws(input$bp_genes)) > 0, input$bp_x)
    genes <- parse_genes(input$bp_genes)
    facet <- if (!is.null(input$bp_facet) && input$bp_facet != "None") input$bp_facet else NULL
    ref_g <- if (input$bp_comp_style == "reference") input$bp_ref_group else NULL
    withProgress(message = "Building boxplots...", {
      tryCatch({
        rv$bp_plot <- plot_boxplots(
          rv$log2cpm, rv$meta,
          genes      = genes, x_by = input$bp_x,
          color_by   = input$bp_color, facet_by = facet,
          test       = input$bp_test,  posthoc  = input$bp_posthoc,
          comp_style = input$bp_comp_style, ref_group = ref_g,
          p_display  = input$bp_p_display,
          color_map  = cmap(input$bp_color)
        )
        rv$params$bp_test    <- input$bp_test
        rv$params$bp_posthoc <- input$bp_posthoc
        updateTextAreaInput(session, "legend_boxplot",
          value = legend_boxplot(genes, input$bp_x, input$bp_test,
                                  input$bp_posthoc, input$bp_comp_style))
      }, error = function(e) {
        showNotification(paste("Boxplot error:", conditionMessage(e)),
                         type="error", duration=8)
      })
    })
  })

  output$boxplot_plot <- renderPlot({ req(rv$bp_plot); rv$bp_plot }, height = 520)

  output$dl_boxplot <- downloadHandler("boxplots.pdf", function(f) {
    req(rv$bp_plot)
    ggplot2::ggsave(f, rv$bp_plot, width=input$bp_width, height=input$bp_height_in, device="pdf")
  })

  # ============================================================
  # PATHWAYS
  # ============================================================

  observeEvent(input$run_path, {
    req(rv$de)
    withProgress(message = paste("Running", input$path_method, "..."), {
      tryCatch({
        log2cpm_use <- rv$log2cpm
        de_src <- if (!is.null(rv$de_mode) && rv$de_mode == "allpairs") {
          pair <- input$path_de_pair_sel %||% names(rv$de)[1]
          req(pair %in% names(rv$de))
          rv$de[[pair]]
        } else rv$de

        # Mouse → human conversion
        if (input$path_species == "mouse2human") {
          setProgress(0.15, detail = "Converting via biomaRt...")
          conv <- convert_mouse_human(rownames(log2cpm_use))
          if (!is.null(conv) && length(conv) > 0) {
            log2cpm_use <- apply_human_symbols(log2cpm_use, conv)
            common <- intersect(rownames(de_src$all), names(conv))
            de_src$all  <- de_src$all[common, ]; rownames(de_src$all) <- conv[common]
            de_src$sig  <- de_src$all[!is.na(de_src$all$FDR) &
              de_src$all$FDR < input$de_fdr & abs(de_src$all$logFC) > input$de_lfc, ]
            de_src$up   <- de_src$all[!is.na(de_src$all$FDR) &
              de_src$all$FDR < input$de_fdr & de_src$all$logFC > input$de_lfc, ]
            de_src$down <- de_src$all[!is.na(de_src$all$FDR) &
              de_src$all$FDR < input$de_fdr & de_src$all$logFC < -input$de_lfc, ]
          } else showNotification("Conversion failed — using original symbols.", type="warning")
        }

        # Gene sets
        gene_sets <- if (input$path_geneset_src == "msigdb") {
          setProgress(0.2, detail = "Loading MSigDB...")
          load_msigdb_sets(input$path_species_msig, input$path_collection)
        } else {
          req(input$gmt_file); fgsea::gmtPathways(input$gmt_file$datapath)
        }

        species_use <- if (input$path_species %in% c("mouse","mouse2human")) "mouse" else "human"

        if (input$path_method == "gsea") {
          setProgress(0.4, detail = "Running fgsea...")
          n_top_gsea <- max(3L, as.integer(input$gsea_n_top %||% 15))
          res          <- run_gsea(de_src$all, gene_sets, pval = input$path_pval,
                                    n_top = n_top_gsea)
          rv$path_res  <- res$results
          rv$path_plot <- res$plot
          if (nrow(res$results) == 0 || is.null(res$plot)) {
            n_ov  <- res$n_overlap %||% "?"
            n_set <- res$n_sets    %||% "?"
            showNotification(
              paste0("GSEA returned no results. ",
                     n_ov, " of your genes matched ", n_set, " gene sets. ",
                     "Verify that species/symbols match the database."),
              type = "warning", duration = 12)
          } else {
            n_sig <- sum(!is.na(res$results$padj) & res$results$padj < input$path_pval, na.rm = TRUE)
            showNotification(
              paste0("GSEA complete: ", nrow(res$results), " pathways tested, ",
                     n_sig, " significant at padj < ", input$path_pval,
                     ". Showing top ", input$gsea_n_top %||% 15, " per direction."),
              type = "message", duration = 6)
          }

        } else if (input$path_method == "gsva") {
          req(input$gsva_de_col, input$gsva_de_a, input$gsva_de_b)
          setProgress(0.4, detail = "Running GSVA...")
          res <- run_gsva(log2cpm_use, gene_sets, rv$meta,
                          input$gsva_de_col, input$gsva_de_a, input$gsva_de_b,
                          fdr = input$path_pval)
          rv$path_res <- res$de
          if (nrow(res$de) > 0) {
            n_top_gsva <- max(3L, as.integer(input$gsva_n_top %||% 25))
            # Show top N up + top N down by logFC (capped to available)
            de_ord  <- res$de[order(res$de$logFC, decreasing = TRUE), , drop = FALSE]
            top_up  <- head(rownames(de_ord[de_ord$logFC > 0, ]), n_top_gsva)
            top_dn  <- tail(rownames(de_ord[de_ord$logFC < 0, ]), n_top_gsva)
            show_pw <- unique(c(top_up, rev(top_dn)))
            show_pw <- intersect(show_pw, rownames(res$scores))
            if (length(show_pw) == 0) show_pw <- rownames(res$scores)

            ann   <- rv$meta[, input$gsva_de_col, drop = FALSE]
            ann_c <- if (!is.null(rv$meta_colors))
              rv$meta_colors[intersect(input$gsva_de_col, names(rv$meta_colors))] else NULL
            rv$path_plot <- pheatmap::pheatmap(
              res$scores[show_pw, , drop = FALSE],
              cluster_cols = FALSE, scale = "row", annotation_col = ann,
              annotation_colors = ann_c, border_color = "black",
              color = colorRampPalette(c(unname(Nour18["darkblue"]), "white",
                                         unname(Nour18["lightred"])))(200),
              silent = TRUE)
          } else {
            rv$path_plot <- NULL
            showNotification("GSVA: no significant gene sets at this threshold.",
                             type = "warning", duration = 8)
          }

        } else if (input$path_method == "ora") {
          setProgress(0.4, detail = "Running ORA...")
          rv$path_res <- run_ora(de_src, species = species_use,
                                  ont = input$ora_ont, pval = input$path_pval)
          if (any(c("up_go","dn_go") %in% names(rv$path_res)) &&
              requireNamespace("enrichplot", quietly = TRUE)) {
            up_go <- rv$path_res[["up_go"]]; dn_go <- rv$path_res[["dn_go"]]
            merged <- if (!is.null(up_go) && !is.null(dn_go))
              clusterProfiler::merge_result(list(Up = up_go, Down = dn_go))
            else up_go %||% dn_go
            rv$path_plot <- if (!is.null(merged))
              enrichplot::dotplot(merged, showCategory = 15) else NULL
          } else rv$path_plot <- NULL
        }

        rv$params$path_method     <- input$path_method
        rv$params$path_use_msigdb <- (input$path_geneset_src == "msigdb")
        rv$params$path_collection <- input$path_collection
        rv$params$path_species    <- input$path_species_msig
        rv$params$path_pval       <- input$path_pval
        rv$params$ora_ont         <- input$ora_ont

        col_val <- if (input$path_geneset_src == "msigdb") input$path_collection else NULL
        updateTextAreaInput(session, "legend_path",
          value = if (input$path_method == "gsea")
                    legend_gsea(col_val, input$path_species_msig, input$path_pval)
                  else if (input$path_method == "gsva")
                    legend_gsva(col_val, input$path_species_msig, input$path_pval)
                  else "")

        showNotification("Pathway analysis complete.", type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Pathway error:", conditionMessage(e)),
                         type = "error", duration = 12)
      })
    })
  })

  path_plot_height <- reactive({
    if (!is.null(input$path_method) && input$path_method == "gsva")
      max(200L, as.integer(input$gsva_hm_height %||% 640L))
    else 640L
  })

  output$path_plot <- renderPlot({
    req(rv$path_plot)
    if (inherits(rv$path_plot, "pheatmap")) {
      grid::grid.newpage(); grid::grid.draw(rv$path_plot$gtable)
    } else print(rv$path_plot)
  }, height = path_plot_height)

  output$path_table <- renderDT({
    req(rv$path_res)
    df <- if (is.data.frame(rv$path_res)) {
      rv$path_res
    } else if (is.list(rv$path_res)) {
      dfs <- lapply(names(rv$path_res), function(nm) {
        d <- tryCatch(as.data.frame(rv$path_res[[nm]]), error = function(e) NULL)
        if (!is.null(d) && nrow(d) > 0) { d$set <- nm; d } else NULL
      })
      dplyr::bind_rows(Filter(Negate(is.null), dfs))
    } else data.frame(message = "No results.")
    if (nrow(df) == 0) df <- data.frame(message = "No significant results at this threshold.")
    .dt_safe(df)
  }, options = list(pageLength = 20, scrollX = TRUE))

  output$dl_path_results <- downloadHandler(
    filename = function() paste0(input$path_method, "_results.csv"),
    content  = function(f) {
      req(rv$path_res)
      df <- if (is.data.frame(rv$path_res)) rv$path_res else {
        dfs <- lapply(names(rv$path_res), function(nm) {
          d <- tryCatch(as.data.frame(rv$path_res[[nm]]), error = function(e) NULL)
          if (!is.null(d) && nrow(d) > 0) { d$set <- nm; d } else NULL
        })
        dplyr::bind_rows(Filter(Negate(is.null), dfs))
      }
      write.csv(.dt_safe(df), f, row.names = FALSE)
    }
  )
  output$dl_path_plot <- downloadHandler(
    filename = function() paste0(input$path_method, "_plot.pdf"),
    content  = function(f) {
      req(rv$path_plot)
      if (inherits(rv$path_plot, "pheatmap")) {
        w <- as.numeric(input$gsva_dl_w %||% 10)
        h <- as.numeric(input$gsva_dl_h %||% 8)
        pdf(f, width = w, height = h)
        grid::grid.newpage(); grid::grid.draw(rv$path_plot$gtable); dev.off()
      } else ggplot2::ggsave(f, rv$path_plot, width = 9, height = 7, device = "pdf")
    }
  )

  # ============================================================
  # METHODS
  # ============================================================

  observeEvent(input$refresh_methods, {
    updateTextAreaInput(session, "methods_text", value = compile_methods(rv$params))
  })

  # ============================================================
  # EXPORT (PPTX)
  # ============================================================

  output$export_summary <- renderText({
    have <- c(
      if (!is.null(rv$dge))          "Library Size, Expression Density" else NULL,
      if (!is.null(rv$pca_plot))     "PCA"          else NULL,
      if (!is.null(rv$volcano_plot)) "Volcano"      else NULL,
      if (!is.null(rv$ma_plot))      "MA Plot"      else NULL,
      if (!is.null(rv$omni_plot))    "Omnibus Plot" else NULL,
      if (!is.null(rv$hm_plot))      "Heatmap"      else NULL,
      if (!is.null(rv$bp_plot))      "Boxplots"     else NULL,
      if (!is.null(rv$path_plot))    "Pathway Plot" else NULL
    )
    if (!length(have)) "No figures generated yet."
    else paste0("Available:\n", paste0("  • ", have, collapse = "\n"))
  })

  observeEvent(input$build_pptx, {
    withProgress(message = "Building PPTX...", {
      col_by <- col_by_qc()
      figs <- list()
      if (input$exp_libsize && !is.null(rv$dge))
        figs <- c(figs, list(
          list(name="Library Size",
               plot=plot_library_size(rv$dge, color_by=col_by, color_map=cmap(col_by),
                                      sort_by_size=isTRUE(input$qc_sort_libsize)),
               legend=input$legend_libsize),
          list(name="Expression Density",
               plot=plot_density(rv$dge,color_by=col_by,color_map=cmap(col_by)),
               legend=input$legend_density)))
      if (input$exp_pca     && !is.null(rv$pca_plot))
        figs <- c(figs, list(list(name="PCA",         plot=rv$pca_plot,     legend=input$legend_pca)))
      if (input$exp_volcano && !is.null(rv$volcano_plot))
        figs <- c(figs, list(list(name="Volcano",     plot=rv$volcano_plot, legend=input$legend_volcano)))
      if (input$exp_ma      && !is.null(rv$ma_plot))
        figs <- c(figs, list(list(name="MA Plot",     plot=rv$ma_plot,      legend=input$legend_ma)))
      if (input$exp_omni    && !is.null(rv$omni_plot))
        figs <- c(figs, list(list(name="Omnibus",     plot=rv$omni_plot,    legend="")))
      if (input$exp_heatmap && !is.null(rv$hm_plot))
        figs <- c(figs, list(list(name="Heatmap",     plot=rv$hm_plot,      legend=input$legend_heatmap)))
      if (input$exp_boxplot && !is.null(rv$bp_plot))
        figs <- c(figs, list(list(name="Boxplots",    plot=rv$bp_plot,      legend=input$legend_boxplot)))
      if (input$exp_pathway && !is.null(rv$path_plot))
        figs <- c(figs, list(list(name="Pathway",     plot=rv$path_plot,    legend=input$legend_path)))

      methods_txt <- if (input$exp_methods) input$methods_text else NULL
      out <- file.path(tempdir(), input$export_filename)
      tryCatch({
        build_bulk_pptx(figs, methods_txt, rv$params, out)
        rv$pptx_path <- out
        showNotification("PPTX ready — click Download.", type="message", duration=6)
      }, error = function(e) {
        showNotification(paste("PPTX error:", conditionMessage(e)), type="error", duration=12)
      })
    })
  })

  output$dl_pptx <- downloadHandler(
    filename = function() input$export_filename,
    content  = function(f) {
      req(rv$pptx_path, file.exists(rv$pptx_path))
      file.copy(rv$pptx_path, f)
    }
  )
}

# ============================================================
shinyApp(ui, server)
