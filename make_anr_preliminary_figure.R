#!/usr/bin/env Rscript

# ==============================================================================
# ANR preliminary-data figure: DMD snRNA-seq + TF activity
#
# Panel A  : real snRNA-seq UMAP
# Panel B  : top 10 candidate TF activity heatmap (donor-aware pseudobulk)
# Panel C  : per-nucleus JUN activity distribution across DMD cell populations
#
# Statistical interpretation:
# - The heatmap is the inferential result: donor-aware pseudobulk DESeq2
#   followed by CollecTRI + decoupleR ULM.
# - The violin plot is an exploratory per-nucleus visualization of JUN activity.
#   It is not used for donor-level hypothesis testing.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------------------------

INSTALL_MISSING <- TRUE

cran_pkgs <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "ggplot2",
  "patchwork",
  "tibble"
)

bioc_pkgs <- c(
  "decoupleR"
)

missing_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  if (!INSTALL_MISSING) {
    stop(
      "Missing CRAN packages: ",
      paste(missing_cran, collapse = ", ")
    )
  }

  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  if (!INSTALL_MISSING) {
    stop("BiocManager is required.")
  }

  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

missing_bioc <- bioc_pkgs[
  !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  if (!INSTALL_MISSING) {
    stop(
      "Missing Bioconductor packages: ",
      paste(missing_bioc, collapse = ", ")
    )
  }

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. CONFIG
# ------------------------------------------------------------------------------

RDS_PATH <- "data/Integration_paravertebral_final.rds"

TF_TABLE <- file.path(
  "results_tf_activity",
  "TF_activity",
  "all_celltypes_CollecTRI_ULM.tsv"
)

OUTDIR <- file.path(
  "results_tf_activity",
  "figures"
)

REDUCTION <- "umap.rpca"
ASSAY <- "RNA"

CELLTYPE_COL <- "cell_type"
CONDITION_COL <- "condition"

CTRL_LABEL <- "CTRL"
DMD_LABEL <- "DMD"

TF_FOCUS <- "JUN"

# Top 10 statistically supported TFs; JUN is always retained for the ANR focus.
TOP_N_TF <- 10

# Display only. Raw ULM scores remain unchanged.
ACTIVITY_DISPLAY_LIMIT <- 4

# Minimum regulon size for per-nucleus exploratory activity inference.
MIN_REGULON_SIZE <- 10

CELLTYPE_ORDER <- c(
  "Myofibers",
  "Regenerating Myofibers",
  "Satellite cells",
  "FAPs",
  "Immune cells",
  "Endothelial",
  "Smooth muscle",
  "Tenocytes",
  "Schwann cells",
  "Neural cells",
  "Adipocytes"
)

CELLTYPE_COLORS <- c(
  "Myofibers" = "#D95F4E",
  "Regenerating Myofibers" = "#F08A84",
  "Satellite cells" = "#59A14F",
  "FAPs" = "#4E79A7",
  "Immune cells" = "#F28E2B",
  "Endothelial" = "#6BAFB0",
  "Smooth muscle" = "#A96F9C",
  "Tenocytes" = "#9C755F",
  "Schwann cells" = "#B07AA1",
  "Neural cells" = "#7BC96F",
  "Adipocytes" = "#E8B92E"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------------------------

message("Reading Seurat object: ", RDS_PATH)
obj <- readRDS(RDS_PATH)

if (!inherits(obj, "Seurat")) {
  stop("RDS_PATH does not contain a Seurat object.")
}

required_meta <- c(
  CELLTYPE_COL,
  CONDITION_COL
)

missing_meta <- setdiff(
  required_meta,
  colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (!file.exists(TF_TABLE)) {
  stop(
    "TF activity table not found: ",
    TF_TABLE,
    "\nRun tf_activity_decoupler.R first."
  )
}

tf_all <- read.delim(
  TF_TABLE,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  dplyr::mutate(
    source = as.character(source),
    celltype = as.character(celltype),
    estimate = as.numeric(estimate),
    padj = as.numeric(padj)
  )

if (!all(c("source", "celltype", "estimate") %in% colnames(tf_all))) {
  stop(
    "TF_TABLE must contain source, celltype and estimate columns."
  )
}

if (!TF_FOCUS %in% unique(tf_all$source)) {
  stop(
    "JUN is not present in the donor-aware TF activity table."
  )
}

# ------------------------------------------------------------------------------
# 3. PANEL A -- CLEAN UMAP
# ------------------------------------------------------------------------------

available_reductions <- Seurat::Reductions(obj)

if (!REDUCTION %in% available_reductions) {

  umap_candidates <- grep(
    "umap",
    available_reductions,
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(umap_candidates) == 0) {
    stop(
      "No UMAP reduction found. Available reductions: ",
      paste(available_reductions, collapse = ", ")
    )
  }

  REDUCTION <- umap_candidates[[1]]

  message(
    "Requested UMAP not found; using: ",
    REDUCTION
  )
}

emb <- as.data.frame(
  Seurat::Embeddings(
    obj,
    reduction = REDUCTION
  )
)

colnames(emb)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

emb$celltype <- as.character(
  obj@meta.data[[CELLTYPE_COL]]
)

present_celltypes <- sort(
  unique(emb$celltype)
)

missing_colors <- setdiff(
  present_celltypes,
  names(CELLTYPE_COLORS)
)

if (length(missing_colors) > 0) {

  extra_cols <- grDevices::hcl.colors(
    length(missing_colors),
    palette = "Dark 3"
  )

  names(extra_cols) <- missing_colors

  CELLTYPE_COLORS <- c(
    CELLTYPE_COLORS,
    extra_cols
  )
}

p_umap <- ggplot2::ggplot(
  emb,
  ggplot2::aes(
    x = UMAP_1,
    y = UMAP_2,
    color = celltype
  )
) +
  ggplot2::geom_point(
    size = 0.20,
    alpha = 0.82,
    stroke = 0
  ) +
  ggplot2::scale_color_manual(
    values = CELLTYPE_COLORS,
    drop = FALSE,
    name = NULL,
    guide = ggplot2::guide_legend(
      override.aes = list(
        size = 3.2,
        alpha = 1
      )
    )
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    tag = "A",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  ggplot2::theme_classic(
    base_size = 9
  ) +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.title = ggplot2::element_text(
      size = 9
    ),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "right",
    legend.text = ggplot2::element_text(
      size = 7.8
    ),
    legend.key.height = grid::unit(
      0.30,
      "cm"
    ),
    legend.key.width = grid::unit(
      0.30,
      "cm"
    ),
    plot.margin = ggplot2::margin(
      6, 6, 6, 6
    )
  )

# ------------------------------------------------------------------------------
# 4. PANEL B -- REPRESENTATIVE DIFFERENTIAL TF ACTIVITIES ACROSS CELL TYPES
# ------------------------------------------------------------------------------

celltype_priority <- CELLTYPE_ORDER[
  CELLTYPE_ORDER %in% unique(tf_all$celltype)
]

extra_celltypes <- setdiff(
  unique(tf_all$celltype),
  celltype_priority
)

celltype_priority <- c(
  celltype_priority,
  sort(extra_celltypes)
)

# Rank TFs independently within each cell population.
# Primary criterion: adjusted P value.
# Secondary criterion: absolute ULM activity effect.
pair_rank <- tf_all %>%
  dplyr::filter(
    is.finite(estimate),
    !is.na(padj),
    padj < 0.05
  ) %>%
  dplyr::mutate(
    abs_activity = abs(estimate),
    celltype = as.character(celltype),
    source = as.character(source)
  ) %>%
  dplyr::group_by(celltype) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs_activity),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    cell_rank = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

# JUN is the regulator of interest. Use its strongest significant cell-type
# association as its representative row. If JUN is not formally significant,
# retain its strongest activity row so the ANR focus remains explicit.
focus_row <- pair_rank %>%
  dplyr::filter(
    source == TF_FOCUS
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs_activity)
  ) %>%
  dplyr::slice_head(
    n = 1
  )

if (nrow(focus_row) == 0) {

  focus_row <- tf_all %>%
    dplyr::filter(
      source == TF_FOCUS,
      is.finite(estimate)
    ) %>%
    dplyr::mutate(
      abs_activity = abs(estimate),
      cell_rank = NA_integer_
    ) %>%
    dplyr::arrange(
      dplyr::desc(abs_activity)
    ) %>%
    dplyr::slice_head(
      n = 1
    )
}

selected_rows <- focus_row
selected_tfs <- unique(
  as.character(selected_rows$source)
)
covered_celltypes <- unique(
  as.character(selected_rows$celltype)
)

# First pass: maximise cell-type coverage by taking the best still-unselected
# significant TF from each population.
for (ct in celltype_priority) {

  if (length(selected_tfs) >= TOP_N_TF) {
    break
  }

  if (ct %in% covered_celltypes) {
    next
  }

  candidate <- pair_rank %>%
    dplyr::filter(
      celltype == ct,
      !source %in% selected_tfs
    ) %>%
    dplyr::slice_head(
      n = 1
    )

  if (nrow(candidate) == 0) {
    next
  }

  selected_rows <- dplyr::bind_rows(
    selected_rows,
    candidate
  )

  selected_tfs <- unique(
    c(
      selected_tfs,
      as.character(candidate$source)
    )
  )

  covered_celltypes <- unique(
    c(
      covered_celltypes,
      ct
    )
  )
}

# Second pass: if fewer than TOP_N_TF unique TFs were obtained because several
# populations share the same top regulator, fill with the best remaining
# within-cell candidates, prioritising rank within each cell population.
if (length(selected_tfs) < TOP_N_TF) {

  n_missing <- TOP_N_TF - length(selected_tfs)

  filler <- pair_rank %>%
    dplyr::filter(
      !source %in% selected_tfs
    ) %>%
    dplyr::arrange(
      cell_rank,
      padj,
      dplyr::desc(abs_activity)
    ) %>%
    dplyr::distinct(
      source,
      .keep_all = TRUE
    ) %>%
    dplyr::slice_head(
      n = n_missing
    )

  selected_rows <- dplyr::bind_rows(
    selected_rows,
    filler
  )

  selected_tfs <- unique(
    c(
      selected_tfs,
      as.character(filler$source)
    )
  )
}

top_tfs <- selected_tfs[
  seq_len(
    min(
      TOP_N_TF,
      length(selected_tfs)
    )
  )
]

# Keep a transparent record of why each TF was selected for the ANR figure.
selection_table <- selected_rows %>%
  dplyr::filter(
    source %in% top_tfs
  ) %>%
  dplyr::distinct(
    source,
    .keep_all = TRUE
  ) %>%
  dplyr::transmute(
    source,
    representative_celltype = as.character(celltype),
    estimate,
    padj,
    within_cell_rank = cell_rank
  )

write.table(
  selection_table,
  file.path(
    OUTDIR,
    "ANR_top10_TF_selection.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Use a biologically readable, fixed cell-type order in the heatmap rather than
# ordering columns by the JUN effect itself.
celltype_order <- celltype_priority

message(
  "ANR TF selection: ",
  paste(
    top_tfs,
    collapse = ", "
  )
)

heat_df <- tf_all %>%
  dplyr::filter(
    source %in% top_tfs,
    celltype %in% celltype_order,
    is.finite(estimate)
  ) %>%
  dplyr::mutate(
    plot_estimate = pmax(
      -ACTIVITY_DISPLAY_LIMIT,
      pmin(
        ACTIVITY_DISPLAY_LIMIT,
        estimate
      )
    ),
    source = factor(
      source,
      levels = rev(top_tfs)
    ),
    celltype = factor(
      celltype,
      levels = celltype_order
    ),
    sig_label = dplyr::case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01 ~ "**",
      !is.na(padj) & padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_heat <- ggplot2::ggplot(
  heat_df,
  ggplot2::aes(
    x = celltype,
    y = source,
    fill = plot_estimate
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sig_label
    ),
    size = 3.0,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -ACTIVITY_DISPLAY_LIMIT,
      ACTIVITY_DISPLAY_LIMIT
    ),
    name = "ULM score"
  ) +
  ggplot2::labs(
    tag = "B",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(
    base_size = 8.5
  ) +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.text.x = ggplot2::element_text(
      angle = 38,
      hjust = 1,
      vjust = 1,
      size = 7.4
    ),
    axis.text.y = ggplot2::element_text(
      size = 8.2,
      face = "bold"
    ),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(
      size = 7.5
    ),
    legend.text = ggplot2::element_text(
      size = 7
    ),
    plot.margin = ggplot2::margin(
      2, 4, 3, 4
    )
  )

# ------------------------------------------------------------------------------
# 5. COLLECTRI JUN REGULON FOR PER-NUCLEUS EXPLORATORY ACTIVITY
# ------------------------------------------------------------------------------

# Avoid an OmnipathR dependency here. In the current conda environment,
# OmnipathR pulls xml2/rvest and can fail to compile because of system zlib
# linkage. CollecTRI provides the signed human regulon directly as a static CSV.
COLLECTRI_URL <- "https://rescued.omnipathdb.org/CollecTRI.csv"
COLLECTRI_CACHE <- file.path(
  OUTDIR,
  "CollecTRI_human.csv"
)

load_collectri_static <- function() {

  if (!file.exists(COLLECTRI_CACHE)) {

    message(
      "Downloading CollecTRI static network..."
    )

    ok <- tryCatch(
      {
        utils::download.file(
          url = COLLECTRI_URL,
          destfile = COLLECTRI_CACHE,
          mode = "wb",
          quiet = FALSE
        )
        TRUE
      },
      error = function(e) {
        message(
          "CollecTRI download failed: ",
          conditionMessage(e)
        )
        FALSE
      }
    )

    if (!ok || !file.exists(COLLECTRI_CACHE)) {
      stop(
        "Could not download CollecTRI from: ",
        COLLECTRI_URL
      )
    }
  }

  raw <- utils::read.csv(
    COLLECTRI_CACHE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "source",
    "target",
    "weight"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(raw)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Unexpected CollecTRI static CSV format. Missing: ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ",
      paste(colnames(raw), collapse = ", ")
    )
  }

  raw %>%
    dplyr::transmute(
      source = as.character(source),
      target = as.character(target),
      mor = as.numeric(weight)
    ) %>%
    dplyr::filter(
      !is.na(source),
      !is.na(target),
      is.finite(mor),
      source != "",
      target != ""
    ) %>%
    dplyr::distinct(
      source,
      target,
      .keep_all = TRUE
    )
}

collectri <- load_collectri_static()

message(
  "CollecTRI loaded: ",
  nrow(collectri),
  " interactions; ",
  dplyr::n_distinct(collectri$source),
  " regulators."
)

focus_net <- collectri %>%
  dplyr::filter(
    source == TF_FOCUS
  )

if (nrow(focus_net) < MIN_REGULON_SIZE) {
  stop(
    "JUN regulon is too small after loading CollecTRI: ",
    nrow(focus_net),
    " interactions."
  )
}

# ------------------------------------------------------------------------------
# 6. PANEL C -- JUN PER-NUCLEUS ACTIVITY VIOLIN
# ------------------------------------------------------------------------------

DefaultAssay(obj) <- ASSAY

expr <- Seurat::GetAssayData(
  obj,
  assay = ASSAY,
  layer = "data"
)

focus_targets <- intersect(
  rownames(expr),
  unique(focus_net$target)
)

focus_net_use <- focus_net %>%
  dplyr::filter(
    target %in% focus_targets
  )

if (length(focus_targets) < MIN_REGULON_SIZE) {
  stop(
    "Too few JUN targets overlap the RNA assay: ",
    length(focus_targets)
  )
}

message(
  "Inferring exploratory per-nucleus JUN activity from ",
  length(focus_targets),
  " targets..."
)

# Restrict to JUN targets before inference to keep memory use manageable.
expr_focus <- expr[
  focus_targets,
  ,
  drop = FALSE
]

focus_cell <- decoupleR::run_ulm(
  mat = expr_focus,
  network = focus_net_use,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = MIN_REGULON_SIZE
)

if ("condition" %in% colnames(focus_cell) &&
    !"sample" %in% colnames(focus_cell)) {
  focus_cell <- focus_cell %>%
    dplyr::rename(
      sample = condition
    )
}

if ("score" %in% colnames(focus_cell) &&
    !"estimate" %in% colnames(focus_cell)) {
  focus_cell <- focus_cell %>%
    dplyr::rename(
      estimate = score
    )
}

if (!all(c("sample", "source", "estimate") %in% colnames(focus_cell))) {
  stop(
    "Unexpected run_ulm() per-nucleus output columns: ",
    paste(
      colnames(focus_cell),
      collapse = ", "
    )
  )
}

cell_meta <- obj@meta.data %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  dplyr::transmute(
    cell_id = as.character(cell_id),
    celltype = as.character(.data[[CELLTYPE_COL]]),
    condition = as.character(.data[[CONDITION_COL]])
  )

violin_df <- focus_cell %>%
  dplyr::filter(
    source == TF_FOCUS
  ) %>%
  dplyr::mutate(
    cell_id = as.character(sample),
    estimate = as.numeric(estimate)
  ) %>%
  dplyr::select(
    -sample
  ) %>%
  dplyr::inner_join(
    cell_meta,
    by = "cell_id"
  ) %>%
  dplyr::filter(
    condition %in% c(
      CTRL_LABEL,
      DMD_LABEL
    ),
    celltype %in% celltype_order,
    is.finite(estimate)
  ) %>%
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = celltype_order
    ),
    condition = factor(
      condition,
      levels = c(
        CTRL_LABEL,
        DMD_LABEL
      )
    )
  )

# Winsorise for visualization only so extreme nuclei do not flatten violins.
violin_limits <- stats::quantile(
  violin_df$estimate,
  probs = c(
    0.01,
    0.99
  ),
  na.rm = TRUE
)

violin_df <- violin_df %>%
  dplyr::mutate(
    plot_estimate = pmax(
      violin_limits[[1]],
      pmin(
        violin_limits[[2]],
        estimate
      )
    )
  )

p_violin <- ggplot2::ggplot(
  violin_df,
  ggplot2::aes(
    x = celltype,
    y = plot_estimate,
    fill = condition,
    group = interaction(
      celltype,
      condition
    )
  )
) +
  ggplot2::geom_violin(
    position = ggplot2::position_dodge(
      width = 0.82
    ),
    scale = "width",
    trim = TRUE,
    color = NA,
    alpha = 0.88
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(
      group = condition
    ),
    fun = stats::median,
    geom = "point",
    position = ggplot2::position_dodge(
      width = 0.82
    ),
    size = 1.15,
    color = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "CTRL" = "#BDBDBD",
      "DMD" = "#D95F4E"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    tag = "C",
    x = NULL,
    y = "JUN activity score"
  ) +
  ggplot2::theme_classic(
    base_size = 8.5
  ) +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.text.x = ggplot2::element_text(
      angle = 38,
      hjust = 1,
      size = 7.2
    ),
    axis.text.y = ggplot2::element_text(
      size = 7
    ),
    axis.title.y = ggplot2::element_text(
      size = 7.8
    ),
    legend.position = "top",
    legend.text = ggplot2::element_text(
      size = 7.8
    ),
    legend.key.width = grid::unit(
      0.45,
      "cm"
    ),
    plot.margin = ggplot2::margin(
      3, 4, 3, 4
    )
  )

# ------------------------------------------------------------------------------
# 7. ASSEMBLE
# ------------------------------------------------------------------------------

right_panel <- (
  p_heat /
    p_violin
) +
  patchwork::plot_layout(
    heights = c(
      1.05,
      0.95
    )
  )

final_plot <- (
  p_umap |
    right_panel
) +
  patchwork::plot_layout(
    widths = c(
      1.05,
      1.45
    )
  ) &
  ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    )
  )

# ------------------------------------------------------------------------------
# 8. EXPORT
# ------------------------------------------------------------------------------

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(OUTDIR)) {
  stop(
    "Could not create output directory: ",
    OUTDIR
  )
}

base_name <- file.path(
  OUTDIR,
  "ANR_preliminary_snRNAseq_TF_activity"
)

safe_output_path <- function(path) {

  if (!file.exists(path)) {
    return(path)
  }

  removed <- suppressWarnings(
    file.remove(path)
  )

  if (isTRUE(removed)) {
    return(path)
  }

  ext <- tools::file_ext(path)
  stem <- sub(
    paste0("\\.", ext, "$"),
    "",
    path
  )

  fallback <- paste0(
    stem,
    "_new.",
    ext
  )

  message(
    "Could not overwrite existing file (possibly open/locked): ",
    path,
    "\nWriting instead to: ",
    fallback
  )

  fallback
}

pdf_path <- safe_output_path(
  paste0(
    base_name,
    ".pdf"
  )
)

png_path <- safe_output_path(
  paste0(
    base_name,
    ".png"
  )
)

svg_path <- safe_output_path(
  paste0(
    base_name,
    ".svg"
  )
)

ggplot2::ggsave(
  pdf_path,
  final_plot,
  width = 12.5,
  height = 5.4,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  png_path,
  final_plot,
  width = 12.5,
  height = 5.4,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggplot2::ggsave(
  svg_path,
  final_plot,
  width = 12.5,
  height = 5.4,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)

message(
  "ANR figure written to:",
  "\n  PDF: ", pdf_path,
  "\n  PNG: ", png_path,
  "\n  SVG: ", svg_path
)
