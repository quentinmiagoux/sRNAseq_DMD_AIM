#!/usr/bin/env Rscript

# ==============================================================================
# ANR preliminary-data figure: DMD snRNA-seq + TF activity
#
# Panel A  : real snRNA-seq UMAP
# Panel B1 : AP1 + top candidate TF activity heatmap (donor-aware pseudobulk)
# Panel B2 : per-nucleus AP1 activity distribution across DMD cell populations
#
# Statistical interpretation:
# - The heatmap is the inferential result: donor-aware pseudobulk DESeq2
#   followed by CollecTRI + decoupleR ULM.
# - The violin plot is an exploratory per-nucleus visualization of AP1 activity.
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

TF_FOCUS <- "AP1"

# AP1 + four contextual candidate TFs
TOP_N_TF <- 5

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
    "AP1 is not present in the donor-aware TF activity table."
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

centroids <- emb %>%
  dplyr::filter(
    !is.na(celltype),
    celltype != ""
  ) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    UMAP_1 = stats::median(
      UMAP_1,
      na.rm = TRUE
    ),
    UMAP_2 = stats::median(
      UMAP_2,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p_umap <- ggplot2::ggplot(
  emb,
  ggplot2::aes(
    x = UMAP_1,
    y = UMAP_2,
    color = celltype
  )
) +
  ggplot2::geom_point(
    size = 0.16,
    alpha = 0.65,
    stroke = 0
  ) +
  ggplot2::geom_label(
    data = centroids,
    ggplot2::aes(
      x = UMAP_1,
      y = UMAP_2,
      label = celltype
    ),
    inherit.aes = FALSE,
    size = 3.1,
    fontface = "bold",
    label.size = 0,
    fill = grDevices::adjustcolor(
      "white",
      alpha.f = 0.78
    ),
    color = "black",
    label.padding = grid::unit(
      0.08,
      "lines"
    )
  ) +
  ggplot2::scale_color_manual(
    values = CELLTYPE_COLORS,
    drop = FALSE
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "A  Previous snRNA-seq analysis identifies\ncell populations in DMD muscle",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  ggplot2::theme_classic(
    base_size = 9
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 11.5,
      hjust = 0
    ),
    axis.title = ggplot2::element_text(
      size = 9
    ),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      6, 6, 6, 6
    )
  )

# ------------------------------------------------------------------------------
# 4. PANEL B1 -- TOP 5 DIFFERENTIAL TF ACTIVITIES
# ------------------------------------------------------------------------------

tf_rank <- tf_all %>%
  dplyr::filter(
    is.finite(estimate)
  ) %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(
    max_abs_activity = max(
      abs(estimate),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(max_abs_activity)
  )

context_tfs <- tf_rank %>%
  dplyr::filter(
    source != TF_FOCUS
  ) %>%
  dplyr::slice_head(
    n = TOP_N_TF - 1
  ) %>%
  dplyr::pull(source) %>%
  as.character()

top_tfs <- c(
  TF_FOCUS,
  context_tfs
)

# Order populations by AP1 differential activity to make the signal obvious.
celltype_order <- tf_all %>%
  dplyr::filter(
    source == TF_FOCUS,
    is.finite(estimate)
  ) %>%
  dplyr::arrange(
    dplyr::desc(estimate)
  ) %>%
  dplyr::pull(celltype) %>%
  unique() %>%
  as.character()

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
    title = "Candidate regulator activity across cell populations",
    subtitle = "AP1 highlighted with four contextual TFs | DMD vs CTRL",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(
    base_size = 8.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.8
    ),
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
# 5. COLLECTRI AP1 REGULON FOR PER-NUCLEUS EXPLORATORY ACTIVITY
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

ap1_net <- collectri %>%
  dplyr::filter(
    source == TF_FOCUS
  )

if (nrow(ap1_net) < MIN_REGULON_SIZE) {
  stop(
    "AP1 regulon is too small after loading CollecTRI: ",
    nrow(ap1_net),
    " interactions."
  )
}

# ------------------------------------------------------------------------------
# 6. PANEL B2 -- AP1 PER-NUCLEUS ACTIVITY VIOLIN
# ------------------------------------------------------------------------------

DefaultAssay(obj) <- ASSAY

expr <- Seurat::GetAssayData(
  obj,
  assay = ASSAY,
  layer = "data"
)

ap1_targets <- intersect(
  rownames(expr),
  unique(ap1_net$target)
)

ap1_net_use <- ap1_net %>%
  dplyr::filter(
    target %in% ap1_targets
  )

if (length(ap1_targets) < MIN_REGULON_SIZE) {
  stop(
    "Too few AP1 targets overlap the RNA assay: ",
    length(ap1_targets)
  )
}

message(
  "Inferring exploratory per-nucleus AP1 activity from ",
  length(ap1_targets),
  " targets..."
)

# Restrict to AP1 targets before inference to keep memory use manageable.
expr_ap1 <- expr[
  ap1_targets,
  ,
  drop = FALSE
]

ap1_cell <- decoupleR::run_ulm(
  mat = expr_ap1,
  network = ap1_net_use,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = MIN_REGULON_SIZE
)

if ("condition" %in% colnames(ap1_cell) &&
    !"sample" %in% colnames(ap1_cell)) {
  ap1_cell <- ap1_cell %>%
    dplyr::rename(
      sample = condition
    )
}

if ("score" %in% colnames(ap1_cell) &&
    !"estimate" %in% colnames(ap1_cell)) {
  ap1_cell <- ap1_cell %>%
    dplyr::rename(
      estimate = score
    )
}

if (!all(c("sample", "source", "estimate") %in% colnames(ap1_cell))) {
  stop(
    "Unexpected run_ulm() per-nucleus output columns: ",
    paste(
      colnames(ap1_cell),
      collapse = ", "
    )
  )
}

cell_meta <- obj@meta.data %>%
  tibble::rownames_to_column(
    "sample"
  ) %>%
  dplyr::transmute(
    sample = as.character(sample),
    celltype = as.character(.data[[CELLTYPE_COL]]),
    condition = as.character(.data[[CONDITION_COL]])
  )

violin_df <- ap1_cell %>%
  dplyr::filter(
    source == TF_FOCUS
  ) %>%
  dplyr::mutate(
    sample = as.character(sample),
    estimate = as.numeric(estimate)
  ) %>%
  dplyr::inner_join(
    cell_meta,
    by = "sample"
  ) %>%
  dplyr::filter(
    condition == DMD_LABEL,
    celltype %in% celltype_order,
    is.finite(estimate)
  ) %>%
  dplyr::mutate(
    celltype = factor(
      celltype,
      levels = celltype_order
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
    fill = celltype
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    color = "white",
    linewidth = 0.25,
    alpha = 0.9
  ) +
  ggplot2::stat_summary(
    fun = stats::median,
    geom = "point",
    size = 1.2,
    color = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = CELLTYPE_COLORS,
    guide = "none"
  ) +
  ggplot2::labs(
    title = "AP1 activity across DMD cell populations",
    subtitle = "Per-nucleus CollecTRI/ULM activity (exploratory visualization)",
    x = NULL,
    y = "AP1 activity score"
  ) +
  ggplot2::theme_classic(
    base_size = 8.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.6
    ),
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
    plot.margin = ggplot2::margin(
      3, 4, 3, 4
    )
  )

# ------------------------------------------------------------------------------
# 7. ASSEMBLE
# ------------------------------------------------------------------------------

panel_b <- (
  p_heat /
    p_violin
) +
  patchwork::plot_layout(
    heights = c(
      0.95,
      1.05
    )
  ) +
  patchwork::plot_annotation(
    title = "B  Preliminary analyses identify AP1 as a candidate regulator",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11.5,
        hjust = 0
      )
    )
  )

final_plot <- (
  p_umap |
    panel_b
) +
  patchwork::plot_layout(
    widths = c(
      1.0,
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

base_name <- file.path(
  OUTDIR,
  "ANR_preliminary_snRNAseq_TF_activity"
)

ggplot2::ggsave(
  paste0(
    base_name,
    ".pdf"
  ),
  final_plot,
  width = 12.5,
  height = 5.4,
  units = "in"
)

ggplot2::ggsave(
  paste0(
    base_name,
    ".png"
  ),
  final_plot,
  width = 12.5,
  height = 5.4,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggplot2::ggsave(
  paste0(
    base_name,
    ".svg"
  ),
  final_plot,
  width = 12.5,
  height = 5.4,
  units = "in",
  device = grDevices::svg
)

message(
  "ANR figure written to: ",
  normalizePath(OUTDIR)
)
