#!/usr/bin/env Rscript

# DMD snRNA-seq: donor-aware TF activity analysis
# Dataset: 3 DMD vs 2 CTRL
#
# Workflow:
#   Seurat RDS
#     -> donor x cell-type pseudobulk
#     -> DESeq2 (DMD vs CTRL)
#     -> CollecTRI TF regulons
#     -> decoupleR ULM on DESeq2 Wald statistics
#     -> per-cell-type TF activity tables + global heatmap
#
# Biological replicates are DONORS, not nuclei.

# ==============================================================================
# 0. DEPENDENCIES
# ==============================================================================

INSTALL_MISSING <- TRUE

cran_pkgs <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tidyr",
  "ggplot2",
  "tibble"
)

bioc_pkgs <- c(
  "DESeq2",
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

  message(
    "Installing missing CRAN packages: ",
    paste(missing_cran, collapse = ", ")
  )

  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  if (!INSTALL_MISSING) {
    stop("BiocManager is required to install Bioconductor packages.")
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

  message(
    "Installing missing Bioconductor packages: ",
    paste(missing_bioc, collapse = ", ")
  )

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(decoupleR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ==============================================================================
# 1. CONFIG
# ==============================================================================

# Your current file
RDS_PATH <- "data/Integration_paravertebral_final.rds"

ASSAY <- "RNA"

# Metadata columns present in Integration_paravertebral_final.rds
DONOR_COL <- "sample"
CONDITION_COL <- "condition"
CELLTYPE_COL <- "cell_type"

CTRL_LABEL <- "CTRL"
DMD_LABEL <- "DMD"

# Minimum nuclei required for one donor x cell-type pseudobulk
MIN_NUCLEI_PER_DONOR_CELLTYPE <- 20

# With 3 DMD vs 2 CTRL, require at least 2 donors in each condition
MIN_DONORS_PER_CONDITION <- 2

# Gene filtering before DESeq2
MIN_COUNT <- 10
MIN_SAMPLES_WITH_COUNT <- 2

# CollecTRI / decoupleR
MIN_REGULON_SIZE <- 10

# Number of TFs displayed in the summary heatmap
TOP_TF_HEATMAP <- 30

OUTDIR <- "results_tf_activity"

# ==============================================================================
# 2. LOAD SEURAT OBJECT
# ==============================================================================

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "DESeq2"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "TF_activity"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "figures"), recursive = TRUE, showWarnings = FALSE)

message("Reading: ", RDS_PATH)
obj <- readRDS(RDS_PATH)

if (!inherits(obj, "Seurat")) {
  stop("RDS_PATH does not contain a Seurat object.")
}

required_meta <- c(DONOR_COL, CONDITION_COL, CELLTYPE_COL)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", "),
    "\nAvailable metadata columns: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  )
}

if (!ASSAY %in% Assays(obj)) {
  stop("Assay '", ASSAY, "' not found in the Seurat object.")
}

DefaultAssay(obj) <- ASSAY

meta <- obj@meta.data %>%
  dplyr::mutate(
    .cell = rownames(obj@meta.data),
    donor = as.character(.data[[DONOR_COL]]),
    condition = as.character(.data[[CONDITION_COL]]),
    celltype = as.character(.data[[CELLTYPE_COL]])
  ) %>%
  dplyr::filter(
    !is.na(donor),
    !is.na(condition),
    !is.na(celltype),
    donor != "",
    condition != "",
    celltype != ""
  )

message(
  "Samples: ",
  paste(sort(unique(meta$donor)), collapse = ", ")
)

message(
  "Conditions: ",
  paste(sort(unique(meta$condition)), collapse = ", ")
)

message(
  "Cell types: ",
  paste(sort(unique(meta$celltype)), collapse = ", ")
)

if (!all(c(CTRL_LABEL, DMD_LABEL) %in% unique(meta$condition))) {
  stop(
    "CTRL_LABEL/DMD_LABEL do not match CONDITION_COL values. Found: ",
    paste(sort(unique(meta$condition)), collapse = ", ")
  )
}

sample_condition <- meta %>%
  dplyr::distinct(donor, condition) %>%
  dplyr::count(donor, name = "n_conditions")

if (any(sample_condition$n_conditions != 1)) {
  stop("At least one donor is associated with more than one condition.")
}

# ==============================================================================
# 3. QC: NUCLEI PER DONOR x CELL TYPE
# ==============================================================================

qc <- meta %>%
  dplyr::count(celltype, donor, condition, name = "n_nuclei") %>%
  dplyr::arrange(celltype, condition, donor)

write.table(
  qc,
  file.path(OUTDIR, "nuclei_per_donor_celltype.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

eligible <- qc %>%
  dplyr::filter(n_nuclei >= MIN_NUCLEI_PER_DONOR_CELLTYPE) %>%
  dplyr::count(celltype, condition, name = "n_donors") %>%
  dplyr::filter(condition %in% c(CTRL_LABEL, DMD_LABEL)) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = n_donors,
    values_fill = 0
  )

if (!CTRL_LABEL %in% colnames(eligible)) eligible[[CTRL_LABEL]] <- 0
if (!DMD_LABEL %in% colnames(eligible)) eligible[[DMD_LABEL]] <- 0

eligible_celltypes <- eligible %>%
  dplyr::filter(
    .data[[CTRL_LABEL]] >= MIN_DONORS_PER_CONDITION,
    .data[[DMD_LABEL]] >= MIN_DONORS_PER_CONDITION
  ) %>%
  dplyr::pull(celltype)

write.table(
  eligible,
  file.path(OUTDIR, "celltype_eligibility.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (length(eligible_celltypes) == 0) {
  stop(
    "No cell type passes the donor/nuclei thresholds. ",
    "Inspect results_tf_activity/nuclei_per_donor_celltype.tsv."
  )
}

message(
  "Eligible cell types (", length(eligible_celltypes), "): ",
  paste(eligible_celltypes, collapse = ", ")
)

# ==============================================================================
# 4. COLLECTRI NETWORK
# ==============================================================================

message("Loading CollecTRI...")

collectri <- decoupleR::get_collectri(
  organism = "human",
  split_complexes = FALSE
)

# Support common decoupleR/CollecTRI column naming variants
if ("weight" %in% colnames(collectri) && !"mor" %in% colnames(collectri)) {
  collectri <- collectri %>%
    dplyr::rename(mor = weight)
}

required_net_cols <- c("source", "target", "mor")
missing_net_cols <- setdiff(required_net_cols, colnames(collectri))

if (length(missing_net_cols) > 0) {
  stop(
    "Unexpected CollecTRI format. Missing columns: ",
    paste(missing_net_cols, collapse = ", "),
    "\nAvailable columns: ",
    paste(colnames(collectri), collapse = ", ")
  )
}

# ==============================================================================
# 5. RAW COUNTS
# ==============================================================================

counts <- Seurat::GetAssayData(
  obj,
  assay = ASSAY,
  layer = "counts"
)

if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

if (nrow(counts) == 0 || ncol(counts) == 0) {
  stop("RNA counts layer is empty.")
}

network_overlap <- length(
  intersect(rownames(counts), unique(collectri$target))
)

message(
  "Gene overlap with CollecTRI targets: ",
  network_overlap
)

if (network_overlap < 100) {
  stop(
    "Very low overlap between RNA feature names and CollecTRI gene symbols. ",
    "Check whether rownames(obj) are Ensembl IDs rather than gene symbols."
  )
}

# ==============================================================================
# 6. PER CELL TYPE:
#    PSEUDOBULK -> DESEQ2 -> COLLECTRI / ULM
# ==============================================================================

all_de <- list()
all_tf <- list()

for (ct in eligible_celltypes) {

  message("\n============================================================")
  message("Cell type: ", ct)
  message("============================================================")

  ct_meta <- meta %>%
    dplyr::filter(
      celltype == ct,
      condition %in% c(CTRL_LABEL, DMD_LABEL)
    )

  valid_donors <- ct_meta %>%
    dplyr::count(donor, condition, name = "n_nuclei") %>%
    dplyr::filter(n_nuclei >= MIN_NUCLEI_PER_DONOR_CELLTYPE)

  ct_meta <- ct_meta %>%
    dplyr::semi_join(
      valid_donors,
      by = c("donor", "condition")
    )

  if (nrow(ct_meta) == 0) {
    warning("No usable nuclei for cell type: ", ct)
    next
  }

  donor_levels <- unique(ct_meta$donor)

  donor_factor <- factor(
    ct_meta$donor,
    levels = donor_levels
  )

  # Sparse nuclei -> donor aggregation matrix
  design <- Matrix::sparse.model.matrix(
    ~ 0 + donor_factor
  )

  colnames(design) <- donor_levels

  ct_counts <- counts[
    ,
    ct_meta$.cell,
    drop = FALSE
  ]

  pseudobulk <- ct_counts %*% design

  donor_meta <- ct_meta %>%
    dplyr::distinct(donor, condition)

  donor_meta <- donor_meta[
    match(colnames(pseudobulk), donor_meta$donor),
    ,
    drop = FALSE
  ]

  rownames(donor_meta) <- donor_meta$donor

  donor_meta$condition <- factor(
    donor_meta$condition,
    levels = c(CTRL_LABEL, DMD_LABEL)
  )

  n_ctrl <- sum(donor_meta$condition == CTRL_LABEL)
  n_dmd <- sum(donor_meta$condition == DMD_LABEL)

  if (
    n_ctrl < MIN_DONORS_PER_CONDITION ||
    n_dmd < MIN_DONORS_PER_CONDITION
  ) {
    warning(
      "Skipping ", ct,
      ": only ", n_ctrl, " CTRL and ", n_dmd, " DMD donors remain."
    )
    next
  }

  # Keep genes expressed at a minimal count in at least N donor pseudobulks
  keep_genes <- Matrix::rowSums(
    pseudobulk >= MIN_COUNT
  ) >= MIN_SAMPLES_WITH_COUNT

  pseudobulk <- pseudobulk[
    keep_genes,
    ,
    drop = FALSE
  ]

  message(
    "Donors: ",
    n_ctrl, " CTRL / ",
    n_dmd, " DMD"
  )

  message(
    "Genes retained: ",
    nrow(pseudobulk)
  )

  if (nrow(pseudobulk) == 0) {
    warning("No genes retained for cell type: ", ct)
    next
  }

  # DESeq2 expects integer counts.
  # Rounding is only relevant if the stored counts layer contains fractional
  # values from an upstream correction procedure.
  pb_dense <- round(as.matrix(pseudobulk))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = pb_dense,
    colData = as.data.frame(donor_meta),
    design = ~ condition
  )

  dds <- DESeq2::DESeq(
    dds,
    quiet = TRUE
  )

  contrast_name <- paste0(
    "condition_",
    DMD_LABEL,
    "_vs_",
    CTRL_LABEL
  )

  if (!contrast_name %in% DESeq2::resultsNames(dds)) {
    stop(
      "DESeq2 contrast not found for ", ct,
      ". Available coefficients: ",
      paste(DESeq2::resultsNames(dds), collapse = ", ")
    )
  }

  res <- DESeq2::results(
    dds,
    name = contrast_name
  )

  de <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      celltype = ct,
      padj = ifelse(is.na(padj), 1, padj)
    ) %>%
    dplyr::select(
      celltype,
      gene,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    ) %>%
    dplyr::arrange(
      padj,
      dplyr::desc(abs(stat))
    )

  safe_ct <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    ct
  )

  write.table(
    de,
    file.path(
      OUTDIR,
      "DESeq2",
      paste0(
        safe_ct,
        "_DMD_vs_CTRL.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  all_de[[ct]] <- de

  # --------------------------------------------------------------------------
  # CollecTRI TF activity
  #
  # We use the DESeq2 Wald statistic rather than log2FC alone.
  # Positive score = regulon more active in DMD relative to CTRL.
  # Negative score = regulon less active in DMD relative to CTRL.
  # --------------------------------------------------------------------------

  stat_tbl <- de %>%
    dplyr::filter(
      is.finite(stat),
      !is.na(gene)
    ) %>%
    dplyr::distinct(
      gene,
      .keep_all = TRUE
    )

  stat_mat <- matrix(
    stat_tbl$stat,
    ncol = 1,
    dimnames = list(
      stat_tbl$gene,
      "DMD_vs_CTRL"
    )
  )

  tf <- decoupleR::run_ulm(
    mat = stat_mat,
    network = collectri,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = MIN_REGULON_SIZE
  )

  # Support output naming differences across decoupleR versions
  if ("condition" %in% colnames(tf) && !"sample" %in% colnames(tf)) {
    tf <- tf %>%
      dplyr::rename(sample = condition)
  }

  if ("score" %in% colnames(tf) && !"estimate" %in% colnames(tf)) {
    tf <- tf %>%
      dplyr::rename(estimate = score)
  }

  if ("p_value" %in% colnames(tf) && !"pvalue" %in% colnames(tf)) {
    tf <- tf %>%
      dplyr::rename(pvalue = p_value)
  }

  if (!all(c("source", "estimate") %in% colnames(tf))) {
    stop(
      "Unexpected run_ulm() output columns: ",
      paste(colnames(tf), collapse = ", ")
    )
  }

  if (!"pvalue" %in% colnames(tf)) {
    tf$pvalue <- NA_real_
  }

  tf <- tf %>%
    dplyr::mutate(
      celltype = ct,
      padj = ifelse(
        is.na(pvalue),
        NA_real_,
        p.adjust(
          pvalue,
          method = "BH"
        )
      )
    ) %>%
    dplyr::select(
      celltype,
      source,
      estimate,
      pvalue,
      padj,
      dplyr::everything()
    ) %>%
    dplyr::arrange(
      dplyr::desc(abs(estimate))
    )

  write.table(
    tf,
    file.path(
      OUTDIR,
      "TF_activity",
      paste0(
        safe_ct,
        "_CollecTRI_ULM.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  all_tf[[ct]] <- tf

  rm(
    dds,
    res,
    pseudobulk,
    pb_dense,
    ct_counts,
    design
  )

  gc()
}

# ==============================================================================
# 7. COMBINED OUTPUTS
# ==============================================================================

if (length(all_de) == 0) {
  stop("No DESeq2 result was generated.")
}

if (length(all_tf) == 0) {
  stop("No TF activity result was generated.")
}

de_all <- dplyr::bind_rows(all_de)
tf_all <- dplyr::bind_rows(all_tf)

write.table(
  de_all,
  file.path(
    OUTDIR,
    "DESeq2",
    "all_celltypes_DMD_vs_CTRL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tf_all,
  file.path(
    OUTDIR,
    "TF_activity",
    "all_celltypes_CollecTRI_ULM.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================================
# 8. GLOBAL TF HEATMAP
# ==============================================================================

top_tfs <- tf_all %>%
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
  ) %>%
  dplyr::slice_head(
    n = min(
      TOP_TF_HEATMAP,
      dplyr::n()
    )
  ) %>%
  dplyr::pull(source)

heat_df <- tf_all %>%
  dplyr::filter(
    source %in% top_tfs
  ) %>%
  dplyr::mutate(
    source = factor(
      source,
      levels = rev(top_tfs)
    ),
    celltype = factor(
      celltype,
      levels = eligible_celltypes
    )
  )

p <- ggplot2::ggplot(
  heat_df,
  ggplot2::aes(
    x = celltype,
    y = source,
    fill = estimate
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    midpoint = 0,
    name = "TF activity\nULM score"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = "Differential TF activity in DMD",
    subtitle = "CollecTRI + decoupleR ULM on DESeq2 Wald statistics"
  ) +
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    OUTDIR,
    "figures",
    "TF_activity_heatmap_top30.pdf"
  ),
  p,
  width = max(
    7,
    0.6 * length(eligible_celltypes) + 4
  ),
  height = 10
)

ggplot2::ggsave(
  file.path(
    OUTDIR,
    "figures",
    "TF_activity_heatmap_top30.png"
  ),
  p,
  width = max(
    7,
    0.6 * length(eligible_celltypes) + 4
  ),
  height = 10,
  dpi = 300
)

# ==============================================================================
# 9. SESSION INFO
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

message("\nDone.")
message(
  "Results: ",
  normalizePath(OUTDIR)
)
