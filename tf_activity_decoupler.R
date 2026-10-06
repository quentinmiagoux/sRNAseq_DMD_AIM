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
  "decoupleR",
  "OmnipathR"
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

# Number of strongest TF activities displayed per cell type
TOP_TF_BARPLOT <- 10

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

# decoupleR::get_collectri() currently fails on some recent OmnipathR versions
# when OmniPath falls back to its static table:
#   Error in if (.keep) . else select(., -!!evs_col)
#
# We first try the official decoupleR accessor. If it fails, we use the same
# OmniPath static CollecTRI table with strict_evidences = FALSE and reproduce
# the CollecTRI formatting performed internally by decoupleR.

load_collectri_safe <- function() {

  message("Loading CollecTRI...")

  net <- tryCatch(
    {
      decoupleR::get_collectri(
        organism = "human",
        split_complexes = FALSE
      )
    },
    error = function(e) {

      message(
        "decoupleR::get_collectri() failed: ",
        conditionMessage(e)
      )

      message(
        "Using OmniPath static CollecTRI fallback ",
        "(strict_evidences = FALSE)..."
      )

      raw <- OmnipathR::static_table(
        query = "interactions",
        resource = "collectri",
        organism = 9606L,
        strict_evidences = FALSE,
        wide = FALSE
      )

      required_raw <- c(
        "source",
        "source_genesymbol",
        "target_genesymbol",
        "is_stimulation",
        "is_inhibition"
      )

      missing_raw <- setdiff(
        required_raw,
        colnames(raw)
      )

      if (length(missing_raw) > 0) {
        stop(
          "Unexpected OmniPath CollecTRI static-table format. Missing: ",
          paste(missing_raw, collapse = ", "),
          "\nAvailable columns: ",
          paste(colnames(raw), collapse = ", ")
        )
      }

      cols <- c(
        "source_genesymbol",
        "target_genesymbol",
        "is_stimulation",
        "is_inhibition"
      )

      is_complex <- grepl(
        "COMPLEX",
        raw$source,
        fixed = TRUE
      )

      interactions <- raw[
        !is_complex,
        cols,
        drop = FALSE
      ]

      complexes <- raw[
        is_complex,
        cols,
        drop = FALSE
      ]

      # Same complex handling used by decoupleR when split_complexes = FALSE:
      # AP-1 family complexes -> AP1
      # NF-kB/REL family complexes -> NFKB
      if (nrow(complexes) > 0) {
        complexes$source_genesymbol <- ifelse(
          grepl(
            "JUN|FOS",
            complexes$source_genesymbol
          ),
          "AP1",
          ifelse(
            grepl(
              "REL|NFKB",
              complexes$source_genesymbol
            ),
            "NFKB",
            complexes$source_genesymbol
          )
        )
      }

      net <- dplyr::bind_rows(
        interactions,
        complexes
      ) %>%
        dplyr::distinct(
          source_genesymbol,
          target_genesymbol,
          .keep_all = TRUE
        ) %>%
        dplyr::mutate(
          mor = dplyr::case_when(
            is_stimulation == 1 ~ 1,
            is_stimulation == 0 ~ -1,
            TRUE ~ NA_real_
          )
        ) %>%
        dplyr::transmute(
          source = source_genesymbol,
          target = target_genesymbol,
          mor = mor
        ) %>%
        dplyr::filter(
          !is.na(source),
          !is.na(target),
          !is.na(mor),
          source != "",
          target != ""
        )

      net
    }
  )

  required_net_cols <- c(
    "source",
    "target",
    "mor"
  )

  missing_net_cols <- setdiff(
    required_net_cols,
    colnames(net)
  )

  if (length(missing_net_cols) > 0) {
    stop(
      "Unexpected CollecTRI format. Missing columns: ",
      paste(missing_net_cols, collapse = ", "),
      "\nAvailable columns: ",
      paste(colnames(net), collapse = ", ")
    )
  }

  net <- net %>%
    dplyr::select(
      source,
      target,
      mor
    ) %>%
    dplyr::distinct()

  message(
    "CollecTRI loaded: ",
    nrow(net),
    " TF-target interactions; ",
    dplyr::n_distinct(net$source),
    " regulators."
  )

  net
}

collectri <- load_collectri_safe()

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
tf_all <- dplyr::bind_rows(all_tf) %>%
  dplyr::mutate(
    source = as.character(source),
    celltype = as.character(celltype),
    estimate = as.numeric(estimate)
  )

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
    n = TOP_TF_HEATMAP
  ) %>%
  dplyr::pull(source) %>%
  as.character()

heat_df <- tf_all %>%
  dplyr::filter(
    source %in% top_tfs
  ) %>%
  dplyr::mutate(
    source = factor(
      as.character(source),
      levels = rev(as.character(top_tfs))
    ),
    celltype = factor(
      as.character(celltype),
      levels = as.character(eligible_celltypes)
    ),
    plot_estimate = pmax(
      -4,
      pmin(4, as.numeric(estimate))
    )
  )

p <- ggplot2::ggplot(
  heat_df,
  ggplot2::aes(
    x = celltype,
    y = source,
    fill = plot_estimate
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-4, 4),
    name = "TF activity\nULM score"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = "Differential TF activity in DMD",
    subtitle = "CollecTRI + decoupleR ULM on DESeq2 Wald statistics | display capped at +/-4"
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
# 9. TOP 10 TF BARPLOTS PER CELL TYPE
# ==============================================================================

for (ct in unique(tf_all$celltype)) {

  top_ct <- tf_all %>%
    dplyr::filter(
      celltype == ct,
      is.finite(estimate)
    ) %>%
    dplyr::arrange(
      dplyr::desc(abs(estimate))
    ) %>%
    dplyr::slice_head(
      n = TOP_TF_BARPLOT
    ) %>%
    dplyr::arrange(estimate) %>%
    dplyr::mutate(
      source = as.character(source),
      direction = ifelse(
        estimate >= 0,
        "Higher in DMD",
        "Lower in DMD"
      ),
      plot_estimate = pmax(
        -4,
        pmin(4, as.numeric(estimate))
      )
    ) %>%
    dplyr::mutate(
      source = stats::reorder(
        source,
        plot_estimate
      )
    )

  if (nrow(top_ct) == 0) {
    next
  }

  p_bar <- ggplot2::ggplot(
    top_ct,
    ggplot2::aes(
      x = source,
      y = plot_estimate,
      fill = direction
    )
  ) +
    ggplot2::geom_col(
      width = 0.75
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "Higher in DMD" = "#B2182B",
        "Lower in DMD" = "#2166AC"
      ),
      name = NULL
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.4
    ) +
    ggplot2::labs(
      x = NULL,
      y = "TF activity (ULM score; capped at +/-4)",
      title = paste0("Top ", TOP_TF_BARPLOT, " TF activities - ", ct),
      subtitle = "DMD vs CTRL | red = higher in DMD, blue = lower in DMD"
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  safe_ct <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    ct
  )

  ggplot2::ggsave(
    file.path(
      OUTDIR,
      "figures",
      paste0(
        "TF_activity_top",
        TOP_TF_BARPLOT,
        "_",
        safe_ct,
        ".pdf"
      )
    ),
    p_bar,
    width = 7,
    height = 5
  )

  ggplot2::ggsave(
    file.path(
      OUTDIR,
      "figures",
      paste0(
        "TF_activity_top",
        TOP_TF_BARPLOT,
        "_",
        safe_ct,
        ".png"
      )
    ),
    p_bar,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ==============================================================================
# 10. SESSION INFO
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
