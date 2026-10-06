#!/usr/bin/env Rscript

# DMD snRNA-seq: donor-aware TF activity analysis
# 3 DMD vs 2 CTRL
#
# Workflow:
#   Seurat RDS
#     -> donor x cell-type pseudobulk
#     -> DESeq2 (DMD vs CTRL)
#     -> CollecTRI TF regulons
#     -> decoupleR ULM using the DESeq2 Wald statistic
#     -> per-cell-type TF activity tables + global heatmap
#
# IMPORTANT:
# Biological replicates are donors, not nuclei.

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
# CONFIG
# ==============================================================================

RDS_PATH <- "data/DMD_snRNAseq.rds"

ASSAY <- "RNA"

# Metadata columns in the Seurat object
DONOR_COL <- "sample_id"
CONDITION_COL <- "condition"
CELLTYPE_COL <- "celltype"

CTRL_LABEL <- "CTRL"
DMD_LABEL <- "DMD"

# Minimum number of nuclei from one donor in one cell type
MIN_NUCLEI_PER_DONOR_CELLTYPE <- 20

# Require at least this many donors from EACH condition for a cell type
MIN_DONORS_PER_CONDITION <- 2

# Gene filtering before DESeq2
MIN_COUNT <- 10
MIN_SAMPLES_WITH_COUNT <- 2

# CollecTRI / decoupleR
MIN_REGULON_SIZE <- 10

OUTDIR <- "results_tf_activity"

# ==============================================================================
# HELPERS
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
    paste(missing_meta, collapse = ", ")
  )
}

if (!ASSAY %in% Assays(obj)) {
  stop("Assay '", ASSAY, "' not found in the Seurat object.")
}

DefaultAssay(obj) <- ASSAY

meta <- obj@meta.data %>%
  mutate(
    .cell = rownames(obj@meta.data),
    donor = as.character(.data[[DONOR_COL]]),
    condition = as.character(.data[[CONDITION_COL]]),
    celltype = as.character(.data[[CELLTYPE_COL]])
  )

if (!all(c(CTRL_LABEL, DMD_LABEL) %in% unique(meta$condition))) {
  stop(
    "CTRL_LABEL/DMD_LABEL do not match CONDITION_COL values. Found: ",
    paste(sort(unique(meta$condition)), collapse = ", ")
  )
}

# ==============================================================================
# 1. QC: nuclei per donor x cell type
# ==============================================================================

qc <- meta %>%
  count(celltype, donor, condition, name = "n_nuclei") %>%
  arrange(celltype, condition, donor)

write.table(
  qc,
  file.path(OUTDIR, "nuclei_per_donor_celltype.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

eligible <- qc %>%
  filter(n_nuclei >= MIN_NUCLEI_PER_DONOR_CELLTYPE) %>%
  count(celltype, condition, name = "n_donors") %>%
  filter(condition %in% c(CTRL_LABEL, DMD_LABEL)) %>%
  pivot_wider(
    names_from = condition,
    values_from = n_donors,
    values_fill = 0
  )

if (!CTRL_LABEL %in% colnames(eligible)) eligible[[CTRL_LABEL]] <- 0
if (!DMD_LABEL %in% colnames(eligible)) eligible[[DMD_LABEL]] <- 0

eligible_celltypes <- eligible %>%
  filter(
    .data[[CTRL_LABEL]] >= MIN_DONORS_PER_CONDITION,
    .data[[DMD_LABEL]] >= MIN_DONORS_PER_CONDITION
  ) %>%
  pull(celltype)

if (length(eligible_celltypes) == 0) {
  stop("No cell type passes the donor/nuclei thresholds.")
}

message(
  "Eligible cell types (", length(eligible_celltypes), "): ",
  paste(eligible_celltypes, collapse = ", ")
)

# ==============================================================================
# 2. CollecTRI network
# ==============================================================================

message("Loading CollecTRI...")
collectri <- decoupleR::get_collectri(
  organism = "human",
  split_complexes = FALSE
)

# Standardize expected column names across decoupleR versions
if ("weight" %in% colnames(collectri) && !"mor" %in% colnames(collectri)) {
  collectri <- collectri %>% rename(mor = weight)
}

required_net_cols <- c("source", "target", "mor")
missing_net_cols <- setdiff(required_net_cols, colnames(collectri))

if (length(missing_net_cols) > 0) {
  stop(
    "Unexpected CollecTRI format. Missing columns: ",
    paste(missing_net_cols, collapse = ", ")
  )
}

# ==============================================================================
# 3. Counts
# ==============================================================================

counts <- GetAssayData(obj, assay = ASSAY, layer = "counts")

if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

# ==============================================================================
# 4. Per-cell-type pseudobulk -> DESeq2 -> decoupleR ULM
# ==============================================================================

all_de <- list()
all_tf <- list()

for (ct in eligible_celltypes) {

  message("\n============================================================")
  message("Cell type: ", ct)
  message("============================================================")

  ct_meta <- meta %>%
    filter(
      celltype == ct,
      condition %in% c(CTRL_LABEL, DMD_LABEL)
    )

  # Exclude donor/cell-type combinations with too few nuclei
  valid_donors <- ct_meta %>%
    count(donor, condition, name = "n_nuclei") %>%
    filter(n_nuclei >= MIN_NUCLEI_PER_DONOR_CELLTYPE)

  ct_meta <- ct_meta %>%
    semi_join(valid_donors, by = c("donor", "condition"))

  if (nrow(ct_meta) == 0) next

  donor_levels <- unique(ct_meta$donor)
  donor_factor <- factor(ct_meta$donor, levels = donor_levels)

  # Sparse aggregation matrix: nuclei -> donors
  design <- sparse.model.matrix(~ 0 + donor_factor)
  colnames(design) <- donor_levels

  ct_counts <- counts[, ct_meta$.cell, drop = FALSE]
  pseudobulk <- ct_counts %*% design

  # Donor metadata
  donor_meta <- ct_meta %>%
    distinct(donor, condition) %>%
    slice(match(colnames(pseudobulk), donor))

  rownames(donor_meta) <- donor_meta$donor
  donor_meta <- donor_meta[colnames(pseudobulk), , drop = FALSE]

  donor_meta$condition <- factor(
    donor_meta$condition,
    levels = c(CTRL_LABEL, DMD_LABEL)
  )

  # Gene filter
  keep_genes <- rowSums(pseudobulk >= MIN_COUNT) >= MIN_SAMPLES_WITH_COUNT
  pseudobulk <- pseudobulk[keep_genes, , drop = FALSE]

  message(
    "Donors: ",
    sum(donor_meta$condition == CTRL_LABEL), " CTRL / ",
    sum(donor_meta$condition == DMD_LABEL), " DMD"
  )
  message("Genes retained: ", nrow(pseudobulk))

  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(pseudobulk)),
    colData = donor_meta,
    design = ~ condition
  )

  dds <- DESeq(dds, quiet = TRUE)

  contrast_name <- paste0("condition_", DMD_LABEL, "_vs_", CTRL_LABEL)

  if (!contrast_name %in% resultsNames(dds)) {
    stop(
      "DESeq2 contrast not found for ", ct,
      ". Available coefficients: ",
      paste(resultsNames(dds), collapse = ", ")
    )
  }

  res <- results(dds, name = contrast_name)

  de <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      celltype = ct,
      padj = ifelse(is.na(padj), 1, padj)
    ) %>%
    select(
      celltype,
      gene,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    ) %>%
    arrange(padj, desc(abs(stat)))

  safe_ct <- gsub("[^A-Za-z0-9._-]+", "_", ct)

  write.table(
    de,
    file.path(
      OUTDIR,
      "DESeq2",
      paste0(safe_ct, "_DMD_vs_CTRL.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  all_de[[ct]] <- de

  # --------------------------------------------------------------------------
  # TF activity:
  # Use DESeq2 Wald statistic rather than log2FC as the enrichment input.
  # Rows = genes, column = DMD_vs_CTRL.
  # --------------------------------------------------------------------------

  stat_tbl <- de %>%
    filter(is.finite(stat), !is.na(gene)) %>%
    distinct(gene, .keep_all = TRUE)

  stat_mat <- matrix(
    stat_tbl$stat,
    ncol = 1,
    dimnames = list(stat_tbl$gene, "DMD_vs_CTRL")
  )

  tf <- decoupleR::run_ulm(
    mat = stat_mat,
    network = collectri,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = MIN_REGULON_SIZE
  )

  # Support common output names used across decoupleR releases
  if ("condition" %in% colnames(tf) && !"sample" %in% colnames(tf)) {
    tf <- tf %>% rename(sample = condition)
  }

  if ("score" %in% colnames(tf) && !"estimate" %in% colnames(tf)) {
    tf <- tf %>% rename(estimate = score)
  }

  if ("p_value" %in% colnames(tf) && !"pvalue" %in% colnames(tf)) {
    tf <- tf %>% rename(pvalue = p_value)
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
    mutate(
      celltype = ct,
      padj = ifelse(
        is.na(pvalue),
        NA_real_,
        p.adjust(pvalue, method = "BH")
      )
    ) %>%
    select(
      celltype,
      source,
      estimate,
      pvalue,
      padj,
      everything()
    ) %>%
    arrange(desc(abs(estimate)))

  write.table(
    tf,
    file.path(
      OUTDIR,
      "TF_activity",
      paste0(safe_ct, "_CollecTRI_ULM.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  all_tf[[ct]] <- tf

  rm(dds, res, pseudobulk, ct_counts, design)
  gc()
}

# ==============================================================================
# 5. Combined outputs
# ==============================================================================

de_all <- bind_rows(all_de)
tf_all <- bind_rows(all_tf)

write.table(
  de_all,
  file.path(OUTDIR, "DESeq2", "all_celltypes_DMD_vs_CTRL.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tf_all,
  file.path(OUTDIR, "TF_activity", "all_celltypes_CollecTRI_ULM.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================================
# 6. Global TF heatmap
# ==============================================================================

if (nrow(tf_all) > 0) {

  # Keep TFs with strongest activity in at least one cell type.
  # This is only a plotting filter; full results remain exported.
  top_tfs <- tf_all %>%
    group_by(source) %>%
    summarise(max_abs_activity = max(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_activity)) %>%
    slice_head(n = min(30, n())) %>%
    pull(source)

  heat_df <- tf_all %>%
    filter(source %in% top_tfs) %>%
    mutate(
      source = factor(source, levels = rev(top_tfs)),
      celltype = factor(celltype, levels = eligible_celltypes)
    )

  p <- ggplot(
    heat_df,
    aes(
      x = celltype,
      y = source,
      fill = estimate
    )
  ) +
    geom_tile() +
    scale_fill_gradient2(
      midpoint = 0,
      name = "TF activity\nULM score"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "Differential TF activity in DMD",
      subtitle = "CollecTRI + decoupleR ULM on DESeq2 Wald statistics"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(
    file.path(OUTDIR, "figures", "TF_activity_heatmap_top30.pdf"),
    p,
    width = max(7, 0.6 * length(eligible_celltypes) + 4),
    height = 10
  )

  ggsave(
    file.path(OUTDIR, "figures", "TF_activity_heatmap_top30.png"),
    p,
    width = max(7, 0.6 * length(eligible_celltypes) + 4),
    height = 10,
    dpi = 300
  )
}

# ==============================================================================
# 7. Session information
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTDIR, "sessionInfo.txt")
)

message("\nDone.")
message("Results: ", normalizePath(OUTDIR))
