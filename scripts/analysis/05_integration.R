#!/usr/bin/env Rscript

# 05_integration.R
# Rebuild the RNA-ATAC, RNA-H3K27me3, and triple-overlap integration tables used
# across the main manuscript figures and supplementary tables.
#
# Required upstream inputs:
#   - outputs/rna_seq/rna_main_all_controls/tables/
#   - outputs/atac_seq/atac_main_local_controls/tables/
#   - outputs/chip_seq/chip_main_all_controls/tables/
#   - config/curated_features.tsv
#
# Main outputs:
#   - outputs/integration/tables/
#
# Run after:
#   - 01_rna_expression.R
#   - 03_atac_accessibility.R
#   - 04_h3k27me3.R

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, ".nrxn1_project_root"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate project root. Ensure the .nrxn1_project_root sentinel file exists.")
    }
    current <- parent
  }
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_base <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

collapse_gene_level <- function(df, gene_col, fc_col, p_col, output_fc_name, output_p_name) {
  df <- df[!is.na(df[[gene_col]]) & nzchar(df[[gene_col]]), ]
  df[[fc_col]] <- as.numeric(df[[fc_col]])
  df[[p_col]] <- as.numeric(df[[p_col]])

  df <- df[order(df[[gene_col]], df[[p_col]], -abs(df[[fc_col]]), na.last = TRUE), ]
  df <- df[!duplicated(df[[gene_col]]), c(gene_col, fc_col, p_col)]
  names(df) <- c("gene_name", output_fc_name, output_p_name)
  df
}

load_curated_features <- function(project_root) {
  cfg <- read_tsv_base(file.path(project_root, "config", "curated_features.tsv"))
  cfg[cfg$feature_type == "gene", "feature_id", drop = TRUE]
}

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required upstream result is missing: ", path)
  }
  path
}

prepare_rna_table <- function(path) {
  df <- read_tsv_base(path)
  out <- df[, c("gene_name", "log2FoldChange", "padj")]
  out$log2FoldChange <- as.numeric(out$log2FoldChange)
  out$padj <- as.numeric(out$padj)
  out <- out[!is.na(out$gene_name) & nzchar(out$gene_name), ]
  out <- out[!duplicated(out$gene_name), ]
  names(out) <- c("gene_name", "rna_log2fc", "rna_padj")
  out
}

prepare_atac_table <- function(path) {
  df <- read_tsv_base(path)
  gene_col <- intersect(c("external_gene_name", "gene_name", "Gene Name"), names(df))[1]
  fc_col <- intersect(c("log2FoldChange", "log2FoldChange_atac"), names(df))[1]
  p_col <- intersect(c("padj", "FDR"), names(df))[1]
  collapse_gene_level(df, gene_col, fc_col, p_col, "atac_log2fc", "atac_padj")
}

prepare_chip_table <- function(path) {
  df <- read_tsv_base(path)
  collapse_gene_level(df, "external_gene_name", "Fold", "FDR", "chip_log2fc", "chip_fdr")
}

correlation_table <- function(df, x_col, y_col, output_label) {
  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  ok <- is.finite(x) & is.finite(y)

  data.frame(
    comparison = output_label,
    n = sum(ok),
    pearson_r = if (sum(ok) > 1) cor(x[ok], y[ok], method = "pearson") else NA_real_,
    spearman_rho = if (sum(ok) > 1) cor(x[ok], y[ok], method = "spearman") else NA_real_,
    stringsAsFactors = FALSE
  )
}

project_root <- find_project_root()
output_dir <- file.path(project_root, "outputs", "integration")
table_dir <- file.path(output_dir, "tables")
log_dir <- file.path(output_dir, "logs")
ensure_dir(table_dir)
ensure_dir(log_dir)

curated_features <- load_curated_features(project_root)

rna_main_path <- require_file(file.path(
  project_root, "outputs", "rna_seq",
  "rna_main_all_controls", "tables", "rna_main_all_controls_gene_results_sig_lfc.tsv"
))
rna_supp_path <- require_file(file.path(
  project_root, "outputs", "rna_seq",
  "rna_supp_female_only", "tables", "rna_supp_female_only_gene_results_sig_lfc.tsv"
))
atac_main_path <- require_file(file.path(
  project_root, "outputs", "atac_seq",
  "atac_main_local_controls", "tables", "promoter_dars_standardized.tsv"
))
atac_supp_path <- require_file(file.path(
  project_root, "outputs", "atac_seq",
  "atac_supp_female_only", "tables", "promoter_dars_standardized.tsv"
))
chip_main_path <- require_file(file.path(
  project_root, "outputs", "chip_seq",
  "chip_main_all_controls", "tables", "promoter_regions_standardized.tsv"
))
chip_supp_path <- require_file(file.path(
  project_root, "outputs", "chip_seq",
  "chip_supp_female_only", "tables", "promoter_regions_standardized.tsv"
))

rna_main <- prepare_rna_table(rna_main_path)
rna_supp <- prepare_rna_table(rna_supp_path)
atac_main <- prepare_atac_table(atac_main_path)
atac_supp <- prepare_atac_table(atac_supp_path)
chip_main <- prepare_chip_table(chip_main_path)
chip_female <- prepare_chip_table(chip_supp_path)

rna_atac_main <- merge(rna_main, atac_main, by = "gene_name")
rna_atac_female <- merge(rna_supp, atac_supp, by = "gene_name")
rna_chip_main <- merge(rna_main, chip_main, by = "gene_name")
rna_chip_female <- merge(rna_supp, chip_female, by = "gene_name")
triple_main <- Reduce(function(x, y) merge(x, y, by = "gene_name"), list(rna_main, atac_main, chip_main))
triple_female <- Reduce(function(x, y) merge(x, y, by = "gene_name"), list(rna_supp, atac_supp, chip_female))

write_tsv_base(rna_atac_main, file.path(table_dir, "rna_atac_main_overlap.tsv"))
write_tsv_base(rna_atac_female, file.path(table_dir, "rna_atac_female_overlap.tsv"))
write_tsv_base(rna_chip_main, file.path(table_dir, "rna_chip_main_overlap.tsv"))
write_tsv_base(rna_chip_female, file.path(table_dir, "rna_chip_female_overlap.tsv"))
write_tsv_base(triple_main, file.path(table_dir, "triple_overlap_main.tsv"))
write_tsv_base(triple_female, file.path(table_dir, "triple_overlap_female.tsv"))

overlap_metrics <- rbind(
  data.frame(metric = "rna_atac_main_overlap_genes", value = nrow(rna_atac_main), stringsAsFactors = FALSE),
  data.frame(metric = "rna_atac_female_overlap_genes", value = nrow(rna_atac_female), stringsAsFactors = FALSE),
  data.frame(metric = "rna_chip_main_overlap_genes", value = nrow(rna_chip_main), stringsAsFactors = FALSE),
  data.frame(metric = "rna_chip_female_overlap_genes", value = nrow(rna_chip_female), stringsAsFactors = FALSE),
  data.frame(metric = "triple_overlap_main_genes", value = nrow(triple_main), stringsAsFactors = FALSE),
  data.frame(metric = "triple_overlap_female_genes", value = nrow(triple_female), stringsAsFactors = FALSE)
)
correlation_metrics <- rbind(
  correlation_table(rna_atac_main, "rna_log2fc", "atac_log2fc", "rna_atac_main"),
  correlation_table(rna_atac_female, "rna_log2fc", "atac_log2fc", "rna_atac_female"),
  correlation_table(rna_chip_main, "rna_log2fc", "chip_log2fc", "rna_chip_main"),
  correlation_table(rna_chip_female, "rna_log2fc", "chip_log2fc", "rna_chip_female"),
  correlation_table(triple_main, "atac_log2fc", "chip_log2fc", "atac_chip_main"),
  correlation_table(triple_female, "atac_log2fc", "chip_log2fc", "atac_chip_female")
)
write_tsv_base(overlap_metrics, file.path(table_dir, "summary_metrics.tsv"))
write_tsv_base(correlation_metrics, file.path(table_dir, "modality_correlations.tsv"))

triple_heatmap_matrix_main <- triple_main[, c("gene_name", "rna_log2fc", "atac_log2fc", "chip_log2fc")]
triple_heatmap_matrix_female <- triple_female[, c("gene_name", "rna_log2fc", "atac_log2fc", "chip_log2fc")]
write_tsv_base(triple_heatmap_matrix_main, file.path(table_dir, "triple_overlap_main_heatmap_matrix.tsv"))
write_tsv_base(triple_heatmap_matrix_female, file.path(table_dir, "triple_overlap_female_heatmap_matrix.tsv"))
write_tsv_base(triple_main, file.path(table_dir, "triple_overlap_fc_table.tsv"))

curated_hits_main <- triple_main[triple_main$gene_name %in% curated_features, ]
curated_hits_female <- triple_female[triple_female$gene_name %in% curated_features, ]
write_tsv_base(curated_hits_main, file.path(table_dir, "curated_triple_overlap_main.tsv"))
write_tsv_base(curated_hits_female, file.path(table_dir, "curated_triple_overlap_female.tsv"))

capture.output(
  list(
    generated_tables = list.files(table_dir),
    triple_overlap_main_n = nrow(triple_main),
    triple_overlap_n = nrow(triple_female)
  ),
  file = file.path(log_dir, "rebuild_log.txt")
)
