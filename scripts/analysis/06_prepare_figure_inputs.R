#!/usr/bin/env Rscript

# 06_prepare_figure_inputs.R
# Snapshot figure-input readiness after the assay-level analyses have run.
#
# Required upstream inputs:
#   - metadata/figure_manifest.tsv
#   - assay outputs produced under outputs/
#
# Main outputs:
#   - outputs/figures/tables/figure_input_status.tsv
#   - outputs/figures/tables/figure_manifest_snapshot.tsv
#
# Run after:
#   - 01_rna_expression.R
#   - 02_splicing.R
#   - 03_atac_accessibility.R
#   - 04_h3k27me3.R
#   - 05_integration.R

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

project_root <- find_project_root()
figure_manifest <- read_tsv_base(file.path(project_root, "metadata", "figure_manifest.tsv"))

output_dir <- file.path(project_root, "outputs", "figures")
table_dir <- file.path(output_dir, "tables")
log_dir <- file.path(output_dir, "logs")
ensure_dir(table_dir)
ensure_dir(log_dir)

expected_inputs <- data.frame(
  figure_id = c(
    "Figure_1", "Figure_2", "Figure_2", "Figure_3", "Figure_4", "Figure_5", "Figure_5"
  ),
  panel = c("C", "A", "D", "D", "B", "B", "C"),
  expected_path = c(
    "outputs/rna_seq/rna_main_all_controls/figures/rna_main_all_controls_pca.pdf",
    "outputs/rna_seq/rna_main_all_controls/tables/rna_main_all_controls_gene_results_sig_lfc.tsv",
    "outputs/splicing/tables/summary_metrics.tsv",
    "outputs/chip_seq/chip_main_all_controls/tables/smad7_promoter_regions.tsv",
    "outputs/atac_seq/atac_main_local_controls/tables/curated_tf_motifs.tsv",
    "outputs/integration/tables/triple_overlap_main_heatmap_matrix.tsv",
    "outputs/integration/tables/modality_correlations.tsv"
  ),
  stringsAsFactors = FALSE
)

expected_inputs$exists <- file.exists(file.path(project_root, expected_inputs$expected_path))
write_tsv_base(expected_inputs, file.path(table_dir, "figure_input_status.tsv"))
write_tsv_base(figure_manifest, file.path(table_dir, "figure_manifest_snapshot.tsv"))

capture.output(
  list(
    missing_inputs = expected_inputs[!expected_inputs$exists, ],
    manifest_rows = nrow(figure_manifest)
  ),
  file = file.path(log_dir, "figure_input_log.txt")
)
