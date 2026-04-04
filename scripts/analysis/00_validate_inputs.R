#!/usr/bin/env Rscript

# 00_validate_inputs.R
# Validate the public manuscript repository before running the workflow.
#
# Purpose:
#   - check that metadata files are present and internally consistent
#   - confirm that documented assay-level input paths exist
#   - write a validation report to outputs/logs/
#
# Inputs:
#   - metadata/sample_sheet.tsv
#   - metadata/comparisons.tsv
#   - metadata/processed_input_manifest.tsv
#   - metadata/figure_manifest.tsv
#
# Outputs:
#   - outputs/logs/input_validation_report.txt

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

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

split_csv_field <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

project_root <- find_project_root()
metadata_dir <- file.path(project_root, "metadata")
log_dir <- file.path(project_root, "outputs", "logs")
ensure_dir(log_dir)

sample_sheet <- read_tsv_base(file.path(metadata_dir, "sample_sheet.tsv"))
comparisons <- read_tsv_base(file.path(metadata_dir, "comparisons.tsv"))
processed_inputs <- read_tsv_base(file.path(metadata_dir, "processed_input_manifest.tsv"))
figure_manifest <- read_tsv_base(file.path(metadata_dir, "figure_manifest.tsv"))

required_sample_cols <- c(
  "sample_id", "modality", "assay", "line_id", "line_type", "clone_id",
  "sex", "density", "genotype_group", "include_manuscript"
)
required_comparison_cols <- c(
  "comparison_id", "modality", "assay", "sample_ids", "design_formula",
  "contrast_column", "numerator_level", "denominator_level", "status"
)

report <- list()
report$missing_sample_columns <- setdiff(required_sample_cols, names(sample_sheet))
report$missing_comparison_columns <- setdiff(required_comparison_cols, names(comparisons))
report$duplicate_sample_modality_pairs <- sample_sheet[duplicated(sample_sheet[, c("sample_id", "modality", "assay")]), c("sample_id", "modality", "assay")]
report$missing_input_counts <- sample_sheet$input_counts_path[
  nzchar(sample_sheet$input_counts_path) &
    !file.exists(file.path(project_root, sample_sheet$input_counts_path))
]
report$missing_peak_paths <- sample_sheet$peak_path[
  nzchar(sample_sheet$peak_path) &
    !file.exists(file.path(project_root, sample_sheet$peak_path))
]
report$missing_processed_inputs <- processed_inputs$source_path[
  !file.exists(file.path(project_root, processed_inputs$source_path))
]

comparison_issues <- list()
for (i in seq_len(nrow(comparisons))) {
  comp <- comparisons[i, ]
  sample_ids <- split_csv_field(comp$sample_ids)
  matching <- sample_sheet[
    sample_sheet$sample_id %in% sample_ids &
      sample_sheet$modality == comp$modality &
      sample_sheet$assay == comp$assay,
  ]
  missing_ids <- setdiff(sample_ids, matching$sample_id)
  if (length(missing_ids)) {
    comparison_issues[[length(comparison_issues) + 1]] <- data.frame(
      comparison_id = comp$comparison_id,
      issue = "missing_samples",
      details = paste(missing_ids, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
}
report$comparison_issues <- if (length(comparison_issues)) do.call(rbind, comparison_issues) else data.frame()
report$figure_manifest_rows <- nrow(figure_manifest)

capture.output(report, file = file.path(log_dir, "input_validation_report.txt"))
