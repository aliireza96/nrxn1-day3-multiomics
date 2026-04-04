#!/usr/bin/env Rscript

# 04_h3k27me3.R
# Rebuild H3K27me3 differential-region and promoter summaries for the manuscript.
#
# Required upstream inputs:
#   - metadata/processed_input_manifest.tsv
#   - H3K27me3 differential-region and promoter-region tables documented there
#
# Main outputs:
#   - outputs/chip_seq/<comparison_id>/tables/
#   - outputs/chip_seq/<comparison_id>/plots/
#
# Run after:
#   - 00_validate_inputs.R

suppressPackageStartupMessages({
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
  }
})

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    sentinel_file <- file.path(current, ".nrxn1_project_root")
    if (file.exists(sentinel_file)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate project root from current working directory.")
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

collapse_chip_annotation <- function(x) {
  x <- as.character(x)
  out <- rep("Other", length(x))
  out[grepl("Promoter", x, ignore.case = TRUE)] <- "Promoter"
  out[grepl("Distal", x, ignore.case = TRUE)] <- "Distal"
  out[grepl("Intergenic", x, ignore.case = TRUE)] <- "Intergenic"
  out[grepl("Intron", x, ignore.case = TRUE)] <- "Intron"
  out[grepl("Exon", x, ignore.case = TRUE)] <- "Exon"
  out[grepl("UTR", x, ignore.case = TRUE) | grepl("^3'|^5'", x)] <- "UTR"
  out
}

load_source_map <- function(project_root) {
  read_tsv_base(file.path(project_root, "metadata", "processed_input_manifest.tsv"))
}

lookup_source <- function(source_map, comparison_id, source_type, project_root) {
  hit <- source_map[source_map$comparison_id == comparison_id & source_map$source_type == source_type, ]
  if (!nrow(hit)) {
    stop("Missing source mapping for comparison_id=", comparison_id, ", source_type=", source_type)
  }
  file.path(project_root, hit$source_path[[1]])
}

standardize_chip_tables <- function(diff_regions_path, promoter_regions_path, comparison_id) {
  diff_regions <- read_tsv_base(diff_regions_path)
  promoter_regions <- read_tsv_base(promoter_regions_path)

  diff_regions$annotation_group <- if ("annotation2" %in% names(diff_regions)) {
    collapse_chip_annotation(diff_regions$annotation2)
  } else {
    collapse_chip_annotation(diff_regions$annotation)
  }
  diff_regions$comparison_id <- comparison_id

  promoter_regions$annotation_group <- "Promoter"
  promoter_regions$comparison_id <- comparison_id

  list(full = diff_regions, promoter = promoter_regions)
}

summarize_chip <- function(full_table, promoter_table) {
  data.frame(
    metric = c(
      "total_significant_regions",
      "positive_fold_regions",
      "negative_fold_regions",
      "promoter_rows",
      "promoter_unique_genes"
    ),
    value = c(
      nrow(full_table),
      sum(as.numeric(full_table$Fold) > 0, na.rm = TRUE),
      sum(as.numeric(full_table$Fold) < 0, na.rm = TRUE),
      nrow(promoter_table),
      length(unique(promoter_table$external_gene_name))
    ),
    stringsAsFactors = FALSE
  )
}

plot_chip_volcano <- function(full_table, output_path, title_text) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }

  plot_df <- full_table
  plot_df$Fold <- as.numeric(plot_df$Fold)
  plot_df$FDR <- as.numeric(plot_df$FDR)
  plot_df$neg_log10_fdr <- -log10(pmax(plot_df$FDR, .Machine$double.xmin))

  colors <- c(
    Promoter = "#7f3c8d",
    Distal = "#11a579",
    Intergenic = "#3969ac",
    Intron = "#f2b701",
    Exon = "#e73f74",
    UTR = "#80ba5a",
    Other = "#9e9e9e"
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Fold, y = neg_log10_fdr, color = annotation_group)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 1.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = title_text,
      x = "H3K27me3 log2 fold-change",
      y = "-log10 FDR",
      color = "Annotation"
    ) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(output_path, plot = p, width = 7, height = 5.5)
}

plot_chip_direction_bar <- function(promoter_table, output_path, title_text) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }

  promoter_table$Fold <- as.numeric(promoter_table$Fold)
  counts <- data.frame(
    direction = c("Higher in patient", "Higher in control"),
    n = c(
      sum(promoter_table$Fold > 0, na.rm = TRUE),
      sum(promoter_table$Fold < 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(counts, ggplot2::aes(x = direction, y = n, fill = direction)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::scale_fill_manual(values = c("Higher in patient" = "#c83e4d", "Higher in control" = "#2f6bff")) +
    ggplot2::labs(title = title_text, x = NULL, y = "Promoter peaks") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(output_path, plot = p, width = 5.5, height = 4.5)
}

write_curated_tables <- function(promoter_table, output_dir) {
  curated_genes <- c("SMAD7", "CDH1", "CDH2", "EZH2", "EED", "SOX9", "OLIG2", "FRAT1", "CNTN4")
  hits <- promoter_table[promoter_table$external_gene_name %in% curated_genes, ]
  write_tsv_base(hits, file.path(output_dir, "curated_promoter_genes.tsv"))

  smad7_hits <- promoter_table[promoter_table$external_gene_name == "SMAD7", ]
  write_tsv_base(smad7_hits, file.path(output_dir, "smad7_promoter_regions.tsv"))
}

project_root <- find_project_root()
source_map <- load_source_map(project_root)
comparisons <- read_tsv_base(file.path(project_root, "metadata", "comparisons.tsv"))

chip_comparisons <- comparisons[
  comparisons$modality == "ChIP_seq" &
    comparisons$assay == "H3K27me3" &
    comparisons$status == "active",
]

for (i in seq_len(nrow(chip_comparisons))) {
  comparison_id <- chip_comparisons$comparison_id[[i]]
  output_dir <- file.path(project_root, "outputs", "chip_seq", comparison_id)
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  log_dir <- file.path(output_dir, "logs")
  ensure_dir(table_dir)
  ensure_dir(figure_dir)
  ensure_dir(log_dir)

  diff_regions_path <- lookup_source(source_map, comparison_id, "diff_regions", project_root)
  promoter_regions_path <- lookup_source(source_map, comparison_id, "promoter_regions", project_root)
  standardized <- standardize_chip_tables(diff_regions_path, promoter_regions_path, comparison_id)

  write_tsv_base(standardized$full, file.path(table_dir, "differential_regions_standardized.tsv"))
  write_tsv_base(standardized$promoter, file.path(table_dir, "promoter_regions_standardized.tsv"))
  write_tsv_base(summarize_chip(standardized$full, standardized$promoter), file.path(table_dir, "summary_metrics.tsv"))
  write_curated_tables(standardized$promoter, table_dir)

  plot_chip_volcano(
    standardized$full,
    file.path(figure_dir, paste0(comparison_id, "_volcano.pdf")),
    paste(comparison_id, "H3K27me3 differential regions")
  )
  plot_chip_direction_bar(
    standardized$promoter,
    file.path(figure_dir, paste0(comparison_id, "_promoter_direction_bar.pdf")),
    paste(comparison_id, "promoter direction summary")
  )

  capture.output(
    list(
      comparison_id = comparison_id,
      source_files = source_map[source_map$comparison_id == comparison_id, c("source_type", "source_path")],
      generated_tables = list.files(table_dir),
      generated_figures = list.files(figure_dir)
    ),
    file = file.path(log_dir, "rebuild_log.txt")
  )
}
