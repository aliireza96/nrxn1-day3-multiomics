#!/usr/bin/env Rscript

# 03_atac_accessibility.R
# Rebuild ATAC differential-accessibility, promoter annotation, and TF-footprinting
# summary tables used in the manuscript.
#
# Required upstream inputs:
#   - metadata/processed_input_manifest.tsv
#   - ATAC differential-accessibility and annotation tables documented there
#
# Main outputs:
#   - outputs/atac_seq/<comparison_id>/tables/
#   - outputs/atac_seq/<comparison_id>/plots/
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

read_table_auto <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }
}

write_tsv_base <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

collapse_atac_annotation <- function(x) {
  x <- as.character(x)
  out <- rep("Other", length(x))
  out[grepl("^promoter", x, ignore.case = TRUE) | grepl("^Promoter", x)] <- "Promoter"
  out[grepl("Intergenic", x, ignore.case = TRUE)] <- "Intergenic"
  out[grepl("Distal", x, ignore.case = TRUE)] <- "Distal"
  out[grepl("Intron", x, ignore.case = TRUE)] <- "Intron"
  out[grepl("Exon", x, ignore.case = TRUE)] <- "Exon"
  out[grepl("UTR", x, ignore.case = TRUE)] <- "UTR"
  out
}

make_peak_id <- function(seqnames, start, end) {
  paste("ID", seqnames, start, end, sep = "_")
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

standardize_main_atac <- function(diff_peaks_path, annotations_path) {
  diff_peaks <- read_tsv_base(diff_peaks_path)
  annotations <- read_tsv_base(annotations_path)

  diff_peaks$PeakID <- make_peak_id(diff_peaks$seqnames, diff_peaks$start, diff_peaks$end)
  merged <- merge(diff_peaks, annotations, by.x = "PeakID", by.y = "PeakID (cmd=annotatePeaks.pl P_vs_ctrl_sig.bed hg38)", all.x = TRUE, sort = FALSE)
  merged$annotation_group <- collapse_atac_annotation(merged$Annotation)
  merged$gene_name <- merged$`Gene Name`
  merged$comparison_id <- "atac_main_local_controls"
  merged
}

standardize_supp_atac <- function(diff_peaks_path, annotations_path, promoter_annotations_path) {
  diff_peaks <- read_tsv_base(diff_peaks_path)
  annotations <- read_tsv_base(annotations_path)
  promoter_annotations <- read_tsv_base(promoter_annotations_path)

  diff_peaks$PeakID <- make_peak_id(diff_peaks$seqnames, diff_peaks$start, diff_peaks$end)
  annotations$PeakID <- if (!"PeakID" %in% names(annotations)) {
    make_peak_id(annotations$seqnames, annotations$start, annotations$end)
  } else {
    annotations$PeakID
  }

  merged <- merge(diff_peaks, annotations, by = "PeakID", all.x = TRUE, sort = FALSE)
  for (base_col in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
    x_col <- paste0(base_col, ".x")
    if (!base_col %in% names(merged) && x_col %in% names(merged)) {
      merged[[base_col]] <- merged[[x_col]]
    }
  }
  merged$annotation_group <- collapse_atac_annotation(if ("annotation2" %in% names(merged)) merged$annotation2 else merged$annotation)
  merged$comparison_id <- "atac_supp_female_only"

  promoter_annotations$annotation_group <- "Promoter"
  promoter_annotations$comparison_id <- "atac_supp_female_only"

  list(full = merged, promoter = promoter_annotations)
}

summarize_atac <- function(full_table, promoter_table) {
  promoter_gene_field <- intersect(c("external_gene_name", "gene_name", "Gene Name"), names(promoter_table))
  promoter_gene_count <- if (length(promoter_gene_field)) {
    length(unique(promoter_table[[promoter_gene_field[[1]]]]))
  } else {
    NA_integer_
  }

  data.frame(
    metric = c(
      "total_significant_dars",
      "positive_log2fc_dars",
      "negative_log2fc_dars",
      "promoter_rows",
      "promoter_unique_genes"
    ),
    value = c(
      nrow(full_table),
      sum(full_table$log2FoldChange > 0, na.rm = TRUE),
      sum(full_table$log2FoldChange < 0, na.rm = TRUE),
      nrow(promoter_table),
      promoter_gene_count
    ),
    stringsAsFactors = FALSE
  )
}

prepare_motif_table <- function(path, comparison_id) {
  motif <- read_table_auto(path)

  change_col <- grep("_change$", names(motif), value = TRUE)[1]
  pvalue_col <- grep("_pvalue$", names(motif), value = TRUE)[1]
  if (is.na(change_col) || is.na(pvalue_col)) {
    stop("Could not identify change/pvalue columns in motif results: ", path)
  }

  motif$binding_change <- as.numeric(motif[[change_col]])
  motif$p_value <- as.numeric(motif[[pvalue_col]])
  motif$neg_log10_p <- -log10(pmax(motif$p_value, .Machine$double.xmin))
  motif$direction <- ifelse(motif$binding_change > 0, "More accessible in patient", "Less accessible in patient")
  motif$comparison_id <- comparison_id
  motif
}

select_motif_labels <- function(motif_table) {
  priority_labels <- c(
    "KLF1", "KLF4", "KLF5", "KLF14", "KLF2", "KLF3",
    "POU3F4", "Pou5f1::Sox2", "SOX9", "SOX13",
    "TEAD1", "TEAD2", "TEAD3", "TEAD4", "OLIG2", "MEIS1"
  )

  present <- motif_table$name %in% priority_labels
  fallback <- order(abs(motif_table$binding_change), motif_table$p_value, decreasing = TRUE)
  fallback_names <- unique(motif_table$name[fallback])[seq_len(min(8, length(unique(motif_table$name[fallback]))))]
  unique(c(priority_labels[priority_labels %in% motif_table$name], fallback_names))
}

plot_atac_volcano <- function(df, output_path, title_text) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }

  colors <- c(
    Promoter = "#7f3c8d",
    Distal = "#11a579",
    Intergenic = "#3969ac",
    Intron = "#f2b701",
    Exon = "#e73f74",
    UTR = "#80ba5a",
    Other = "#9e9e9e"
  )

  plot_df <- df
  p_col <- intersect(c("padj", "FDR", "pvalue", "PValue"), names(plot_df))[1]
  if (is.na(p_col)) {
    stop("Could not identify an adjusted p-value column for ATAC volcano plotting.")
  }
  plot_df[[p_col]] <- as.numeric(plot_df[[p_col]])
  plot_df$neg_log10_padj <- -log10(pmax(plot_df[[p_col]], .Machine$double.xmin))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = annotation_group)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = title_text,
      x = "ATAC log2 fold-change",
      y = "-log10 adjusted p-value",
      color = "Annotation"
    ) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(output_path, plot = p, width = 7, height = 5.5)
}

plot_motif_volcano <- function(motif_table, output_path, title_text) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(invisible(NULL))
  }

  labels <- select_motif_labels(motif_table)
  motif_table$label <- ifelse(motif_table$name %in% labels, motif_table$name, "")

  p <- ggplot2::ggplot(
    motif_table,
    ggplot2::aes(x = binding_change, y = neg_log10_p, color = direction)
  ) +
    ggplot2::geom_point(alpha = 0.55, size = 1.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = c(
      "More accessible in patient" = "#2f6bff",
      "Less accessible in patient" = "#d1495b"
    )) +
    ggplot2::labs(
      title = title_text,
      x = "Differential binding score",
      y = "-log10 p-value",
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 12)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = motif_table[motif_table$label != "", ],
      ggplot2::aes(label = label),
      size = 3,
      max.overlaps = 50
    )
  }

  ggplot2::ggsave(output_path, plot = p, width = 8, height = 5.5)
}

write_curated_tables <- function(promoter_table, motif_table, output_dir) {
  promoter_gene_field <- intersect(c("external_gene_name", "gene_name", "Gene Name"), names(promoter_table))
  if (length(promoter_gene_field)) {
    genes <- promoter_gene_field[[1]]
    curated_genes <- c("SMAD7", "CDH1", "CDH2", "SOX9", "OLIG2", "FZD10", "CNTN4", "FRAT1")
    gene_hits <- promoter_table[promoter_table[[genes]] %in% curated_genes, ]
    write_tsv_base(gene_hits, file.path(output_dir, "curated_promoter_genes.tsv"))
  }

  motif_hits <- motif_table[motif_table$name %in% select_motif_labels(motif_table), ]
  write_tsv_base(motif_hits, file.path(output_dir, "curated_tf_motifs.tsv"))
}

project_root <- find_project_root()
source_map <- load_source_map(project_root)
comparisons <- read_tsv_base(file.path(project_root, "metadata", "comparisons.tsv"))

atac_comparisons <- comparisons[comparisons$modality == "ATAC_seq" & comparisons$status == "active", ]

for (i in seq_len(nrow(atac_comparisons))) {
  comparison_id <- atac_comparisons$comparison_id[[i]]
  output_dir <- file.path(project_root, "outputs", "atac_seq", comparison_id)
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  log_dir <- file.path(output_dir, "logs")
  ensure_dir(table_dir)
  ensure_dir(figure_dir)
  ensure_dir(log_dir)

  motif_path <- lookup_source(source_map, comparison_id, "motif_results", project_root)
  motif_table <- prepare_motif_table(motif_path, comparison_id)
  write_tsv_base(motif_table, file.path(table_dir, "motif_results_standardized.tsv"))
  plot_motif_volcano(motif_table, file.path(figure_dir, paste0(comparison_id, "_tobias_volcano.pdf")), paste(comparison_id, "TF footprinting"))

  if (comparison_id == "atac_main_local_controls") {
    diff_path <- lookup_source(source_map, comparison_id, "diff_peaks", project_root)
    annot_path <- lookup_source(source_map, comparison_id, "annotations", project_root)
    full_table <- standardize_main_atac(diff_path, annot_path)
    promoter_table <- full_table[full_table$annotation_group == "Promoter" & !is.na(full_table$gene_name), ]
  } else if (comparison_id == "atac_supp_female_only") {
    diff_path <- lookup_source(source_map, comparison_id, "diff_peaks", project_root)
    annot_path <- lookup_source(source_map, comparison_id, "annotations", project_root)
    promoter_path <- lookup_source(source_map, comparison_id, "promoter_annotations", project_root)
    standardized <- standardize_supp_atac(diff_path, annot_path, promoter_path)
    full_table <- standardized$full
    promoter_table <- standardized$promoter
  } else {
    next
  }

  write_tsv_base(full_table, file.path(table_dir, "differential_accessibility_standardized.tsv"))
  write_tsv_base(promoter_table, file.path(table_dir, "promoter_dars_standardized.tsv"))
  write_tsv_base(summarize_atac(full_table, promoter_table), file.path(table_dir, "summary_metrics.tsv"))
  write_curated_tables(promoter_table, motif_table, table_dir)

  plot_atac_volcano(full_table, file.path(figure_dir, paste0(comparison_id, "_volcano.pdf")), paste(comparison_id, "differential accessibility"))

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
