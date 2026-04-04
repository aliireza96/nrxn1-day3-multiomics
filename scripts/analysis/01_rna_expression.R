#!/usr/bin/env Rscript

# 01_rna_expression.R
# Rebuild RNA differential-expression outputs for the manuscript comparisons.
#
# Required upstream inputs:
#   - metadata/sample_sheet.tsv
#   - metadata/comparisons.tsv
#   - data/processed_inputs/rna/stringtie_counts/<sample_id>/{e,i,t}_data.ctab
#
# Main outputs:
#   - outputs/rna_seq/<comparison_id>/tables/
#   - outputs/rna_seq/<comparison_id>/plots/
#
# Run after:
#   - 00_validate_inputs.R

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
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

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

as_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  toupper(as.character(x)) == "TRUE"
}

write_tsv_base <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

split_csv_field <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

parse_args <- function(args) {
  parsed <- list(comparison_ids = NULL)

  for (arg in args) {
    if (startsWith(arg, "--comparison_ids=")) {
      value <- sub("^--comparison_ids=", "", arg)
      parsed$comparison_ids <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      parsed$comparison_ids <- parsed$comparison_ids[nzchar(parsed$comparison_ids)]
    } else {
      stop("Unrecognized argument: ", arg)
    }
  }

  parsed
}

make_tx2gene <- function(ctab_path) {
  ctab <- read.delim(ctab_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  unique(ctab[, c("t_name", "gene_name")])
}

safe_results <- function(dds, contrast_vector, use_ihw = TRUE) {
  if (use_ihw && requireNamespace("IHW", quietly = TRUE)) {
    return(results(dds, contrast = contrast_vector, filterFun = IHW::ihw))
  }
  results(dds, contrast = contrast_vector)
}

aggregate_transcript_results <- function(transcript_df) {
  if (!requireNamespace("aggregation", quietly = TRUE)) {
    message("Package 'aggregation' is not available; skipping Lancaster transcript aggregation.")
    return(NULL)
  }

  transcript_df <- transcript_df[!is.na(transcript_df$pvalue) & !is.na(transcript_df$gene_name), ]
  if (!nrow(transcript_df)) {
    return(NULL)
  }

  split_by_gene <- split(transcript_df, transcript_df$gene_name)
  aggregated <- lapply(split_by_gene, function(df) {
    data.frame(
      gene_name = df$gene_name[[1]],
      n_transcripts = nrow(df),
      min_log2FoldChange = min(df$log2FoldChange, na.rm = TRUE),
      max_log2FoldChange = max(df$log2FoldChange, na.rm = TRUE),
      n_log2FC_negative = sum(df$log2FoldChange < 0, na.rm = TRUE),
      n_log2FC_positive = sum(df$log2FoldChange > 0, na.rm = TRUE),
      lancaster_p = aggregation::lancaster(df$pvalue, df$baseMean),
      fisher_p = aggregation::fisher(df$pvalue),
      sidak_p = aggregation::sidak(df$pvalue),
      stringsAsFactors = FALSE
    )
  })

  aggregated <- do.call(rbind, aggregated)
  aggregated[order(aggregated$lancaster_p, aggregated$gene_name), ]
}

build_pca_plot <- function(vsd, coldata, title_text, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' is not available; skipping PCA plot.")
    return(invisible(NULL))
  }

  pca_data <- plotPCA(vsd, intgroup = c("genotype_group", "sex", "density"), returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  pca_data$sample_id <- rownames(pca_data)

  plot_obj <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x = PC1, y = PC2, color = genotype_group, shape = sex)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = title_text,
      x = paste0("PC1: ", percent_var[[1]], "% variance"),
      y = paste0("PC2: ", percent_var[[2]], "% variance"),
      color = "Group",
      shape = "Sex"
    ) +
    ggplot2::theme_bw(base_size = 12)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot_obj <- plot_obj +
      ggrepel::geom_text_repel(ggplot2::aes(label = sample_id), size = 3, max.overlaps = 30)
  }

  ggplot2::ggsave(output_path, plot = plot_obj, width = 6.5, height = 5.5)
}

run_rna_comparison <- function(comparison_row, sample_sheet, txi_tx, txi_gene, tx2gene, project_root) {
  comparison_id <- comparison_row$comparison_id
  sample_ids <- split_csv_field(comparison_row$sample_ids)
  coldata <- sample_sheet[sample_sheet$sample_id %in% sample_ids, ]
  coldata <- coldata[match(sample_ids, coldata$sample_id), ]

  if (any(is.na(coldata$sample_id))) {
    stop("Comparison ", comparison_id, " references sample IDs that are missing from the sample sheet.")
  }

  coldata$genotype_group <- factor(coldata$genotype_group, levels = c("control", "patient"))
  coldata$sex <- factor(coldata$sex, levels = c("male", "female"))
  coldata$density <- factor(coldata$density, levels = c("low", "high"))
  rownames(coldata) <- coldata$sample_id

  tx_index <- match(sample_ids, colnames(txi_tx$counts))
  gene_index <- match(sample_ids, colnames(txi_gene$counts))

  txi_tx_sub <- txi_tx
  txi_tx_sub$counts <- txi_tx$counts[, tx_index, drop = FALSE]
  txi_tx_sub$abundance <- txi_tx$abundance[, tx_index, drop = FALSE]
  txi_tx_sub$length <- txi_tx$length[, tx_index, drop = FALSE]

  txi_gene_sub <- txi_gene
  txi_gene_sub$counts <- txi_gene$counts[, gene_index, drop = FALSE]
  txi_gene_sub$abundance <- txi_gene$abundance[, gene_index, drop = FALSE]
  txi_gene_sub$length <- txi_gene$length[, gene_index, drop = FALSE]

  design_formula <- as.formula(comparison_row$design_formula)
  contrast_vector <- c(
    comparison_row$contrast_column,
    comparison_row$numerator_level,
    comparison_row$denominator_level
  )

  output_dir <- file.path(project_root, "outputs", "rna_seq", comparison_id)
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  log_dir <- file.path(output_dir, "logs")
  ensure_dir(table_dir)
  ensure_dir(figure_dir)
  ensure_dir(log_dir)

  write_tsv_base(coldata, file.path(table_dir, paste0(comparison_id, "_sample_table.tsv")))

  dds_tx <- DESeqDataSetFromTximport(txi_tx_sub, colData = coldata, design = design_formula)
  dds_tx <- dds_tx[rowSums(counts(dds_tx)) > 1, ]
  dds_tx <- DESeq(dds_tx)
  tx_results <- as.data.frame(safe_results(dds_tx, contrast_vector, as_flag(comparison_row$use_ihw)))
  tx_results$transcript_id <- rownames(tx_results)
  tx_results <- merge(tx_results, tx2gene, by.x = "transcript_id", by.y = "t_name", all.x = TRUE)
  tx_results <- tx_results[order(tx_results$padj, tx_results$pvalue), ]
  write_tsv_base(tx_results, file.path(table_dir, paste0(comparison_id, "_transcript_results.tsv")))

  tx_aggregated <- NULL
  if (as_flag(comparison_row$aggregate_transcripts)) {
    tx_aggregated <- aggregate_transcript_results(tx_results)
    if (!is.null(tx_aggregated)) {
      write_tsv_base(tx_aggregated, file.path(table_dir, paste0(comparison_id, "_aggregated_transcript_results.tsv")))
    }
  }

  dds_gene <- DESeqDataSetFromTximport(txi_gene_sub, colData = coldata, design = design_formula)
  dds_gene <- dds_gene[rowSums(counts(dds_gene)) > 1, ]
  dds_gene <- DESeq(dds_gene)
  gene_results <- as.data.frame(safe_results(dds_gene, contrast_vector, as_flag(comparison_row$use_ihw)))
  gene_results$gene_name <- rownames(gene_results)
  gene_results <- gene_results[order(gene_results$padj, gene_results$pvalue), ]
  write_tsv_base(gene_results, file.path(table_dir, paste0(comparison_id, "_gene_results.tsv")))

  lfc_threshold <- as.numeric(comparison_row$lfc_threshold)
  sig_gene_results <- gene_results[
    !is.na(gene_results$padj) &
      gene_results$padj < 0.05 &
      abs(gene_results$log2FoldChange) >= lfc_threshold,
  ]
  write_tsv_base(sig_gene_results, file.path(table_dir, paste0(comparison_id, "_gene_results_sig_lfc.tsv")))

  if (!is.null(tx_aggregated)) {
    sig_aggregated <- tx_aggregated[
      !is.na(tx_aggregated$lancaster_p) &
        tx_aggregated$lancaster_p < 0.05,
    ]
    write_tsv_base(sig_aggregated, file.path(table_dir, paste0(comparison_id, "_aggregated_transcript_results_sig.tsv")))
  }

  vsd <- vst(dds_gene, blind = FALSE)
  saveRDS(vsd, file = file.path(output_dir, paste0(comparison_id, "_vsd.rds")))
  build_pca_plot(
    vsd = vsd,
    coldata = transform(coldata, sample_id = rownames(coldata)),
    title_text = comparison_id,
    output_path = file.path(figure_dir, paste0(comparison_id, "_pca.pdf"))
  )

  capture.output(
    sessionInfo(),
    file = file.path(log_dir, paste0(comparison_id, "_session_info.txt"))
  )
}

project_root <- find_project_root()
metadata_dir <- file.path(project_root, "metadata")
args <- parse_args(commandArgs(trailingOnly = TRUE))

sample_sheet <- read_tsv_base(file.path(metadata_dir, "sample_sheet.tsv"))
comparisons <- read_tsv_base(file.path(metadata_dir, "comparisons.tsv"))

rna_samples <- sample_sheet[
  sample_sheet$modality == "RNA_seq" &
    sample_sheet$assay == "RNA" &
    as_flag(sample_sheet$include_manuscript),
]

tx_files <- file.path(project_root, rna_samples$input_counts_path)
names(tx_files) <- rna_samples$sample_id

missing_files <- tx_files[!file.exists(tx_files)]
if (length(missing_files)) {
  stop("Missing RNA count files:\n", paste(missing_files, collapse = "\n"))
}

tx2gene <- make_tx2gene(tx_files[[1]])
txi_tx <- tximport(tx_files, type = "stringtie", tx2gene = tx2gene, txOut = TRUE)
txi_gene <- tximport(tx_files, type = "stringtie", tx2gene = tx2gene, txOut = FALSE)

rna_comparisons <- comparisons[
  comparisons$modality == "RNA_seq" &
    comparisons$assay == "RNA" &
    comparisons$status == "active",
]

if (!is.null(args$comparison_ids)) {
  requested_ids <- unique(args$comparison_ids)
  missing_ids <- setdiff(requested_ids, rna_comparisons$comparison_id)

  if (length(missing_ids)) {
    stop(
      "Requested comparison IDs were not found among active RNA comparisons: ",
      paste(missing_ids, collapse = ", ")
    )
  }

  rna_comparisons <- rna_comparisons[match(requested_ids, rna_comparisons$comparison_id), , drop = FALSE]
}

for (i in seq_len(nrow(rna_comparisons))) {
  run_rna_comparison(
    comparison_row = rna_comparisons[i, ],
    sample_sheet = rna_samples,
    txi_tx = txi_tx,
    txi_gene = txi_gene,
    tx2gene = tx2gene,
    project_root = project_root
  )
}
