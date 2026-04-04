#!/usr/bin/env Rscript

# 09_density_pca_figure.R
# Build the multi-panel PCA summary used in the main manuscript Figure 1C.
#
# Required upstream inputs:
#   - outputs/rna_seq/*/objects/*.rds
#   - processed ATAC and H3K27me3 objects or matrices in the documented public layout
#   - metadata/sample_sheet.tsv
#
# Main outputs:
#   - outputs/figures/main/Figure_1C_density_pca_threepanel.pdf
#   - outputs/figures/main/Figure_1C_density_pca_threepanel.png

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2",             quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("patchwork",           quietly = TRUE)) stop("Package 'patchwork' is required.")
  if (!requireNamespace("cowplot",             quietly = TRUE)) stop("Package 'cowplot' is required.")
  if (!requireNamespace("GenomicRanges",       quietly = TRUE)) stop("Package 'GenomicRanges' is required.")
  if (!requireNamespace("IRanges",             quietly = TRUE)) stop("Package 'IRanges' is required.")
  if (!requireNamespace("S4Vectors",           quietly = TRUE)) stop("Package 'S4Vectors' is required.")
  if (!requireNamespace("BiocGenerics",        quietly = TRUE)) stop("Package 'BiocGenerics' is required.")
  if (!requireNamespace("SummarizedExperiment",quietly = TRUE)) stop("Package 'SummarizedExperiment' is required.")
  library(ggplot2)
  library(patchwork)
})

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, ".nrxn1_project_root"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate project root containing .nrxn1_project_root")
    }
    current <- parent
  }
}

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

require_file_strict <- function(path, label) {
  if (!file.exists(path)) stop("Missing required input file for ", label, ": ", path)
  path
}

has_flag <- function(args, flags) {
  any(args %in% flags)
}

save_plot_pair <- function(plot_obj, pdf_path, png_path = NULL, width = 7, height = 5, dpi = 300) {
  ensure_dir(dirname(pdf_path))
  ggplot2::ggsave(pdf_path, plot = plot_obj, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  message("Saved: ", pdf_path)
  if (!is.null(png_path)) {
    ensure_dir(dirname(png_path))
    ggplot2::ggsave(png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = dpi, bg = "white")
    message("Saved: ", png_path)
  }
}

line_color_values <- c(
  "ctrl10" = "#4472C4",
  "ctrl14" = "#ED7D31",
  "ctrl7" = "#70AD47",
  "patient clone I" = "#FF4B4B",
  "patient clone II" = "#C00000",
  "patient clone III" = "#7A0000"
)

density_size_values <- c("low" = 2.5, "high" = 4)
genotype_shape_values <- c("control" = 16, "patient" = 17)

make_display_line_id <- function(meta_df) {
  clone_map <- c("PI" = "patient clone I", "PII" = "patient clone II", "PIII" = "patient clone III")
  out <- meta_df$line_id
  use_clone <- meta_df$genotype_group == "patient" & !is.na(meta_df$clone_id) & nzchar(meta_df$clone_id)
  out[use_clone] <- clone_map[meta_df$clone_id[use_clone]]
  out[is.na(out)] <- meta_df$line_id[is.na(out)]
  out
}

make_pair_id <- function(meta_df) {
  make_display_line_id(meta_df)
}

compute_row_vars <- function(mat) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::rowVars(mat, na.rm = TRUE))
  }
  apply(mat, 1, stats::var, na.rm = TRUE)
}

finalize_pca_plot_df <- function(plot_df) {
  plot_df$genotype_group <- factor(as.character(plot_df$genotype_group), levels = c("control", "patient"))
  plot_df$density <- factor(as.character(plot_df$density), levels = c("low", "high"))
  if (!"display_line_id" %in% names(plot_df) || all(is.na(plot_df$display_line_id))) {
    plot_df$display_line_id <- make_display_line_id(plot_df)
  }
  plot_df$display_line_id <- factor(as.character(plot_df$display_line_id), levels = names(line_color_values))
  if (!"pair_id" %in% names(plot_df) || all(is.na(plot_df$pair_id))) {
    plot_df$pair_id <- make_pair_id(plot_df)
  }
  plot_df$density_rank <- ifelse(plot_df$density == "low", 1, 2)
  plot_df[order(plot_df$pair_id, plot_df$density_rank, plot_df$sample_id), , drop = FALSE]
}

build_pca_frame <- function(mat, meta_df, top_n, sample_order = NULL) {
  if (!is.null(sample_order)) {
    sample_order <- sample_order[sample_order %in% colnames(mat)]
    mat <- mat[, sample_order, drop = FALSE]
  }
  gene_var <- compute_row_vars(mat)
  gene_var[is.na(gene_var)] <- -Inf
  top_idx <- head(order(gene_var, decreasing = TRUE), min(top_n, nrow(mat)))
  pca <- stats::prcomp(t(mat[top_idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pct_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  plot_df <- merge(meta_df, pca_df, by = "sample_id", all.x = FALSE, all.y = TRUE, sort = FALSE)
  plot_df <- finalize_pca_plot_df(plot_df)
  list(data = plot_df, pct_var = pct_var)
}

write_pca_cache <- function(cache_dir, cache_prefix, pca_res, panel_title) {
  coords_path <- file.path(cache_dir, paste0(cache_prefix, "_coords.tsv"))
  meta_path <- file.path(cache_dir, paste0(cache_prefix, "_meta.tsv"))

  coords_df <- pca_res$data
  factor_cols <- vapply(coords_df, is.factor, logical(1))
  coords_df[factor_cols] <- lapply(coords_df[factor_cols], as.character)
  write.table(coords_df, coords_path, sep = "\t", quote = FALSE, row.names = FALSE)

  meta_df <- data.frame(
    panel_title = panel_title,
    pct_var1 = unname(pca_res$pct_var[[1]]),
    pct_var2 = unname(pca_res$pct_var[[2]]),
    method = if (!is.null(pca_res$method)) pca_res$method else "",
    stringsAsFactors = FALSE
  )
  write.table(meta_df, meta_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_pca_cache <- function(cache_dir, cache_prefix, panel_title_fallback = NULL) {
  coords_path <- file.path(cache_dir, paste0(cache_prefix, "_coords.tsv"))
  meta_path <- file.path(cache_dir, paste0(cache_prefix, "_meta.tsv"))
  coords_df <- read_tsv_base(coords_path)
  if (file.exists(meta_path)) {
    meta_df <- read_tsv_base(meta_path)
    pct_var <- c(meta_df$pct_var1[[1]], meta_df$pct_var2[[1]])
    panel_title <- meta_df$panel_title[[1]]
    method <- meta_df$method[[1]]
  } else {
    pct_var <- c(NA_real_, NA_real_)
    panel_title <- if (!is.null(panel_title_fallback)) panel_title_fallback else cache_prefix
    method <- "cached PCA coordinates only"
  }
  list(
    data = finalize_pca_plot_df(coords_df),
    pct_var = pct_var,
    method = method,
    panel_title = panel_title
  )
}

has_pca_cache <- function(cache_dir, cache_prefix) {
  coords_path <- file.path(cache_dir, paste0(cache_prefix, "_coords.tsv"))
  file.exists(coords_path)
}

load_or_build_pca <- function(cache_dir, cache_prefix, panel_title, build_fun, force_recompute = FALSE, cache_only = FALSE) {
  if (!force_recompute && has_pca_cache(cache_dir, cache_prefix)) {
    message(panel_title, ": using cached PCA tables.")
    return(read_pca_cache(cache_dir, cache_prefix, panel_title_fallback = panel_title))
  }
  if (cache_only) {
    stop(panel_title, ": cache-only mode requested, but cache is incomplete for prefix ", cache_prefix)
  }
  message(panel_title, ": cache missing or recompute requested; rebuilding PCA.")
  pca_res <- build_fun()
  write_pca_cache(cache_dir, cache_prefix, pca_res, panel_title)
  pca_res
}

# Read a MACS2 peak file and return a GRanges with signal scores attached.
# Column 7 (signalValue = fold enrichment) is used as the continuous score.
# Falls back to column 5 (integer score) if col 7 is missing or all-zero.
read_peak_table_with_scores <- function(path, label) {
  require_file_strict(path, label)
  peak_df <- tryCatch(
    read.delim(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) stop("Could not read peak file for ", label, ": ", path, "\n", conditionMessage(e))
  )
  if (ncol(peak_df) < 3) stop("Peak file has fewer than 3 columns for ", label, ": ", path)
  peak_df <- peak_df[
    !is.na(peak_df[[1]]) & !is.na(peak_df[[2]]) & !is.na(peak_df[[3]]),
    , drop = FALSE
  ]
  # Prefer col 7 (signalValue / fold enrichment); fall back to col 5 (score)
  if (ncol(peak_df) >= 7 && is.numeric(peak_df[[7]]) && any(peak_df[[7]] > 0, na.rm = TRUE)) {
    signal <- as.numeric(peak_df[[7]])
  } else {
    signal <- as.numeric(peak_df[[5]])
  }
  signal[is.na(signal) | signal < 0] <- 0
  gr <- GenomicRanges::GRanges(
    seqnames = peak_df[[1]],
    ranges   = IRanges::IRanges(
      start = as.integer(peak_df[[2]]) + 1L,
      end   = as.integer(peak_df[[3]])
    )
  )
  gr$signal <- signal
  gr
}

# Build a continuous signal matrix (peaks x samples) from MACS2 peak files.
# Strategy:
#   1. Form a union (consensus) peak set across all samples via reduce().
#   2. For each sample, find overlaps with the union set and assign that
#      sample's fold-enrichment score to each overlapping union peak.
#      Peaks not detected in a sample get score 0.
#   3. log2(score + 1) transform to reduce skew before PCA.
# This gives a continuous matrix that captures real signal variation —
# far superior to the binary presence/absence matrix for PCA purposes.
build_peak_score_matrix <- function(meta_df, project_root, label) {
  message(label, ": building union-peak signal score matrix from MACS2 output.")
  peak_list <- lapply(seq_len(nrow(meta_df)), function(i) {
    peak_path <- file.path(project_root, meta_df$peak_path[[i]])
    read_peak_table_with_scores(peak_path, paste0(label, " / ", meta_df$sample_id[[i]]))
  })
  names(peak_list) <- meta_df$sample_id

  # Union peak set
  all_peaks <- suppressWarnings(Reduce(c, peak_list))
  consensus  <- suppressWarnings(GenomicRanges::reduce(all_peaks, ignore.strand = TRUE))
  message(label, ": union peak set contains ", length(consensus), " regions across ", nrow(meta_df), " samples.")

  # Score matrix: rows = consensus peaks, cols = samples
  score_matrix <- vapply(
    peak_list,
    function(gr) {
      scores <- rep(0, length(consensus))
      hits <- GenomicRanges::findOverlaps(consensus, gr, ignore.strand = TRUE)
      q_hits <- S4Vectors::queryHits(hits)
      s_hits <- S4Vectors::subjectHits(hits)
      # If multiple peaks overlap one consensus region, take the max signal
      for (idx in seq_along(q_hits)) {
        scores[q_hits[idx]] <- max(scores[q_hits[idx]], gr$signal[s_hits[idx]])
      }
      scores
    },
    numeric(length(consensus))
  )
  colnames(score_matrix) <- meta_df$sample_id
  rownames(score_matrix) <- paste0(
    as.character(GenomicRanges::seqnames(consensus)), ":",
    BiocGenerics::start(consensus), "-",
    BiocGenerics::end(consensus)
  )
  # log2(x + 1) transform — reduces right-skew from high fold-enrichment outliers
  log2(score_matrix + 1)
}

build_peak_pca_from_metadata <- function(meta_df, project_root, label, top_n = 5000) {
  if (nrow(meta_df) < 4) {
    warning(label, ": fewer than 4 samples (n=", nrow(meta_df), "). Plot will still be produced.")
  }
  mat <- build_peak_score_matrix(meta_df, project_root, label)
  pca_res <- build_pca_frame(mat, meta_df, top_n = top_n, sample_order = meta_df$sample_id)
  pca_res$method <- "MACS2 fold-enrichment score matrix (log2 + union peaks)"
  pca_res
}

build_modality_pca_plot <- function(plot_df, pct_var, panel_title) {
  line_df <- plot_df[duplicated(plot_df$pair_id) | duplicated(plot_df$pair_id, fromLast = TRUE), , drop = FALSE]
  x_label <- if (length(pct_var) >= 1 && is.finite(pct_var[[1]])) {
    sprintf("PC1 (%.1f%% variance)", pct_var[[1]])
  } else {
    "PC1"
  }
  y_label <- if (length(pct_var) >= 2 && is.finite(pct_var[[2]])) {
    sprintf("PC2 (%.1f%% variance)", pct_var[[2]])
  } else {
    "PC2"
  }
  ggplot2::ggplot(plot_df, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_path(
      data = line_df,
      ggplot2::aes(group = pair_id, color = display_line_id),
      linewidth = 0.8,
      alpha = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = display_line_id, shape = genotype_group, size = density),
      alpha = 0.95,
      stroke = 0.3
    ) +
    ggplot2::scale_color_manual(values = line_color_values, drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = genotype_shape_values, name = NULL, drop = FALSE) +
    ggplot2::scale_size_manual(
      values = density_size_values,
      breaks = c("low", "high"),
      labels = c("low density", "high density"),
      name = NULL,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = panel_title,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 8.5),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8),
      legend.position = "right"
    )
}

build_reference_legend_plot <- function() {
  sample_df <- data.frame(
    x = 0.08,
    y = c(0.925, 0.885, 0.845, 0.790, 0.750, 0.710),
    label = c("ctrl10", "ctrl14", "ctrl7", "Patient I", "Patient II", "Patient III"),
    display_line_id = factor(
      c("ctrl10", "ctrl14", "ctrl7", "patient clone I", "patient clone II", "patient clone III"),
      levels = names(line_color_values)
    ),
    genotype_group = factor(c("control", "control", "control", "patient", "patient", "patient"),
                            levels = c("control", "patient"))
  )
  density_df <- data.frame(
    x = 0.08,
    y = c(0.535, 0.470),
    label = c("Low density", "High density"),
    density = factor(c("low", "high"), levels = c("low", "high"))
  )

  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.02, y = 0.975, label = "Samples", fontface = "bold", hjust = 0, size = 2.7) +
    ggplot2::geom_point(
      data = sample_df,
      ggplot2::aes(x = x, y = y, color = display_line_id, shape = genotype_group),
      size = 1.9,
      stroke = 0.3
    ) +
    ggplot2::geom_text(
      data = sample_df,
      ggplot2::aes(x = 0.145, y = y, label = label),
      hjust = 0,
      size = 2.25
    ) +
    ggplot2::annotate("text", x = 0.02, y = 0.605, label = "Density", fontface = "bold", hjust = 0, size = 2.7) +
    ggplot2::geom_point(
      data = density_df,
      ggplot2::aes(x = x, y = y, size = density),
      shape = 16,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = density_df,
      ggplot2::aes(x = 0.145, y = y, label = label),
      hjust = 0,
      size = 2.25
    ) +
    ggplot2::scale_color_manual(values = line_color_values, guide = "none") +
    ggplot2::scale_shape_manual(values = genotype_shape_values, guide = "none") +
    ggplot2::scale_size_manual(values = density_size_values, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0.12, 1), expand = FALSE, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))
}

project_root <- find_project_root()
script_args <- commandArgs(trailingOnly = TRUE)
force_recompute <- has_flag(script_args, c("--force-recompute", "--recompute"))
cache_only <- has_flag(script_args, c("--cache-only"))
if (force_recompute && cache_only) {
  stop("Choose either --force-recompute or --cache-only, not both.")
}
res_dir <- file.path(project_root, "outputs")
fig_main <- file.path(res_dir, "figures", "main")
fig_tables <- file.path(res_dir, "figures", "tables")
ensure_dir(fig_main)
ensure_dir(fig_tables)

sample_sheet_path <- require_file_strict(file.path(project_root, "metadata", "sample_sheet.tsv"), "sample sheet")
sample_sheet <- read_tsv_base(sample_sheet_path)

# RNA panel
message("--- Figure 1D Panel A: RNA-seq PCA ---")
rna_meta <- sample_sheet[sample_sheet$modality == "RNA_seq", , drop = FALSE]
rna_pca <- load_or_build_pca(
  cache_dir = fig_tables,
  cache_prefix = "Figure_1D_rna_pca",
  panel_title = "RNA-seq",
  force_recompute = force_recompute,
  cache_only = cache_only,
  build_fun = function() {
    rna_vsd_path <- require_file_strict(file.path(res_dir, "rna_seq", "rna_main_all_controls", "rna_main_all_controls_vsd.rds"), "RNA VST object")
    rna_vsd <- readRDS(rna_vsd_path)
    rna_mat <- SummarizedExperiment::assay(rna_vsd)
    rna_pca_res <- build_pca_frame(rna_mat, rna_meta, top_n = 500, sample_order = rna_meta$sample_id)
    rna_pca_res$method <- "DESeq2 variance-stabilized expression; top 500 variable genes"
    rna_pca_res
  }
)
rna_plot <- build_modality_pca_plot(rna_pca$data, rna_pca$pct_var, rna_pca$panel_title)

# ATAC panel
message("--- Figure 1D Panel B: ATAC-seq PCA ---")
atac_meta <- sample_sheet[sample_sheet$modality == "ATAC_seq", , drop = FALSE]
atac_pca <- load_or_build_pca(
  cache_dir = fig_tables,
  cache_prefix = "Figure_1D_atac_pca",
  panel_title = "ATAC-seq",
  force_recompute = force_recompute,
  cache_only = cache_only,
  build_fun = function() {
    build_peak_pca_from_metadata(atac_meta, project_root, label = "ATAC-seq", top_n = 5000)
  }
)
atac_plot <- build_modality_pca_plot(atac_pca$data, atac_pca$pct_var, atac_pca$panel_title)

# ChIP panel
message("--- Figure 1D Panel C: H3K27me3 ChIP-seq PCA ---")
chip_meta <- sample_sheet[
  sample_sheet$modality == "ChIP_seq" &
    sample_sheet$assay == "H3K27me3" &
    grepl("ME", sample_sheet$sample_id),
  ,
  drop = FALSE
]
chip_pca <- load_or_build_pca(
  cache_dir = fig_tables,
  cache_prefix = "Figure_1D_chip_pca",
  panel_title = "H3K27me3 ChIP-seq",
  force_recompute = force_recompute,
  cache_only = cache_only,
  build_fun = function() {
    build_peak_pca_from_metadata(chip_meta, project_root, label = "H3K27me3 ChIP-seq", top_n = 5000)
  }
)
chip_plot <- build_modality_pca_plot(chip_pca$data, chip_pca$pct_var, chip_pca$panel_title)

# Combined figure
message("--- Figure 1D Combined three-panel PCA ---")
legend_plot <- build_reference_legend_plot()
panel_row <- cowplot::plot_grid(
  rna_plot + ggplot2::theme(legend.position = "none"),
  atac_plot + ggplot2::theme(legend.position = "none"),
  chip_plot + ggplot2::theme(legend.position = "none"),
  nrow = 1,
  align = "h"
)
combined_plot <- cowplot::plot_grid(
  cowplot::plot_grid(panel_row, legend_plot, nrow = 1, rel_widths = c(1, 0.075)),
  ncol = 1
)

save_plot_pair(
  rna_plot + ggplot2::theme(legend.position = "none"),
  file.path(fig_main, "Figure_1D_rna_pca.pdf"),
  NULL,
  width = 4.5,
  height = 4.5
)

save_plot_pair(
  atac_plot + ggplot2::theme(legend.position = "none"),
  file.path(fig_main, "Figure_1D_atac_pca.pdf"),
  NULL,
  width = 4.5,
  height = 4.5
)

save_plot_pair(
  chip_plot + ggplot2::theme(legend.position = "none"),
  file.path(fig_main, "Figure_1D_chip_pca.pdf"),
  NULL,
  width = 4.5,
  height = 4.5
)

save_plot_pair(
  combined_plot,
  file.path(fig_main, "Figure_1D_density_pca_threepanel.pdf"),
  file.path(fig_main, "Figure_1D_density_pca_threepanel.png"),
  width = 13,
  height = 4.5
)

save_plot_pair(
  combined_plot,
  file.path(fig_main, "Figure_1C_density_pca_threepanel.pdf"),
  file.path(fig_main, "Figure_1C_density_pca_threepanel.png"),
  width = 13,
  height = 4.5
)

message("ATAC PCA method: ", atac_pca$method)
message("ChIP PCA method: ", chip_pca$method)
