#!/usr/bin/env Rscript

# 09s_build_supp_qc_panels.R
# Build the Supplementary Figure S1 QC panels from compact QC summary tables.
#
# Required upstream inputs:
#   - data/processed_inputs/qc/Table_S1_metadata.csv
#   - data/processed_inputs/qc/Table_S2_seq_qc.csv
#   - data/processed_inputs/qc/Table_S3_atac_qc.csv
#
# Optional richer QC inputs:
#   - data/processed_inputs/qc/Table_S3b_atac_fragment_profiles.tsv.gz
#   - data/processed_inputs/qc/Table_S3c_atac_tss_profiles.tsv.gz
#
# Main outputs:
#   - outputs/figures/supp/

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(SummarizedExperiment)
  library(grid)
})

locate_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    sentinel <- file.path(current, ".nrxn1_project_root")
    if (file.exists(sentinel)) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate project root via .nrxn1_project_root sentinel.")
    }
    current <- parent
  }
}

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

find_external_root <- function() {
  # Optional large ATAC inputs (fragment/TSS BAMs and the computeMatrix output)
  # can live outside the repository. Set NRXN1_EXTERNAL_ROOT to that directory to
  # enable the extended ATAC profile panels. When it is unset or missing, the
  # workflow skips those optional panels.
  root <- Sys.getenv("NRXN1_EXTERNAL_ROOT", unset = NA_character_)
  if (is.na(root) || !nzchar(root) || !file.exists(root)) {
    return(NA_character_)
  }
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

theme_qc <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#f3f3f3", colour = "#d0d0d0"),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(face = "bold"),
      axis.title.y = element_blank(),
      plot.subtitle = element_text(size = rel(0.9)),
      legend.margin = margin(t = 2, r = 0, b = 0, l = 0)
    )
}

plot_chip_frip <- function(chip_qc, metadata) {
  dt <- merge(copy(chip_qc), metadata[, .(sample_name, sex, seeding_density)], by = "sample_name", all.x = TRUE)
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  dt[, has_matched_input := factor(has_matched_input, levels = c("yes", "no"), labels = c("Matched input", "No matched input"))]
  dt[, sample_name := factor(sample_name, levels = dt[order(frip_score, sample_name), sample_name])]

  ggplot(dt, aes(x = frip_score, y = sample_name, colour = genotype_label, shape = has_matched_input)) +
    geom_segment(aes(x = 0, xend = frip_score, y = sample_name, yend = sample_name), inherit.aes = FALSE, colour = "#d7d7d7", linewidth = 0.8) +
    geom_vline(xintercept = 0.05, colour = "#bfbfbf", linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 2.7) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    labs(
      title = "H3K27me3 ChIP-seq FRiP score",
      subtitle = "Counts, mapping rates, and matched-input details are provided in supplementary QC tables",
      x = "FRiP score",
      colour = NULL,
      shape = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    theme_qc(10)
}

plot_chip_frip_lollipop <- function(chip_qc, metadata) {
  dt <- merge(copy(chip_qc), metadata[, .(sample_name, sex, seeding_density)], by = "sample_name", all.x = TRUE)
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  dt[, has_matched_input := factor(has_matched_input, levels = c("yes", "no"), labels = c("Matched input", "No matched input"))]
  dt[, sample_name := factor(sample_name, levels = dt[order(frip_score, sample_name), sample_name])]

  ggplot(dt, aes(x = frip_score, y = sample_name, colour = genotype_label, shape = has_matched_input)) +
    geom_segment(aes(x = 0, xend = frip_score, y = sample_name, yend = sample_name), inherit.aes = FALSE, colour = "#d7d7d7", linewidth = 1.0) +
    geom_vline(xintercept = 0.05, colour = "#bfbfbf", linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 3.0) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = "H3K27me3 ChIP-seq FRiP score",
      subtitle = "Ranked lollipop view",
      x = "FRiP score",
      colour = NULL,
      shape = NULL
    ) +
    theme_qc(10)
}

plot_atac_fragment_summary <- function(atac_qc) {
  dt <- melt(
    copy(atac_qc),
    id.vars = c("sample_name", "genotype", "sex"),
    measure.vars = c("sub_nucleosomal_pct", "mono_nucleosomal_pct"),
    variable.name = "fragment_class",
    value.name = "pct"
  )
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  dt[, fragment_class := factor(
    fragment_class,
    levels = c("sub_nucleosomal_pct", "mono_nucleosomal_pct"),
    labels = c("Sub-nucleosomal", "Mono-nucleosomal")
  )]
  dt[, sample_name := factor(sample_name, levels = rev(atac_qc[order(genotype, sample_name), sample_name]))]

  ggplot(dt, aes(x = pct, y = sample_name, colour = genotype_label, shape = fragment_class)) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    scale_shape_manual(values = c(16, 17)) +
    labs(
      title = "ATAC fragment composition",
      subtitle = "Sub-nucleosomal fragments dominate, with preserved nucleosomal signal across all samples",
      x = "Fragment fraction (%)",
      colour = NULL,
      shape = NULL
    ) +
    theme_qc(10)
}

build_atac_profile_data <- function(external_root, metadata) {
  project_root <- locate_project_root()
  cache_dir <- file.path(project_root, "data", "processed_inputs", "qc")
  frag_cache <- file.path(cache_dir, "Table_S3b_atac_fragment_profiles.tsv.gz")
  tss_cache <- file.path(cache_dir, "Table_S3c_atac_tss_profiles.tsv.gz")

  if (file.exists(frag_cache) && file.exists(tss_cache)) {
    return(list(
      frag_hist = fread(frag_cache),
      profiles = fread(tss_cache)
    ))
  }

  atac <- copy(metadata[assay == "ATAC-seq"])
  atac[, bam_file := file.path(external_root, "ATAC_processed", "ATAC_seq", paste0(sample_name, ".bam"))]
  atac <- atac[file.exists(bam_file)]
  if (!nrow(atac)) {
    stop("No ATAC BAM files found on mounted external drive.")
  }

  frag_hist <- rbindlist(lapply(seq_len(nrow(atac)), function(i) {
    bam <- atac$bam_file[i]
    hist_cmd <- sprintf(
      "samtools view %s | awk 'function abs(x){return x<0?-x:x} $9!=0 {l=abs($9); if (l<=800) {c[l]++; n++; if (n>=400000) exit}} END {for (i=1;i<=800;i++) if (c[i]>0) print i, c[i]}'",
      shQuote(bam)
    )
    frag <- fread(cmd = hist_cmd, header = FALSE, col.names = c("length", "count"))
    frag[, .(sample_name = atac$sample_name[i], length, density = count / sum(count))]
  }))

  matrix_path <- file.path(external_root, "ATAC_processed", "ATAC_seq", "ComputeMatrix", "matrix1_reference.gz")
  if (!file.exists(matrix_path)) {
    stop("ATAC TSS matrix1_reference.gz not found on mounted external drive.")
  }
  con <- gzfile(matrix_path, "rt")
  header <- sub("^@", "", readLines(con, n = 1))
  close(con)
  matrix_meta <- fromJSON(header)
  matrix_dt <- fread(cmd = sprintf("gzip -dc %s | tail -n +2", shQuote(matrix_path)), header = FALSE)
  sample_labels <- matrix_meta$sample_labels
  sample_bounds <- matrix_meta$sample_boundaries
  bin_size <- as.numeric(matrix_meta[["bin size"]][1])
  x_pos <- seq(-3000 + bin_size / 2, 3000 - bin_size / 2, by = bin_size)

  profiles <- rbindlist(lapply(seq_along(sample_labels), function(i) {
    cols <- (7 + sample_bounds[i]):(6 + sample_bounds[i + 1])
    mat <- as.matrix(matrix_dt[, ..cols])
    data.table(sample_name = sample_labels[i], pos = x_pos, signal = colMeans(mat, na.rm = TRUE))
  }))

  fwrite(frag_hist, frag_cache)
  fwrite(profiles, tss_cache)

  list(frag_hist = frag_hist, profiles = profiles)
}

plot_atac_fragment_profile <- function(frag_hist, metadata) {
  dt <- merge(copy(frag_hist), metadata[, .(sample_name, genotype)], by = "sample_name", all.x = TRUE)
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]

  ggplot(dt[length <= 800], aes(length, density, colour = genotype_label, group = sample_name)) +
    geom_line(alpha = 0.72, linewidth = 0.55) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    labs(
      title = "ATAC fragment size distribution",
      subtitle = "Sub-nucleosomal and mono-nucleosomal peaks are preserved across all libraries",
      x = "Fragment length (bp)",
      y = "Relative density",
      colour = NULL
    ) +
    theme_qc(10)
}

plot_atac_tss_profile <- function(profiles, metadata) {
  dt <- merge(copy(profiles), metadata[, .(sample_name, genotype)], by = "sample_name", all.x = TRUE)
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  dt <- dt[abs(pos) <= 2000]
  grp <- dt[, .(
    mean_signal = mean(signal, na.rm = TRUE),
    se = sd(signal, na.rm = TRUE) / sqrt(.N)
  ), by = .(genotype_label, pos)]

  ggplot(grp, aes(pos, mean_signal, colour = genotype_label, fill = genotype_label)) +
    geom_ribbon(aes(ymin = mean_signal - se, ymax = mean_signal + se), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.85) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    scale_fill_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    labs(
      title = "ATAC TSS enrichment",
      subtitle = "Aggregate profiles are comparable between control and NRXN1alpha deletion libraries",
      x = "Distance from TSS (bp)",
      y = "Mean signal",
      colour = NULL,
      fill = NULL
    ) +
    theme_qc(10)
}

plot_atac_tss_score <- function(atac_qc) {
  dt <- copy(atac_qc)
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  dt[, sample_name := factor(sample_name, levels = rev(dt[order(genotype, sample_name), sample_name]))]

  ggplot(dt, aes(x = tss_enrichment_score, y = sample_name, colour = genotype_label, shape = sex)) +
    geom_vline(xintercept = 7, colour = "#bfbfbf", linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    labs(
      title = "ATAC TSS enrichment",
      subtitle = "All libraries exceed a typical TSS enrichment threshold",
      x = "TSS enrichment score",
      colour = NULL,
      shape = "Sex"
    ) +
    theme_qc(10)
}

plot_atac_frip <- function(atac_qc) {
  dt <- copy(atac_qc)
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  dt[, sample_name := factor(sample_name, levels = dt[order(frip_score, sample_name), sample_name])]

  ggplot(dt, aes(x = frip_score, y = sample_name, colour = genotype_label, shape = sex)) +
    geom_segment(aes(x = 0, xend = frip_score, y = sample_name, yend = sample_name), inherit.aes = FALSE, colour = "#d7d7d7", linewidth = 0.8) +
    geom_vline(xintercept = 0.2, colour = "#bfbfbf", linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    labs(
      title = "ATAC FRiP score",
      subtitle = "All libraries exceed a typical minimum FRiP threshold",
      x = "FRiP score",
      colour = NULL,
      shape = "Sex"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    theme_qc(10)
}

plot_atac_frip_lollipop <- function(atac_qc) {
  dt <- copy(atac_qc)
  dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  dt[, sample_name := factor(sample_name, levels = dt[order(frip_score, sample_name), sample_name])]

  ggplot(dt, aes(x = frip_score, y = sample_name, colour = genotype_label, shape = sex)) +
    geom_segment(aes(x = 0, xend = frip_score, y = sample_name, yend = sample_name), inherit.aes = FALSE, colour = "#d7d7d7", linewidth = 1.0) +
    geom_vline(xintercept = 0.2, colour = "#bfbfbf", linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 3.0) +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = "ATAC FRiP score",
      subtitle = "Ranked lollipop view",
      x = "FRiP score",
      colour = NULL,
      shape = "Sex"
    ) +
    theme_qc(10)
}

plot_combined_frip_panel <- function(chip_qc, atac_qc, metadata) {
  chip_dt <- merge(copy(chip_qc), metadata[, .(sample_name, sex)], by = "sample_name", all.x = TRUE)
  chip_dt[, assay := "H3K27me3 ChIP-seq"]
  chip_dt[, threshold := 0.05]
  chip_dt[, shape_group := fifelse(has_matched_input == "yes", "Matched input", "No matched input")]
  chip_dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  chip_dt <- chip_dt[, .(sample_name, assay, frip_score, threshold, shape_group, genotype_label)]

  atac_dt <- copy(atac_qc)
  atac_dt[, assay := "ATAC-seq"]
  atac_dt[, threshold := 0.2]
  atac_dt[, shape_group := ifelse(sex == "female", "Female", "Male")]
  atac_dt[, genotype_label := fifelse(genotype == "control", "Control", "NRXN1alpha deletion")]
  atac_dt <- atac_dt[, .(sample_name, assay, frip_score, threshold, shape_group, genotype_label)]

  dt <- rbindlist(list(chip_dt, atac_dt), use.names = TRUE)
  dt[, assay := factor(assay, levels = c("H3K27me3 ChIP-seq", "ATAC-seq"))]
  dt[, sample_key := paste0(as.character(assay), "|||", sample_name)]
  ordered_levels <- unlist(lapply(split(dt, dt$assay), function(sub) sub[order(frip_score, sample_name), sample_key]))
  dt[, sample_key := factor(sample_key, levels = ordered_levels)]

  ggplot(dt, aes(x = frip_score, y = sample_key, colour = genotype_label, shape = shape_group)) +
    geom_segment(aes(x = 0, xend = frip_score, y = sample_key, yend = sample_key), inherit.aes = FALSE, colour = "#d7d7d7", linewidth = 0.9) +
    geom_vline(aes(xintercept = threshold), colour = "#bfbfbf", linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 2.8) +
    facet_wrap(~ assay, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = c("Control" = "#4c78a8", "NRXN1alpha deletion" = "#e45756")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
    labs(
      title = "FRiP score comparison across assays",
      subtitle = "Assay-specific thresholds shown as dashed lines",
      x = "FRiP score",
      colour = NULL,
      shape = NULL
    ) +
    theme_qc(10)
}

plot_rna_correlation <- function(corr_mat) {
  ord <- hclust(as.dist(1 - corr_mat), method = "average")$order
  corr_ord <- corr_mat[ord, ord]
  dt <- as.data.table(as.table(corr_ord))
  setnames(dt, c("sample_row", "sample_col", "corr"))
  dt[, sample_row := factor(sample_row, levels = rev(rownames(corr_ord)))]
  dt[, sample_col := factor(sample_col, levels = colnames(corr_ord))]

  ggplot(dt, aes(x = sample_col, y = sample_row, fill = corr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0.985) +
    coord_fixed() +
    labs(
      title = "RNA-seq sample correlation",
      subtitle = "Hierarchical ordering of Pearson correlations across all RNA-seq libraries",
      x = NULL,
      y = NULL,
      fill = "r"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    )
}

main <- function() {
  project_root <- locate_project_root()
  qc_table_dir <- file.path(project_root, "data", "processed_inputs", "qc")
  fig_dir <- file.path(manuscript_root, "results", "figures", "supp", "qc_panels")
  ensure_dir(fig_dir)

  atac_qc <- fread(file.path(qc_table_dir, "Table_S3_atac_qc.csv"))
  chip_qc <- fread(file.path(qc_table_dir, "Table_S4_chip_qc_EXCLUDED.csv"))
  meta <- fread(file.path(qc_table_dir, "Table_S1_metadata.csv"))
  vsd <- readRDS(file.path(manuscript_root, "results", "rna_seq", "rna_main_all_controls", "rna_main_all_controls_vsd.rds"))
  corr_mat <- cor(assay(vsd), method = "pearson")
  external_root <- find_external_root()

  p_corr <- plot_rna_correlation(corr_mat)
  p_chip <- plot_chip_frip(chip_qc, meta)
  if (!is.na(external_root)) {
    atac_profiles <- build_atac_profile_data(external_root, meta)
    p_frag <- plot_atac_fragment_profile(atac_profiles$frag_hist, meta[assay == "ATAC-seq"])
    p_tss <- plot_atac_tss_profile(atac_profiles$profiles, meta[assay == "ATAC-seq"])
  } else {
    p_frag <- plot_atac_fragment_summary(atac_qc)
    p_tss <- plot_atac_tss_score(atac_qc)
  }
  p_frip <- plot_atac_frip(atac_qc)
  p_chip_lollipop <- plot_chip_frip_lollipop(chip_qc, meta)
  p_atac_lollipop <- plot_atac_frip_lollipop(atac_qc)
  p_frip_combined <- plot_combined_frip_panel(chip_qc, atac_qc, meta)

  out_corr <- file.path(fig_dir, "Figure_S1A_rna_sample_correlation.pdf")
  out_chip <- file.path(fig_dir, "Figure_S1B_chip_frip_scores.pdf")
  out_frag <- file.path(fig_dir, "Figure_S1C_atac_fragment_composition.pdf")
  out_tss <- file.path(fig_dir, "Figure_S1D_atac_tss_enrichment_scores.pdf")
  out_frip <- file.path(fig_dir, "Figure_S1E_atac_frip_scores.pdf")
  out_chip_lollipop <- file.path(fig_dir, "Figure_S1B_chip_frip_lollipop.pdf")
  out_atac_lollipop <- file.path(fig_dir, "Figure_S1E_atac_frip_lollipop.pdf")
  out_frip_combined <- file.path(fig_dir, "Figure_S1F_combined_frip_panel.pdf")

  ggsave(out_corr, p_corr, width = 6.4, height = 5.6, useDingbats = FALSE)
  ggsave(out_chip, p_chip, width = 6.2, height = 4.4, useDingbats = FALSE)
  ggsave(out_frag, p_frag, width = 6.2, height = 4.4, useDingbats = FALSE)
  ggsave(out_tss, p_tss, width = 6.0, height = 4.4, useDingbats = FALSE)
  ggsave(out_frip, p_frip, width = 6.0, height = 4.4, useDingbats = FALSE)
  ggsave(out_chip_lollipop, p_chip_lollipop, width = 6.2, height = 4.4, useDingbats = FALSE)
  ggsave(out_atac_lollipop, p_atac_lollipop, width = 6.0, height = 4.4, useDingbats = FALSE)
  ggsave(out_frip_combined, p_frip_combined, width = 9.0, height = 4.6, useDingbats = FALSE)

  cat("Saved:\n")
  cat("  ", out_corr, "\n")
  cat("  ", out_chip, "\n")
  cat("  ", out_frag, "\n")
  cat("  ", out_tss, "\n")
  cat("  ", out_frip, "\n")
  cat("  ", out_chip_lollipop, "\n")
  cat("  ", out_atac_lollipop, "\n")
  cat("  ", out_frip_combined, "\n")
}

main()
