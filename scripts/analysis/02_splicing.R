#!/usr/bin/env Rscript

# 02_splicing.R
# Rebuild the manuscript alternative-splicing summary from rMATS JCEC tables.
#
# Required upstream inputs:
#   - data/processed_inputs/splicing/SE.MATS.JCEC.txt
#   - data/processed_inputs/splicing/RI.MATS.JCEC.txt
#   - data/processed_inputs/splicing/A3SS.MATS.JCEC.txt
#   - data/processed_inputs/splicing/A5SS.MATS.JCEC.txt
#   - data/processed_inputs/splicing/MXE.MATS.JCEC.txt
#   - config/curated_features.tsv
#
# Main outputs:
#   - outputs/splicing/tables/
#
# Run after:
#   - 00_validate_inputs.R

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

load_curated_features <- function(project_root) {
  config_path <- file.path(project_root, "config", "curated_features.tsv")
  cfg <- read_tsv_base(config_path)
  cfg[cfg$feature_type == "gene", "feature_id", drop = TRUE]
}

load_rmat_events <- function(base_dir, mode = "JCEC") {
  event_files <- c("SE", "RI", "A3SS", "A5SS", "MXE")
  out <- lapply(event_files, function(event_type) {
    path <- file.path(base_dir, paste0(event_type, ".MATS.", mode, ".txt"))
    df <- read_tsv_base(path)
    df$event_type <- event_type
    df
  })
  names(out) <- event_files
  out
}

bind_significant_events <- function(event_list, fdr_cutoff = 0.05) {
  filtered <- lapply(event_list, function(df) {
    df$FDR <- as.numeric(df$FDR)
    df$IncLevelDifference <- as.numeric(df$IncLevelDifference)
    df[df$FDR < fdr_cutoff, ]
  })
  all_cols <- unique(unlist(lapply(filtered, names)))
  filtered <- lapply(filtered, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols)) {
      for (col in missing_cols) {
        df[[col]] <- NA
      }
    }
    df[, all_cols, drop = FALSE]
  })
  combined <- do.call(rbind, filtered)
  combined$IncLevelDifference <- as.numeric(combined$IncLevelDifference)
  combined
}

make_splicing_summary <- function(events_df, mode_label) {
  event_counts <- as.data.frame(table(events_df$event_type), stringsAsFactors = FALSE)
  names(event_counts) <- c("event_type", "n_events")

  data.frame(
    mode = mode_label,
    metric = c(
      "total_significant_events",
      "unique_genes",
      "positive_inclusion_events",
      "negative_inclusion_events"
    ),
    value = c(
      nrow(events_df),
      length(unique(events_df$geneSymbol)),
      sum(events_df$IncLevelDifference > 0, na.rm = TRUE),
      sum(events_df$IncLevelDifference < 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

write_optional_go <- function(genes, output_path) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(invisible(NULL))
  }

  mapped <- clusterProfiler::bitr(
    unique(genes),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  if (!nrow(mapped)) {
    return(invisible(NULL))
  }

  ego <- clusterProfiler::enrichGO(
    gene = mapped$ENTREZID,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  if (!is.null(ego) && nrow(ego@result)) {
    write_tsv_base(as.data.frame(ego@result), output_path)
  }
}

project_root <- find_project_root()
base_dir <- file.path(project_root, "data", "processed_inputs", "splicing")
output_dir <- file.path(project_root, "outputs", "splicing")
table_dir <- file.path(output_dir, "tables")
log_dir <- file.path(output_dir, "logs")
ensure_dir(table_dir)
ensure_dir(log_dir)

curated_features <- load_curated_features(project_root)
summary_file <- file.path(base_dir, "summary.txt")
summary_table <- read_tsv_base(summary_file)
write_tsv_base(summary_table, file.path(table_dir, "rmats_summary.tsv"))

events_jcec <- bind_significant_events(load_rmat_events(base_dir, "JCEC"))

write_tsv_base(events_jcec, file.path(table_dir, "significant_events_jcec.tsv"))

event_counts_jcec <- as.data.frame(table(events_jcec$event_type), stringsAsFactors = FALSE)
names(event_counts_jcec) <- c("event_type", "n_events")
write_tsv_base(event_counts_jcec, file.path(table_dir, "event_type_counts_jcec.tsv"))

summary_metrics <- make_splicing_summary(events_jcec, "JCEC")
write_tsv_base(summary_metrics, file.path(table_dir, "summary_metrics.tsv"))

curated_hits <- events_jcec[events_jcec$geneSymbol %in% curated_features, ]
write_tsv_base(curated_hits, file.path(table_dir, "curated_splicing_events.tsv"))

directional_gene_sets <- data.frame(
  direction = c("positive_inclusion", "negative_inclusion"),
  n_unique_genes = c(
    length(unique(events_jcec$geneSymbol[events_jcec$IncLevelDifference > 0])),
    length(unique(events_jcec$geneSymbol[events_jcec$IncLevelDifference < 0]))
  ),
  stringsAsFactors = FALSE
)
write_tsv_base(directional_gene_sets, file.path(table_dir, "directional_gene_set_sizes.tsv"))

write_optional_go(unique(events_jcec$geneSymbol), file.path(table_dir, "splicing_gene_go_bp.tsv"))

rna_sig_path <- file.path(
  project_root,
  "outputs", "rna_seq", "rna_main_all_controls",
  "tables", "rna_main_all_controls_gene_results_sig_lfc.tsv"
)
if (file.exists(rna_sig_path)) {
  rna_sig <- read_tsv_base(rna_sig_path)
  overlap <- intersect(unique(events_jcec$geneSymbol), unique(rna_sig$gene_name))
  write_tsv_base(
    data.frame(gene_name = overlap, stringsAsFactors = FALSE),
    file.path(table_dir, "deg_splicing_overlap.tsv")
  )
}

capture.output(
  list(
    generated_tables = list.files(table_dir),
    curated_gene_count = length(unique(curated_hits$geneSymbol))
  ),
  file = file.path(log_dir, "rebuild_log.txt")
)
