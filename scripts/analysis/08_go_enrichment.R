#!/usr/bin/env Rscript

# 08_go_enrichment.R
# Run the topGO enrichment analyses used in the RNA, RNA-ATAC, RNA-H3K27me3,
# and triple-overlap manuscript panels.
#
# Required upstream inputs:
#   - outputs/rna_seq/rna_main_all_controls/tables/
#   - outputs/integration/tables/
#
# Main outputs:
#   - outputs/go_enrichment/tables/
#
# Run after:
#   - 01_rna_expression.R
#   - 05_integration.R

suppressPackageStartupMessages({
  req_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  if (!req_pkg("topGO")) stop("Package 'topGO' is required for script 08.")
  if (!req_pkg("org.Hs.eg.db")) stop("Package 'org.Hs.eg.db' is required for script 08.")
  if (!req_pkg("GO.db")) stop("Package 'GO.db' is required for script 08.")
  library(topGO)
  library(GO.db)
})

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, ".nrxn1_project_root"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate project root. Ensure .nrxn1_project_root exists.")
    }
    current <- parent
  }
}

read_tsv_base <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_base <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

parse_topgo_pvalue <- function(x) {
  x <- gsub("^<\\s*", "", trimws(as.character(x)))
  suppressWarnings(as.numeric(x))
}

collect_common_genes <- function(go_data, go_ids, query_genes) {
  term_genes <- topGO::genesInTerm(go_data, go_ids)
  vapply(
    term_genes,
    function(sig) paste(sort(intersect(query_genes, sig)), collapse = " "),
    character(1)
  )
}

run_topgo_bp <- function(gene_vec, label, node_size = 5) {
  gene_vec <- sort(unique(gene_vec[!is.na(gene_vec) & nzchar(gene_vec)]))
  if (!length(gene_vec)) {
    stop("No genes available for topGO: ", label)
  }

  background_map <- topGO::annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
  background_genes <- sort(unique(unlist(background_map)))
  gene_factor <- factor(ifelse(background_genes %in% gene_vec, 1L, 0L))
  names(gene_factor) <- background_genes

  go_data <- methods::new(
    "topGOdata",
    ontology = "BP",
    allGenes = gene_factor,
    nodeSize = node_size,
    annot = topGO::annFUN.org,
    mapping = "org.Hs.eg.db",
    ID = "symbol"
  )

  result_classic <- topGO::runTest(go_data, algorithm = "classic", statistic = "fisher")
  result_weight01 <- topGO::runTest(go_data, algorithm = "weight01", statistic = "fisher")

  all_res <- topGO::GenTable(
    go_data,
    classicFisher = result_classic,
    weight01Fisher = result_weight01,
    orderBy = "weight01Fisher",
    topNodes = length(topGO::usedGO(go_data)),
    numChar = 1000
  )
  all_res <- all_res[all_res$classicFisher != "1.00000", , drop = FALSE]
  all_res$common <- collect_common_genes(go_data, all_res$GO.ID, gene_vec)
  all_res$classicFisher_num <- parse_topgo_pvalue(all_res$classicFisher)
  all_res$weight01Fisher_num <- parse_topgo_pvalue(all_res$weight01Fisher)
  all_res$gene_ratio <- as.numeric(all_res$Significant) / pmax(as.numeric(all_res$Annotated), 1)
  all_res$enrichment_score <- -log10(pmax(all_res$weight01Fisher_num, .Machine$double.xmin))
  all_res$classic_enrichment_score <- -log10(pmax(all_res$classicFisher_num, .Machine$double.xmin))
  all_res <- all_res[order(all_res$weight01Fisher_num, all_res$classicFisher_num, all_res$Term), , drop = FALSE]

  list(
    gene_count = length(gene_vec),
    background_gene_count = length(background_genes),
    table = all_res
  )
}

project_root <- find_project_root()
result_dir <- file.path(project_root, "outputs", "go_enrichment")
table_dir <- file.path(result_dir, "tables")
log_dir <- file.path(result_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

go_jobs <- list(
  list(
    comparison_id = "rna_main_up",
    label = "Main RNA upregulated genes",
    job_type = "rna_directional",
    input_path = file.path(project_root, "outputs", "rna_seq", "rna_main_all_controls", "tables", "rna_main_all_controls_gene_results_sig_lfc.tsv"),
    direction = "up"
  ),
  list(
    comparison_id = "rna_main_down",
    label = "Main RNA downregulated genes",
    job_type = "rna_directional",
    input_path = file.path(project_root, "outputs", "rna_seq", "rna_main_all_controls", "tables", "rna_main_all_controls_gene_results_sig_lfc.tsv"),
    direction = "down"
  ),
  list(
    comparison_id = "rna_female_up",
    label = "Female-only RNA upregulated genes",
    job_type = "rna_directional",
    input_path = file.path(project_root, "outputs", "rna_seq", "rna_supp_female_only", "tables", "rna_supp_female_only_gene_results_sig_lfc.tsv"),
    direction = "up"
  ),
  list(
    comparison_id = "rna_female_down",
    label = "Female-only RNA downregulated genes",
    job_type = "rna_directional",
    input_path = file.path(project_root, "outputs", "rna_seq", "rna_supp_female_only", "tables", "rna_supp_female_only_gene_results_sig_lfc.tsv"),
    direction = "down"
  ),
  list(
    comparison_id = "rna_chip_main",
    label = "Main RNA-ChIP overlap",
    job_type = "overlap",
    input_path = file.path(project_root, "outputs", "integration", "tables", "rna_chip_main_overlap.tsv")
  ),
  list(
    comparison_id = "rna_chip_female",
    label = "Female-only RNA-ChIP overlap",
    job_type = "overlap",
    input_path = file.path(project_root, "outputs", "integration", "tables", "rna_chip_female_overlap.tsv")
  ),
  list(
    comparison_id = "triple_main",
    label = "Main RNA-ChIP-ATAC triple overlap",
    job_type = "overlap",
    input_path = file.path(project_root, "outputs", "integration", "tables", "triple_overlap_main.tsv")
  ),
  list(
    comparison_id = "triple_female",
    label = "Female-only RNA-ChIP-ATAC triple overlap",
    job_type = "overlap",
    input_path = file.path(project_root, "outputs", "integration", "tables", "triple_overlap_female.tsv")
  ),
  # ── Figure 4C: RNA-ATAC concordant genes (both up) ──────────────────────────
  list(
    comparison_id = "rna_atac_main_up",
    label = "RNA-ATAC concordant: both increased in patient",
    job_type = "rna_atac_directional",
    input_path = file.path(project_root, "outputs", "integration", "tables", "rna_atac_main_overlap.tsv"),
    direction_filter = "up"
  ),
  # ── Figure 4C: RNA-ATAC concordant genes (both down) ────────────────────────
  list(
    comparison_id = "rna_atac_main_down",
    label = "RNA-ATAC concordant: both decreased in patient",
    job_type = "rna_atac_directional",
    input_path = file.path(project_root, "outputs", "integration", "tables", "rna_atac_main_overlap.tsv"),
    direction_filter = "down"
  )
)

summary_rows <- list()

for (job in go_jobs) {
  if (!file.exists(job$input_path)) {
    warning("Skipping missing input: ", job$input_path)
    next
  }

  input_df <- read_tsv_base(job$input_path)
  if (!"gene_name" %in% names(input_df)) {
    warning("Skipping ", job$comparison_id, ": gene_name column missing.")
    next
  }

  gene_vec <- input_df$gene_name
  if (identical(job$job_type, "rna_directional")) {
    if (!"log2FoldChange" %in% names(input_df)) {
      warning("Skipping ", job$comparison_id, ": log2FoldChange column missing.")
      next
    }
    input_df$log2FoldChange <- as.numeric(input_df$log2FoldChange)
    if (identical(job$direction, "up")) {
      gene_vec <- input_df$gene_name[is.finite(input_df$log2FoldChange) & input_df$log2FoldChange > 0]
    } else {
      gene_vec <- input_df$gene_name[is.finite(input_df$log2FoldChange) & input_df$log2FoldChange < 0]
    }
  }
  if (identical(job$job_type, "rna_atac_directional")) {
    if (!all(c("rna_log2fc", "atac_log2fc") %in% names(input_df))) {
      warning("Skipping ", job$comparison_id, ": rna_log2fc or atac_log2fc column missing.")
      next
    }
    input_df$rna_log2fc  <- as.numeric(input_df$rna_log2fc)
    input_df$atac_log2fc <- as.numeric(input_df$atac_log2fc)
    if (identical(job$direction_filter, "up")) {
      gene_vec <- input_df$gene_name[
        is.finite(input_df$rna_log2fc) & input_df$rna_log2fc > 0 &
        is.finite(input_df$atac_log2fc) & input_df$atac_log2fc > 0
      ]
    } else {
      gene_vec <- input_df$gene_name[
        is.finite(input_df$rna_log2fc) & input_df$rna_log2fc < 0 &
        is.finite(input_df$atac_log2fc) & input_df$atac_log2fc < 0
      ]
    }
  }

  gene_vec <- sort(unique(gene_vec[!is.na(gene_vec) & nzchar(gene_vec)]))
  if (!length(gene_vec)) {
    warning("Skipping ", job$comparison_id, ": no genes after filtering.")
    next
  }

  message("Running topGO for ", job$comparison_id, " (", length(gene_vec), " genes).")
  topgo_res <- run_topgo_bp(gene_vec, job$label, node_size = 5)
  full_table <- topgo_res$table
  top40_table <- utils::head(full_table, 40)

  write_tsv_base(data.frame(gene_name = gene_vec, stringsAsFactors = FALSE),
                 file.path(table_dir, paste0(job$comparison_id, "_input_gene_set.tsv")))
  write_tsv_base(full_table, file.path(table_dir, paste0(job$comparison_id, "_topgo_bp_full.tsv")))
  write_tsv_base(top40_table, file.path(table_dir, paste0(job$comparison_id, "_topgo_bp_top40.tsv")))

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    comparison_id = job$comparison_id,
    job_type = job$job_type,
    input_rows = nrow(input_df),
    unique_genes = topgo_res$gene_count,
    background_genes = topgo_res$background_gene_count,
    enriched_terms = nrow(full_table),
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows)) {
  summary_df <- do.call(rbind, summary_rows)
  write_tsv_base(summary_df, file.path(table_dir, "summary_metrics.tsv"))
}

writeLines(capture.output(sessionInfo()), file.path(log_dir, "session_info.txt"))
message("Done. Results written to: ", result_dir)
