#!/usr/bin/env Rscript

# 10_splicing_go.R
# Run the topGO enrichment analysis for genes carrying significant splicing events.
#
# Required upstream inputs:
#   - outputs/splicing/tables/significant_events_jcec.tsv
#
# Main outputs:
#   - outputs/splicing_go/tables/
#   - outputs/figures/supp/splicing_go/
#
# Run after:
#   - 02_splicing.R

suppressPackageStartupMessages({
  req_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  if (!req_pkg("topGO")) stop("Package 'topGO' is required.")
  if (!req_pkg("org.Hs.eg.db")) stop("Package 'org.Hs.eg.db' is required.")
  if (!req_pkg("GO.db")) stop("Package 'GO.db' is required.")
  if (!req_pkg("ggplot2")) stop("Package 'ggplot2' is required.")
  library(topGO)
  library(GO.db)
  library(ggplot2)
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

save_plot <- function(p, path, width = 7.8, height = 5.6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = p, width = width, height = height)
  message("Saved: ", path)
}

theme_manuscript <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

wrap_text_vec <- function(x, width = 36) {
  vapply(x, function(val) paste(strwrap(as.character(val), width = width), collapse = "\n"), character(1))
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

run_topgo_bp <- function(gene_vec, node_size = 5) {
  gene_vec <- sort(unique(gene_vec[!is.na(gene_vec) & nzchar(gene_vec)]))
  if (!length(gene_vec)) stop("No genes available for topGO.")

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
  all_res$enrichment_score <- -log10(pmax(all_res$weight01Fisher_num, .Machine$double.xmin))
  all_res <- all_res[order(all_res$weight01Fisher_num, all_res$classicFisher_num, all_res$Term), , drop = FALSE]
  all_res
}

build_topgo_plot <- function(df, title_text, top_n = 10) {
  plot_df <- utils::head(df, top_n)
  plot_df$score <- as.numeric(plot_df$enrichment_score)
  plot_df$count_n <- as.numeric(plot_df$Significant)
  plot_df$term_label <- wrap_text_vec(plot_df$Term, width = 34)
  plot_df$count_label <- format(plot_df$count_n, big.mark = ",", trim = TRUE)
  plot_df$term_label <- factor(plot_df$term_label, levels = rev(unique(plot_df$term_label)))
  x_upper <- max(plot_df$score, na.rm = TRUE)
  label_offset <- x_upper * 0.045

  ggplot2::ggplot(plot_df, ggplot2::aes(x = score, y = term_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = score, yend = term_label),
      linewidth = 1.5,
      color = "#deebf7",
      lineend = "round"
    ) +
    ggplot2::geom_point(size = 3.1, shape = 21, stroke = 0.3, fill = "#08519c", color = "#08519c") +
    ggplot2::geom_label(
      ggplot2::aes(x = score + label_offset, label = count_label),
      hjust = 0,
      size = 2.9,
      label.size = 0.12,
      label.padding = grid::unit(0.08, "lines"),
      fill = grDevices::adjustcolor("white", alpha.f = 0.92),
      color = "#2c3e50"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_upper * 1.22),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = title_text,
      x = expression(-log[10]("weight01 Fisher")),
      y = NULL
    ) +
    theme_manuscript(base_size = 11) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8.8, lineheight = 0.9, margin = ggplot2::margin(r = 4)),
      plot.margin = ggplot2::margin(t = 7, r = 18, b = 7, l = 7)
    )
}

project_root <- find_project_root()
splicing_path <- file.path(project_root, "outputs", "splicing", "tables", "significant_events_jcec.tsv")
if (!file.exists(splicing_path)) stop("Missing splicing event table: ", splicing_path)

events_df <- read_tsv_base(splicing_path)
if (!"geneSymbol" %in% names(events_df)) stop("Splicing event table lacks geneSymbol column.")
gene_vec <- sort(unique(events_df$geneSymbol[!is.na(events_df$geneSymbol) & nzchar(events_df$geneSymbol)]))
if (!length(gene_vec)) stop("No valid spliced genes found.")

message("Running topGO for spliced genes (", length(gene_vec), " genes).")
topgo_df <- run_topgo_bp(gene_vec, node_size = 5)
top40_df <- utils::head(topgo_df, 40)

result_dir <- file.path(project_root, "outputs", "splicing_go")
table_dir <- file.path(result_dir, "tables")
fig_dir <- file.path(project_root, "outputs", "figures", "supp", "splicing_go")

write_tsv_base(data.frame(gene_name = gene_vec, stringsAsFactors = FALSE), file.path(table_dir, "splicing_input_gene_set.tsv"))
write_tsv_base(topgo_df, file.path(table_dir, "splicing_topgo_bp_full.tsv"))
write_tsv_base(top40_df, file.path(table_dir, "splicing_topgo_bp_top40.tsv"))
write_tsv_base(
  data.frame(
    gene_count = length(gene_vec),
    enriched_terms = nrow(topgo_df),
    stringsAsFactors = FALSE
  ),
  file.path(table_dir, "summary_metrics.tsv")
)

p <- build_topgo_plot(topgo_df, "Spliced genes topGO BP", top_n = 10)
save_plot(p, file.path(fig_dir, "Supplement_topgo_splicing_genes.pdf"), width = 7.8, height = 5.8)
save_plot(p, file.path(fig_dir, "Supplement_topgo_splicing_genes.png"), width = 7.8, height = 5.8)

message("Done. Splicing topGO analysis written to: ", result_dir)
