packages <- scan("envs/r_packages.txt", what = character(), quiet = TRUE)
cran_packages <- c(
  "tximport", "IHW", "aggregation", "ggplot2", "ggrepel", "dplyr", "tidyr",
  "readxl", "writexl", "data.table", "jsonlite"
)
bioc_packages <- c("DESeq2", "SummarizedExperiment", "clusterProfiler", "org.Hs.eg.db", "topGO", "GO.db")

cran_to_install <- intersect(packages, cran_packages)
bioc_to_install <- intersect(packages, bioc_packages)

install_if_missing <- function(pkg_vec, installer) {
  missing <- pkg_vec[!vapply(pkg_vec, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    installer(missing)
  }
}

install_if_missing(cran_to_install, function(pkgs) install.packages(pkgs))

if (length(bioc_to_install)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  install_if_missing(bioc_to_install, function(pkgs) BiocManager::install(pkgs, ask = FALSE, update = FALSE))
}
