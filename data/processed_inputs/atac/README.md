## ATAC processed inputs

Place the processed ATAC-seq assay inputs needed by `scripts/analysis/03_atac_accessibility.R`
in this directory.

Expected files:

- `P_vs_ctrl_sig.txt`
- `P_vs_ctrl_sig_annot.txt`
- `bindetect_results.txt`
- `P_vs_ctrl14_sig.txt`
- `anno_P_vs_ctrl14_sig_df.txt`
- `anno_P_vs_ctrl14_sig_promoter_geneid.txt`
- `ctrl14_vs_allP_bindetect_results.csv`

Optional larger inputs:

- `footprinting/` for Figure 4E footprint aggregation inputs
- `macs2_peaks/` if you want to keep assay-specific peak sets alongside the
  public workflow
