# Figure Provenance

This document maps each locked manuscript figure to the public scripts and inputs used to generate it.

## Main figures

| Figure | Status | Public generation path |
|---|---|---|
| Figure 1A | fully scripted | `07k_plot_figure_1A_study_design.py` |
| Figure 1B | fully scripted | `07_plot_figures.R` |
| Figure 1C | fully scripted | `09_density_pca_figure.R` plus `assemble_complete_figures.py` |
| Figure 2A-D | fully scripted | `07_plot_figures.R` |
| Figure 2E | fully scripted | `10_splicing_go.R` plus `07n_plot_figure_2E_splicing_go.py` |
| Figure 2F | fully scripted | `07f2_plot_figure_2F_chromatin_panel.py` |
| Figure 3A | fully scripted | `07_plot_figures.R` plus `07p2_plot_figure_3A_chip_annotation_diverging.py` |
| Figure 3B | fully scripted | `05_integration.R` plus `07o_plot_figure_3B_rna_chip_scatter.py` |
| Figure 3C | fully scripted | `07_plot_figures.R` plus `07c_plot_figure_3C_go.py` |
| Figure 3D | fully scripted | `plot_shared_manuscript_panels.py` |
| Figure 4A | fully scripted | `03_atac_accessibility.R` plus `07g2_plot_figure_4A_promoter_dars.py` |
| Figure 4B | fully scripted | `05_integration.R` plus `07i_plot_figure_4B_rna_atac_scatter.py` |
| Figure 4C | fully scripted | `08_go_enrichment.R` plus `07s_plot_figure_4C_rna_atac_go.py` |
| Figure 4D | fully scripted | `07_plot_figures.R` plus `07t_plot_figure_4D_tobias_volcano.py` |
| Figure 4E | optional large inputs required | `07u_plot_figure_4E_tobias_mechanistic.py` with additional footprint signal inputs |
| Figure 5A | maintained final asset | final manuscript asset maintained outside the plotting scripts |
| Figure 5B | fully scripted | `05_integration.R` plus `plot_figure_5B_integration_heatmap.py` |
| Figure 5C | maintained final asset | final manuscript asset maintained outside the plotting scripts |

## Supplementary figures

| Figure | Status | Public generation path |
|---|---|---|
| Figure S1 | fully scripted | `09s_build_supp_qc_panels.R` plus `assemble_complete_figures.py` |
| Figure S2 | fully scripted | `07_plot_figures.R` plus `assemble_complete_figures.py` |
| Figure S3 | fully scripted | `plot_figure_S3_qpcr_validation.py` plus `assemble_complete_figures.py` |

## Notes

- The public figure workflow writes panel outputs under `outputs/figures/`.
- Figure 5A and Figure 5C are intentionally documented as maintained final assets because the locked manuscript versions are not rebuilt by a dedicated public plotting script in this repository.
- Figure 4E is script-generated, but reproducibility requires optional large processed footprint inputs not bundled by default.
