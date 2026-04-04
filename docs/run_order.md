# Run Order

## Full rebuild

Run the workflow from the repository root.

1. Validate inputs

```bash
Rscript scripts/analysis/00_validate_inputs.R
```

2. RNA differential expression

```bash
Rscript scripts/analysis/01_rna_expression.R --comparison_ids=rna_main_all_controls,rna_supp_female_only,rna_density_controls,rna_density_patient
```

3. Alternative splicing from rMATS JCEC tables

```bash
Rscript scripts/analysis/02_splicing.R
```

4. ATAC accessibility and TOBIAS summary tables

```bash
Rscript scripts/analysis/03_atac_accessibility.R
```

5. H3K27me3 differential-region processing

```bash
Rscript scripts/analysis/04_h3k27me3.R
```

6. Multi-omic integration

```bash
Rscript scripts/analysis/05_integration.R
```

7. GO enrichment

```bash
Rscript scripts/analysis/08_go_enrichment.R
Rscript scripts/analysis/10_splicing_go.R
```

8. Figure-input status snapshot

```bash
Rscript scripts/analysis/06_prepare_figure_inputs.R
```

9. Main and supplementary panel generation

```bash
Rscript scripts/figures/07_plot_figures.R
Rscript scripts/figures/09_density_pca_figure.R
Rscript scripts/figures/09s_build_supp_qc_panels.R
python3 scripts/figures/python/07k_plot_figure_1A_study_design.py
python3 scripts/figures/python/07f2_plot_figure_2F_chromatin_panel.py
python3 scripts/figures/python/07c_plot_figure_3C_go.py
python3 scripts/figures/python/07p2_plot_figure_3A_chip_annotation_diverging.py
python3 scripts/figures/python/07o_plot_figure_3B_rna_chip_scatter.py
python3 scripts/figures/python/07g2_plot_figure_4A_promoter_dars.py
python3 scripts/figures/python/07i_plot_figure_4B_rna_atac_scatter.py
python3 scripts/figures/python/07s_plot_figure_4C_rna_atac_go.py
python3 scripts/figures/python/07t_plot_figure_4D_tobias_volcano.py
python3 scripts/figures/python/07u_plot_figure_4E_tobias_mechanistic.py
python3 scripts/figures/python/07n_plot_figure_2E_splicing_go.py
python3 scripts/figures/python/plot_shared_manuscript_panels.py
python3 scripts/figures/python/plot_figure_5B_integration_heatmap.py
python3 scripts/figures/python/plot_figure_S3_qpcr_validation.py
```

10. Assemble complete figures

```bash
python3 scripts/figures/assemble_complete_figures.py
```

11. Build supplementary tables

```bash
python3 scripts/figures/build_supplementary_tables.py
```

## Figure-only rebuild

If assay-level outputs already exist under [`outputs`](../outputs), the shorter rebuild path is:

1. `Rscript scripts/analysis/06_prepare_figure_inputs.R`
2. run the figure scripts needed for the target panel(s)
3. `python3 scripts/figures/assemble_complete_figures.py`
4. `python3 scripts/figures/build_supplementary_tables.py`

## Notes

- Before running the workflow, unpack the documented assay-level inputs into `data/processed_inputs/`.
- `01_rna_expression.R` expects processed StringTie count directories under `data/processed_inputs/rna/stringtie_counts/`.
- `07u_plot_figure_4E_tobias_mechanistic.py` expects additional large footprint signal inputs under `data/processed_inputs/atac/footprinting/` if full Figure 4E regeneration is desired.
- Figure 5A and Figure 5C are maintained final manuscript assets rather than fully rebuilt plotting products in this public workflow.
