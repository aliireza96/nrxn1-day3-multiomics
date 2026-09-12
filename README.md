# NRXN1 Day-3 Multi-Omics Manuscript Repository

This repository is the public code companion for the NRXN1alpha day-3 neural induction manuscript. It combines the assay-level analysis logic used to rebuild the manuscript results with the final figure-generation scripts used for the locked main and supplementary figures.

## Repository scope

This repo starts from `processed assay inputs`, not from raw FASTQ alignment pipelines and not from manuscript-only frozen figure tables. The goal is to show how the manuscript analyses and figures were generated while keeping the repository small enough to share and review.

Included here:
- canonical assay-level rebuild scripts in [`scripts/analysis`](scripts/analysis)
- final figure and supplementary-table builders in [`scripts/figures`](scripts/figures)
- metadata and manuscript-facing manifests in [`metadata`](metadata)
- placeholder input directories plus modality-specific setup notes in [`data/processed_inputs`](data/processed_inputs)

Not included here:
- raw sequencing files
- BAMs, BigWigs, or large peak/motif work directories
- assay-level processed input files
- manuscript draft `.docx` files and co-author packaging artifacts
- generated outputs committed to Git

Generated files should be written under [`outputs`](outputs) and shared as release assets rather than tracked in the repository history. Assay-level processed inputs should be distributed separately through a release archive, Zenodo record, institutional repository, or another linked location and unpacked into the documented paths before rerunning the workflow.

## Top-level layout

```text
config/                  curated feature definitions used across analyses
data/processed_inputs/   processed assay inputs for public rebuilds
docs/                    workflow, data, and figure-provenance notes
envs/                    dependency manifests
metadata/                sample sheet, comparisons, figure manifest, input manifest
outputs/                 runtime output directory (gitignored except .gitkeep)
scripts/analysis/        assay-level rebuild scripts
scripts/figures/         figure builders, assemblers, and supplementary-table export
```

## Quick start

1. Install the R packages listed in [`envs/r_packages.txt`](envs/r_packages.txt).
2. Install the Python packages listed in [`envs/python-figures-requirements.txt`](envs/python-figures-requirements.txt).
3. Populate [`data/processed_inputs`](data/processed_inputs) using:
   - [`metadata/processed_input_manifest.tsv`](metadata/processed_input_manifest.tsv)
   - the README file inside each modality folder
   - the header comments in the relevant scripts
4. Validate the public input layer:

```bash
Rscript scripts/analysis/00_validate_inputs.R
```

5. Run the analysis workflow in the order documented in [`docs/run_order.md`](docs/run_order.md).
6. Build figures and supplementary tables after the assay-level outputs are present.

## Main entrypoints

Analysis:
- [`00_validate_inputs.R`](scripts/analysis/00_validate_inputs.R)
- [`01_rna_expression.R`](scripts/analysis/01_rna_expression.R)
- [`02_splicing.R`](scripts/analysis/02_splicing.R)
- [`03_atac_accessibility.R`](scripts/analysis/03_atac_accessibility.R)
- [`04_h3k27me3.R`](scripts/analysis/04_h3k27me3.R)
- [`05_integration.R`](scripts/analysis/05_integration.R)
- [`06_prepare_figure_inputs.R`](scripts/analysis/06_prepare_figure_inputs.R)
- [`08_go_enrichment.R`](scripts/analysis/08_go_enrichment.R)
- [`10_splicing_go.R`](scripts/analysis/10_splicing_go.R)

Figures:
- [`07_plot_figures.R`](scripts/figures/07_plot_figures.R)
- [`09_density_pca_figure.R`](scripts/figures/09_density_pca_figure.R)
- [`09s_build_supp_qc_panels.R`](scripts/figures/09s_build_supp_qc_panels.R)
- [`assemble_complete_figures.py`](scripts/figures/assemble_complete_figures.py)
- [`build_supplementary_tables.py`](scripts/figures/build_supplementary_tables.py)

Each main analysis script documents:
- what the script does
- which inputs it expects
- which outputs it writes
- which earlier workflow steps it depends on
