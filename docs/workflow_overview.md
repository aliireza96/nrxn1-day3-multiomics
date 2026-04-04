# Workflow Overview

This repository is organized as a two-layer workflow.

## 1. Assay-level analysis layer

The scripts in [`scripts/analysis`](../scripts/analysis) rebuild the manuscript result tables from processed assay inputs.

Layer outputs:
- RNA differential expression and PCA-ready objects
- significant alternative splicing events from rMATS JCEC tables
- ATAC differential accessibility, promoter summaries, and motif-footprinting tables
- H3K27me3 differential-region and promoter-region summaries
- cross-modality overlap and concordance tables
- GO enrichment tables for manuscript figure panels

These scripts write to assay-specific folders under [`outputs`](../outputs).

## 2. Figure and supplementary-table layer

The scripts in [`scripts/figures`](../scripts/figures) consume assay-level outputs and build the manuscript figures and supplementary tables.

This layer includes:
- R figure drivers for the main RNA / QC / PCA panels
- standalone Python panel builders for locked manuscript panels
- complete-figure assembly into multi-panel PDFs
- supplementary-table workbook export

## Public workflow start point

The public workflow begins from `data/processed_inputs`, not from raw sequencing alignments and not from manuscript-only frozen figure tables. That keeps the repository shareable while still exposing the analysis logic used to generate the results.
