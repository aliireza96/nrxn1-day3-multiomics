# Data Availability

## Repository start point

This repository begins from processed assay inputs expected under [`data/processed_inputs`](../data/processed_inputs). The repository documents those inputs, but does not commit the assay-level files themselves.

Documented processed inputs:
- StringTie transcript abundance directories for the manuscript RNA sample set
- rMATS JCEC alternative-splicing event tables
- ATAC differential-accessibility and promoter-annotation tables
- ATAC motif-footprinting summary tables
- H3K27me3 differential-region and promoter-region tables
- qPCR expression and summary-statistics tables
- compact QC tables used for Supplementary Figure S1
- curated integration display-order tables

The exact expected paths are listed in [`metadata/processed_input_manifest.tsv`](../metadata/processed_input_manifest.tsv), and each modality directory contains a short README describing the required filenames and how they are used.

## Not included in Git

The repository does not track:
- raw FASTQ files
- aligned BAMs
- BigWig signal tracks
- full TOBIAS work directories
- full peak-calling work directories
- assay-level processed input files
- manuscript PDFs, PNGs, and supplementary workbooks as committed outputs

## Raw data and final assets

- Raw sequencing data should be provided through the manuscript data-accession statement or archive record when public.
- Assay-level processed inputs should be distributed separately from the Git repository through a release archive, Zenodo record, institutional repository, or another linked location.
- Final figure PDFs, supplementary figures, and supplementary tables are better distributed as release assets or through a Zenodo archive than committed into the code tree.

## Larger optional inputs

One panel uses optional larger inputs that are not bundled by default:
- Figure 4E footprint aggregation expects processed footprint bigWigs and selected BED site lists under `data/processed_inputs/atac/footprinting/`.

The rest of the public workflow is designed to run once the documented assay-level inputs have been unpacked into the expected locations.
