## RNA processed inputs

This directory is the public landing point for the processed RNA assay inputs used by
`scripts/analysis/01_rna_expression.R`.

Expected layout:

- `stringtie_counts/<sample_id>/e_data.ctab`
- `stringtie_counts/<sample_id>/i_data.ctab`
- `stringtie_counts/<sample_id>/t_data.ctab`

The full `stringtie_counts` tree is not committed to Git because it is too large for a
lightweight manuscript code repository. Distribute it through a release archive, Zenodo
record, or another linked data object and unpack it here before rerunning the RNA analysis.
