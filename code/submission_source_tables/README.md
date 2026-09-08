# Submission source tables

These typed tables are the panel-specific outputs assembled by the upstream R analyses for the final round-2 figures. `manifest.json` fixes the submission sheet order, table positions, and descriptive/header rows. They are kept outside `source_data/` because they are derived analysis outputs, not public accession files.

`code/R/scripts/08_export_round2_supplementary_data.R` converts these tables into `Supplementary_Data_1_round2.xlsx`. Set `NPJAGING_SUPPLEMENTARY_OUTPUT` to write the workbook to a different path.

