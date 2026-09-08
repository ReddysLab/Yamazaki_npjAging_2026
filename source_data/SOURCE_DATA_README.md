# Source data

This directory contains only public/source inputs and accession-linked download instructions. Generated analysis results are written to `code/results/`.

`SOURCE_MANIFEST.csv` records the dataset, accession/reference, role, bundled path, and SHA-256 hash. Large objects that cannot be committed to an ordinary GitHub repository are represented by a README and downloaded separately.

The final healthy-blood aging branch uses JenAge (GSE103232 and GSE75337), not Dutch 500FG/GSE134080. The directory name `GSE134080_or_JenAge` is retained only for compatibility with the established package layout.
