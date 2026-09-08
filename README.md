# npj Aging round-2 reproducibility package

This package contains the R workflow and public/source inputs for the final round-2 analyses. Generated files are written only under `code/results/`; bundled files in `source_data/` are not modified.

## Run

```sh
export NPJAGING_PROJECT_ROOT="$(pwd)"
export NPJAGING_SCHAUM_H5AD="/path/to/tabula_muris_senis.h5ad"
Rscript code/R/run_all.R
```

`NPJAGING_PROJECT_ROOT` is the only project-root setting. `NPJAGING_SCHAUM_H5AD` points to the downloaded Schaum et al. (GSE132040) single-cell object, which is not committed because of its size. Source accessions and bundled-file hashes are listed in `source_data/SOURCE_MANIFEST.csv`.

The tested R and package versions are recorded in `environment/R_VERSION.txt` and `environment/package_versions.tsv`.

## Final human sleep-restriction model

The Moller-Levet GSE39445 analysis uses

```r
expression ~ condition + sin(2*pi*phase/24) + cos(2*pi*phase/24) +
  condition:sin(2*pi*phase/24) + condition:cos(2*pi*phase/24)
```

with `limma::duplicateCorrelation(..., block = subject)` and participant identity as the repeated-measures block. The analyzed matrix contains 19,541 genes, 427 samples, and 26 participants.

The final human blood SR-down set requires all four criteria:

- participant-aware phase-adjusted log2FC < 0;
- phase-model BH FDR < 0.05;
- delta MESOR < 0;
- MESOR BH FDR < 0.05.

This definition yields 2,163 genes before panel-specific shared-universe restriction. It is used by the final downstream exporter for Fig. 3e-f, Fig. 4c/e-g, Supplementary Fig. 3b, Supplementary Fig. 4e-g, and Supplementary Fig. 5a-f.

Reactome networks display representative terms with nominal Reactome P < 0.05 and at least two overlapping genes. Full raw P values and BH FDR values are exported for every tested term; no HSF1/proteostasis term is forced into a network.

## Code order

- `01b_sleep_restriction_participant_aware_phase_model.R`: repeated-measures/cosinor gene and module models.
- `01c_propagate_participant_aware_moller_levet_dependencies.R`: Fig. 1 and Supplementary Fig. 1 source tables.
- `01d_update_participant_aware_fig1f_fig1g.R`: participant-aware Fig. 1f-g enrichment outputs.
- `01e_update_participant_aware_fig3de.R`: Schaum neuronal aging slopes and shared-universe membership.
- `04a_multiomics_rna_protein_aging_decline.R`: tissue RNA/protein aging-decline definitions; verifies brain n = 1,634.
- `01j_export_bothFDR_human_SRdown_matched_panels.R`: final dual-FDR downstream panels, full gene lists, Fisher tests, Reactome tables, and node/edge tables.
- `08_export_round2_supplementary_data.R`: deterministic export of the 56 panel-specific sheets plus seven reviewer-requested completeness/curation sheets to `Supplementary_Data_1_round2.xlsx`.
- `09_export_sf6e_exact_fisher_venn.R`: Supplementary Fig. 6e Venn panel using the exact Fisher result reported in sheet `SF6e`.
- `02`-`07`: retained source-audit modules for the acute SD/chronic SF, SACS, young-versus-old SD, multi-omics, spatial, and SEA-AD branches. They are not part of the participant-aware downstream rerun in `run_all.R`.

The final participant-aware outputs are written to `code/results/final_dual_fdr/`. The submission workbook is regenerated as `Supplementary_Data_1_round2.xlsx`; set `NPJAGING_SUPPLEMENTARY_OUTPUT` to choose another output path. Its typed panel tables and 63-sheet manifest are retained under `code/submission_source_tables/` as derived outputs, separate from the public inputs in `source_data/`. Reviewer-requested completeness sheets include the full participant-aware GSE39445 selection statistics, complete Hallmark/Reactome memberships, an initial-file mapping, Takasugi source metadata, and a curated tissue-gene multi-omics meta-analysis table. The 3.63-million-row replicate table remains compressed under `source_data/Takasugi_2024/` and is not duplicated in the submission workbook.

## Data notes

- GSE39445: Moller-Levet human crossover sleep restriction.
- GSE113754 and PRJNA757396: acute SD and chronic sleep-fragmentation inputs.
- GSE128770: young-versus-old acute SD response.
- GSE132040: Schaum et al./Tabula Muris Senis aging analysis.
- JenAge: healthy-blood aging analysis derived from GSE103232 and GSE75337. GSE134080/500FG is not used in the final Supplementary Fig. 4 analysis.
- Takasugi et al. 2024: mouse tissue RNA/protein aging data; the bundled long table is gzip-compressed.
- SEA-AD: public middle temporal gyrus mean-expression and AD-effect matrices; access instructions are provided in the dataset README.

## Reproducibility notes

Fisher tests use the panel-specific intersected tested-gene universe reported in each output table. Gene-level and pathway-level multiple-testing corrections use Benjamini-Hochberg FDR. SACS is defined as acute SD induction plus chronic SF attenuation plus aging decline. HSF1 activity is inferred from associated gene programs and was not directly measured.
