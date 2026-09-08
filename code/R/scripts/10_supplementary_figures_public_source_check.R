# Supplementary figure public/source-data reproducibility audit.
# This module maps the actual Supplementary Figs. 1-6 in the npj Aging
# resubmission figure deck to public/source-data scripts and copied result files.

supp_figures <- data.frame(
  supplementary_figure = c(
    "Supplementary Fig. 1",
    "Supplementary Fig. 2",
    "Supplementary Fig. 3",
    "Supplementary Fig. 4",
    "Supplementary Fig. 5",
    "Supplementary Fig. 6"
  ),
  title = c(
    "Phase-resolved and cosinor-style analyses of HSF1-associated genes in the human sleep-restriction dataset",
    "HSF1-associated acute sleep-deprivation responses exhibit preferential chronic attenuation in mouse brain",
    "Neuronal subtype-specific aging trajectories and inhibitory-neuron convergence network",
    "Human sleep-restriction and healthy-blood aging signatures share transcriptional features enriched for HSF1-associated proteostasis pathways",
    "Multi-omics convergence of aging and sleep restriction reveals systemic HSF1-associated proteostasis decline",
    "Aging attenuates acute SD-responsive transcriptional programs in mouse brain"
  ),
  public_source_script = c(
    "scripts/36_moller_levet_mesor_confirmed_suppression_definition.R",
    "scripts/10_PRJNA757396_DESeq2_sensitivity.R; scripts/17_acute_SD_Hsf1_status_GSE211088_GSE113754.R; scripts/18/44/59-64 PRJNA/GSE113754 transition checks",
    "GitHub_ready/Figure3 source tables plus neuronal-label audit outputs copied under source_materials",
    "human blood sleep-restriction/aging overlap reviewer-response scripts; Moller-Levet GSE39445 public data; healthy blood aging source outputs copied under source_materials",
    "GitHub_ready/Figure4 multi-omics tables plus RNA/protein concordance and tissue Reactome outputs",
    "scripts/05_age_attenuated_sd_response.R plus mouse cortex SD-aging public-source scripts for GSE128770"
  ),
  source_accessions_or_resources = c(
    "GSE39445; Reactome/Hallmark gene-set resources; copied derived source tables",
    "GSE113754; PRJNA757396; source mouse aging-reference outputs; Reactome resources",
    "Public neuronal aging/source resources used for Figure 3; copied Figure3 source tables and neuronal-label audit tables",
    "GSE39445; healthy human blood aging source outputs used in reviewer-response package; Reactome resources",
    "Public aging multi-omics/source resources used for Figure 4; copied Figure4 source tables; Reactome resources",
    "GSE128770; mouse acute-to-chronic-to-aging convergence source outputs; Reactome resources"
  ),
  public_source_status = c(
    "Public accession/source resources confirmed for plotted inputs",
    "Public accession/source resources confirmed for plotted inputs",
    "Public/source resources represented by copied source-data tables used for plotted inputs",
    "Public/source resources confirmed for plotted inputs; healthy-blood aging source outputs are included in the final source manifest",
    "Public/source resources represented by copied source-data tables used for plotted inputs",
    "Public accession/source resources confirmed for plotted inputs"
  ),
  main_thresholds = c(
    "phase-aware/cosinor-style SR-control MESOR, amplitude, and acrophase summaries; matched-size random gene-set nulls",
    "acute SD induction FDR < 0.05; chronic SF transition-down/attenuation; aging slope < 0; Reactome ORA/FDR where indicated",
    "neuronal subtype definitions and label/gene-slope audits as exported in source tables; statistics as reported in Figure3 source data",
    "SR-down nominal P < 0.05 or phase-aware thresholds as stated; aging-down negative age direction; Reactome/Fisher using intersected expressed background",
    "RNA/protein concordance and tissue-specific Reactome analyses using exported Figure4 source-data thresholds and background universes",
    "young SD induction log2FC > 0 and FDR < 0.05; old attenuation by young-old contrast/trajectory; overlap with SACS/acute-chronic-aging set; Fisher/Reactome tests"
  ),
  stringsAsFactors = FALSE
)

rel_paths <- list(
  c(
    "Aging_revision_upload_ready/additional_visuals/Additional_Figure_1_MollerLevet_cosinor_summary.png",
    "Aging_revision_upload_ready/additional_visuals/MollerLevet_predefined_HSF_proteostasis_modules_MESOR_amp_phase_random.xlsx",
    "Aging_revision_upload_ready/additional_visuals/MollerLevet_Hallmark_modules_MESOR_amp_phase_summary.csv",
    "Aging_revision_upload_ready/additional_visuals/MollerLevet_phase_adjusted_plus_cosinor_gene_level_suppression_definitions.csv"
  ),
  c(
    "Aging_revision_upload_ready/PRJNA757396_sleep_fragmentation/PRJNA757396_Figure_GSE113754_HSF_HSP_acute_chronic_transition_examples.png",
    "Aging_revision_upload_ready/PRJNA757396_sleep_fragmentation/PRJNA757396_Figure_GSE113754_acuteFDRup_transition_counts.png",
    "Aging_revision_upload_ready/PRJNA757396_sleep_fragmentation/PRJNA757396_59_GSE113754_acuteFDRup_to_chronic_transition_all_genes.csv",
    "Aging_revision_upload_ready/PRJNA757396_sleep_fragmentation/PRJNA757396_60_GSE113754_acuteFDRup_transition_summary.csv",
    "Aging_revision_upload_ready/PRJNA757396_sleep_fragmentation/PRJNA757396_63_GSE113754_acuteFDRup_transition_down_FDR05_chronic_down_aging_down_genes.csv",
    "SACS_reactome_update_outputs/SACS_triple_convergence_Reactome_network_fdr_plus_hsf_labels_labeled_terms.csv"
  ),
  c(
    "Aging_revision_upload_ready/additional_visuals/Figure3D_neuronal_label_definitions.csv",
    "Aging_revision_upload_ready/additional_visuals/Figure3D_neuronal_label_gene_slope_audit.csv",
    "Aging_revision_upload_ready/additional_visuals/Figure3D_neuronal_label_summary_by_gene_class.csv",
    "Aging_revision_upload_ready/additional_visuals/Figure3D_neuronal_binning_group_definitions.csv",
    "GitHub_ready/Figure3/Tables/Figure_3D_source_data.csv",
    "GitHub_ready/Figure3/Tables/Figure_3E_source_data.csv",
    "GitHub_ready/Figure3/Tables/Figure_3F_source_data.csv"
  ),
  c(
    "Aging_revision_upload_ready/additional_visuals/MollerLevet_Reactome_heat_circadian_IL_ROS_AGE_modules.png",
    "Aging_revision_upload_ready/additional_visuals/MollerLevet_Reactome_heat_circadian_IL_ROS_AGE_modules.xlsx",
    "Aging_revision_upload_ready/additional_visuals/MollerLevet_SR_down_nominal_phase_adjusted_and_MESOR_down_genes_full.csv",
    "Aging_revision_upload_ready/additional_visuals/MollerLevet_Reactome_heat_circadian_IL_ROS_AGE_gene_classification_summary.csv"
  ),
  c(
    "GitHub_ready/Figure4/Tables/Figure_4A_multiomics_long_table.csv",
    "GitHub_ready/Figure4/Tables/Figure_4_multiomics_long.csv",
    "GitHub_ready/Figure4/Tables/Figure_4B_stats.csv",
    "GitHub_ready/Figure4/Tables/Figure_4C_source_data.csv",
    "GitHub_ready/Figure4/Tables/Figure_4D_source_data.csv",
    "GitHub_ready/Figure4/Tables/Figure_4D_enrichment.csv",
    "GitHub_ready/Figure4/Tables/Figure_4D_enrichment_all.csv",
    "Aging_revision_upload_ready/additional_visuals/SR_aging_RNA_protein_concordance_original_scale_source.csv",
    "Aging_revision_upload_ready/additional_visuals/SR_aging_HSF_RNAprotein_enrichment_main_term_summary.csv"
  ),
  c(
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/Mouse_Figure_2_age_attenuated_SD_response.png",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/Mouse_Figure_8_GSE128770_aging_attenuated_SD_responding_top_genes.png",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/Mouse_Figure_9_GSE128770_aging_attenuated_SD_responding_HSF_HSP_genes.png",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/Mouse_Figure_10_GSE128770_aging_attenuated_overlap_triple_venn.png",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/Mouse_Figure_12_GSE128770_aging_attenuated_SD_Reactome_network.png",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/mouse_19_GSE128770_aging_attenuated_SD_responding_by_duration.csv",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/mouse_26_GSE128770_aging_attenuated_SD_Reactome_enrichment.csv",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/mouse_28_GSE128770_aging_attenuated_SD_Reactome_network_nodes.csv",
    "Aging_revision_upload_ready/mouse_cortex_sd_aging/mouse_29_GSE128770_aging_attenuated_SD_Reactome_network_edges.csv"
  )
)

rows <- do.call(rbind, lapply(seq_len(nrow(supp_figures)), function(i) {
  full_paths <- file.path(paths$source_materials, rel_paths[[i]])
  data.frame(
    supplementary_figure = supp_figures$supplementary_figure[i],
    title = supp_figures$title[i],
    public_source_script = supp_figures$public_source_script[i],
    source_accessions_or_resources = supp_figures$source_accessions_or_resources[i],
    public_source_status = supp_figures$public_source_status[i],
    main_thresholds = supp_figures$main_thresholds[i],
    relative_source_path = rel_paths[[i]],
    present = file.exists(full_paths),
    bytes = ifelse(file.exists(full_paths), file.info(full_paths)$size, NA_real_),
    stringsAsFactors = FALSE
  )
}))

write.csv(rows, file.path(paths$tables, "supplementary_figure_public_source_manifest.csv"), row.names = FALSE)

summary <- aggregate(present ~ supplementary_figure + title + public_source_script + source_accessions_or_resources + public_source_status + main_thresholds,
                     data = rows, FUN = function(x) all(x))
names(summary)[names(summary) == "present"] <- "all_required_sources_present"
write.csv(summary, file.path(paths$tables, "supplementary_figure_public_source_summary.csv"), row.names = FALSE)

if (!all(rows$present)) {
  missing <- rows[!rows$present, c("supplementary_figure", "relative_source_path")]
  write.csv(missing, file.path(paths$tables, "supplementary_figure_missing_sources.csv"), row.names = FALSE)
  stop("Missing supplementary figure source files. See supplementary_figure_missing_sources.csv.", call. = FALSE)
}
