# Build the round2 submission source-data workbook from figure-source CSVs.
# Each worksheet name matches the corresponding figure-panel label.

suppressPackageStartupMessages({
  library(openxlsx)
  library(readr)
  library(dplyr)
  library(stringr)
})

rel_path <- function(...) file.path(...)
abs_path <- function(rel) file.path(project_root, rel)

read_table_safe <- function(rel) {
  path <- abs_path(rel)
  require_file(path, rel)
  if (file.info(path)$size == 0) {
    return(data.frame(note = paste("Source file is empty:", rel), stringsAsFactors = FALSE))
  }
  read_csv(path, show_col_types = FALSE, progress = FALSE)
}

filter_tissue <- function(dat, tissue) {
  if (!"tissue" %in% names(dat)) return(dat)
  dat %>% filter(tolower(.data$tissue) == tolower(tissue))
}

filter_rows <- function(dat, expr) {
  if (is.null(expr)) return(dat)
  expr(dat)
}

section <- function(rel, description, filter = NULL) {
  list(rel = rel, description = description, filter = filter)
}

add_sheet <- function(wb, sheet, sections) {
  addWorksheet(wb, sheet, gridLines = FALSE)
  row <- 1
  for (sec in sections) {
    dat <- filter_rows(read_table_safe(sec$rel), sec$filter)
    writeData(wb, sheet, x = data.frame(field = "Source", value = sec$rel), startRow = row, colNames = FALSE)
    addStyle(wb, sheet, createStyle(textDecoration = "bold", fgFill = "#D9EAF7"), rows = row, cols = 1:2, gridExpand = TRUE)
    row <- row + 1
    writeData(wb, sheet, x = data.frame(field = "Description", value = sec$description), startRow = row, colNames = FALSE)
    row <- row + 2
    if (nrow(dat) == 0) dat <- data.frame(note = "No rows after panel-specific filtering.", stringsAsFactors = FALSE)
    writeData(wb, sheet, x = dat, startRow = row, colNames = TRUE)
    addStyle(wb, sheet, createStyle(textDecoration = "bold", fgFill = "#EDEDED"), rows = row, cols = seq_len(ncol(dat)), gridExpand = TRUE)
    row <- row + nrow(dat) + 3
  }
  freezePane(wb, sheet, firstActiveRow = 4)
  setColWidths(wb, sheet, cols = 1:60, widths = 16)
  invisible(wb)
}

fix_xlsx_dimensions <- function(path) {
  tmp <- tempfile("xlsx_dim_fix_")
  dir.create(tmp)
  old_wd <- getwd()
  on.exit({
    setwd(old_wd)
    unlink(tmp, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  unzip(path, exdir = tmp)
  sheet_files <- list.files(file.path(tmp, "xl", "worksheets"), pattern = "^sheet[0-9]+\\.xml$", full.names = TRUE)
  col_to_num <- function(x) {
    chars <- utf8ToInt(x)
    Reduce(function(acc, ch) acc * 26 + ch - utf8ToInt("A") + 1, chars, init = 0)
  }
  num_to_col <- function(n) {
    out <- character()
    while (n > 0) {
      n <- n - 1
      out <- c(intToUtf8((n %% 26) + utf8ToInt("A")), out)
      n <- n %/% 26
    }
    paste(out, collapse = "")
  }
  for (sf in sheet_files) {
    txt <- paste(readLines(sf, warn = FALSE), collapse = "\n")
    refs <- unlist(regmatches(txt, gregexpr(' r="[A-Z]+[0-9]+"', txt, perl = TRUE)))
    refs <- gsub('^ r="|"$', "", refs)
    if (!length(refs)) next
    cols <- gsub("[0-9]", "", refs)
    rows <- as.integer(gsub("[A-Z]", "", refs))
    max_col <- max(vapply(cols, col_to_num, numeric(1)), na.rm = TRUE)
    max_row <- max(rows, na.rm = TRUE)
    dim_ref <- paste0("A1:", num_to_col(max_col), max_row)
    txt <- sub('<dimension ref="[^"]*"/>', paste0('<dimension ref="', dim_ref, '"/>'), txt, perl = TRUE)
    writeLines(txt, sf, useBytes = TRUE)
  }
  patched <- tempfile(fileext = ".xlsx")
  setwd(tmp)
  utils::zip(patched, files = list.files(tmp, all.files = TRUE, no.. = TRUE), flags = "-qr9X")
  file.copy(patched, path, overwrite = TRUE)
  invisible(path)
}

gh <- "source_materials/GitHub_ready"
stage <- "code/results/tables/participant_aware_staging"
ml <- "code/results/tables/moller_levet_participant_aware"
sacs <- "source_materials/SACS_reactome_update_outputs"
addv <- "source_materials/Aging_revision_upload_ready/additional_visuals"
prj <- "source_materials/Aging_revision_upload_ready/PRJNA757396_sleep_fragmentation"
mouse <- "source_materials/Aging_revision_upload_ready/mouse_cortex_sd_aging"

panel_order <- c(
  "F1b","F1d","F1e","F1f","F1g",
  "F2a","F2b","F2c","F2d",
  "F3b","F3c","F3d","F3e","F3f",
  "F4b","F4c","F4d","F4e","F4f","F4g",
  "F5b","F5c","F5d","F5e","F5f","F5g","F5h",
  "SF1a","SF1b","SF1c","SF1d","SF1e",
  "SF2c","SF2d","SF2e","SF2f",
  "SF3a","SF3b",
  "SF4b","SF4c","SF4d","SF4e","SF4f","SF4g",
  "SF5a","SF5b","SF5c","SF5d","SF5e","SF5f",
  "SF6a","SF6b","SF6d","SF6e","SF6f","SF6g"
)

specs <- list(
  F1b = list(section(rel_path(gh, "Figure1/Tables/Figure_1B_source_data.csv"), "Acute SD-responsive proteostasis genes plotted in Fig. 1b.")),
  F1d = list(
    section(rel_path(stage, "F1d_participant_aware_source_data.csv"), "Fig. 1d plotted time-course means/SEM for HSPH1 and HSPA5 with participant-aware inferential statistics."),
    section(rel_path(ml, "MollerLevet_participant_aware_HSF_HSP_gene_audit.csv"), "Participant-aware HSF/HSP gene-level model audit.")
  ),
  F1e = list(section(rel_path(stage, "F1e_SF1c_SF1d_SF1e_participant_aware_selected_reactome_modules.csv"), "Fig. 1e participant-aware predefined Reactome module MESOR-like effects, empirical P values, FDR values, amplitude and timing summaries.")),
  F1f = list(
    section(rel_path(gh, "Figure1/Tables/Figure_1F_nodes_source_data.csv"), "Fig. 1f Reactome network node source data."),
    section(rel_path(gh, "Figure1/Tables/Figure_1F_edges_source_data.csv"), "Fig. 1f Reactome network edge source data."),
    section(rel_path(gh, "Figure1/Tables/Figure_1F_CAMERA_reactome_full.csv"), "Full Reactome enrichment/CAMERA table underlying Fig. 1f.")
  ),
  F1g = list(section(rel_path(gh, "Figure1/Tables/Figure_1E_1F_gene_level_phase_adjusted_limma.csv"), "Gene-level phase-adjusted sleep-restriction results supporting Fig. 1 source-data panels.")),

  F2a = list(section(rel_path(gh, "Figure2/Tables/Figure_2A_source_data.csv"), "Fig. 2a sample/age source data.")),
  F2b = list(section(rel_path(gh, "Figure2/Tables/Figure_2B_source_data.csv"), "Fig. 2b gene-level aging direction source data.")),
  F2c = list(section(rel_path(gh, "Figure2/Tables/Figure_2C_source_data.csv"), "Fig. 2c gene trajectory cluster/source data.")),
  F2d = list(
    section(rel_path(gh, "Figure2/Tables/Figure_2D_source_data.csv"), "Fig. 2d representative gene trajectory plotting data."),
    section(rel_path(addv, "Figure2D_representative_gene_selection_statistics.csv"), "Selection statistics for representative Fig. 2d genes.")
  ),

  F3b = list(section(rel_path(gh, "Figure3/Tables/Figure_3B_source_data.csv"), "Fig. 3b aging-associated pathway absolute-change summary.")),
  F3c = list(section(rel_path(gh, "Figure3/Tables/Figure_3C_source_data.csv"), "Fig. 3c tissue-by-pathway aging heatmap source data.")),
  F3d = list(
    section(rel_path(addv, "Figure3D_neuronal_label_definitions.csv"), "Fig. 3d neuronal/cell-type label definitions."),
    section(rel_path(addv, "Figure3D_neuronal_label_gene_slope_audit.csv"), "Fig. 3d gene-by-cell-type aging slopes.")
  ),
  F3e = list(
    section(rel_path(gh, "Figure3/Tables/Figure_3E_source_data.csv"), "Fig. 3e SR-aging overlap Venn/Fisher source data."),
    section(rel_path(gh, "Figure3/Tables/Figure_3F_aligned_decline_genes.csv"), "Full aligned-decline/convergent gene list supporting the Fig. 3e overlap.")
  ),
  F3f = list(section(rel_path(gh, "Figure3/Tables/Figure_3F_source_data.csv"), "Fig. 3f Reactome ORA/network source data for excitatory-neuron SR-aging convergence.")),

  F4b = list(section(rel_path(gh, "Figure4/Tables/Figure_4B_stats.csv"), "Fig. 4b RNA-protein directional coupling Fisher statistics.")),
  F4c = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_gene_sets.csv"), "Fig. 4c multi-omics SR-aging convergence gene-set source table.")),
  F4d = list(section(rel_path(gh, "Figure4/Tables/Figure_4C_source_data.csv"), "Fig. 4d representative multi-omics aging trajectory source data.")),
  F4e = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Fig. 4e brain Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Brain"))),
  F4f = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Fig. 4f liver Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Liver"))),
  F4g = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_source_data.csv"), "Fig. 4g recurrent Reactome enrichment heatmap source data.")),

  F5b = list(
    section(rel_path(gh, "Figure5/Tables/Figure_5B_young_reference_spots.csv"), "Fig. 5b young reference spatial spots."),
    section(rel_path(gh, "Figure5/Tables/Figure_5B_middle_spots.csv"), "Fig. 5b middle-age spatial spots."),
    section(rel_path(gh, "Figure5/Tables/Figure_5B_old_spots.csv"), "Fig. 5b old spatial spots."),
    section(rel_path(gh, "Figure5/Tables/Figure_5B_region_age_summary_detail.csv"), "Fig. 5b regional age summary.")
  ),
  F5c = list(section(rel_path(gh, "Figure5/Tables/Figure_5C_section_summary.csv"), "Fig. 5c section-level SACS score summary.")),
  F5d = list(
    section(rel_path(gh, "Figure5/Tables/Figure_5D_anatomy_age_summary.csv"), "Fig. 5d anatomy-by-age SACS score summary."),
    section(rel_path(sacs, "SACS_Fig5D_region_stats_raw_score.csv"), "Fig. 5d statistical tests by region.")
  ),
  F5e = list(section(rel_path(gh, "Figure5/Tables/Figure_5E_source_data.csv"), "Fig. 5e SEA-AD donor age/disease-group source data.")),
  F5f = list(section(rel_path(gh, "Figure5/Tables/Figure_5F_source_data.csv"), "Fig. 5f baseline SACS enrichment across SEA-AD cell types.")),
  F5g = list(
    section(rel_path(gh, "Figure5/Tables/Figure_5G_source_data.csv"), "Fig. 5g baseline enrichment versus AD-associated SACS shift."),
    section(rel_path(gh, "Figure5/Tables/Figure_5G_correlation_summary.csv"), "Fig. 5g Spearman correlation summary.")
  ),
  F5h = list(
    section(rel_path(gh, "Figure5/Tables/Figure_5H_source_data.csv"), "Fig. 5h projected AD-vulnerability and dominant cell-type spatial source data."),
    section(rel_path(gh, "Figure5/Tables/Figure_5H_AD_weights_by_celltype.csv"), "Cell-type weights used for Fig. 5h projection.")
  ),

  SF1a = list(section(rel_path(stage, "SF1a_participant_aware_source_data.csv"), "Supplementary Fig. 1a descriptive time-course source data with participant-aware inferential statistics.")),
  SF1b = list(section(rel_path(stage, "SF1b_participant_aware_reactome_wide_MESOR_ranking.csv"), "Supplementary Fig. 1b participant-aware Reactome-wide MESOR-like module ranking.")),
  SF1c = list(section(rel_path(stage, "F1e_SF1c_SF1d_SF1e_participant_aware_selected_reactome_modules.csv"), "Supplementary Fig. 1c rhythmic-amplitude component for predefined modules.")),
  SF1d = list(section(rel_path(stage, "F1e_SF1c_SF1d_SF1e_participant_aware_selected_reactome_modules.csv"), "Supplementary Fig. 1d acrophase/timing component for predefined modules.")),
  SF1e = list(
    section(rel_path(stage, "F1e_SF1c_SF1d_SF1e_participant_aware_selected_reactome_modules.csv"), "Supplementary Fig. 1e observed participant-aware module MESOR effects."),
    section(rel_path(stage, "SF1e_participant_aware_reactome_random_nulls.csv"), "Supplementary Fig. 1e matched random-set null distributions.")
  ),

  SF2c = list(section(rel_path(prj, "PRJNA757396_73_GSE113754_HSF_HSP_acute_chronic_transition_examples.csv"), "Supplementary Fig. 2c HSF/HSP acute-SD versus chronic-SF example genes.")),
  SF2d = list(section(rel_path(sacs, "SACS_acuteSD_chronicSF_attenuation_all_gene_scores.csv"), "Supplementary Fig. 2d acute SD and chronic SF effect-size/attenuation score gene-level table.")),
  SF2e = list(section(rel_path(sacs, "SACS_acuteSD_up_x_positive_attenuation_venn_stats.csv"), "Supplementary Fig. 2e acute SD-up versus positive chronic-attenuation Venn/Fisher statistics.")),
  SF2f = list(
    section(rel_path(sacs, "SACS_acuteSD_chronicSF_attenuation_statistics.csv"), "Supplementary Fig. 2f attenuation-score statistics."),
    section(rel_path(sacs, "SACS_triple_convergence_Reactome_network_fdr_plus_hsf_labels_labeled_terms.csv"), "Supplementary Fig. 2f Reactome enrichment/network terms.")
  ),

  SF3a = list(
    section(rel_path(addv, "Figure3D_neuronal_binning_group_definitions.csv"), "Supplementary Fig. 3a alternative neuronal subtype group definitions."),
    section(rel_path(addv, "Figure3D_neuronal_binning_gene_slopes.csv"), "Supplementary Fig. 3a gene-by-neuronal-subtype aging slopes.")
  ),
  SF3b = list(section(rel_path(gh, "Figure3/Tables/Figure_3F_source_data.csv"), "Supplementary Fig. 3b inhibitory-neuron convergence Reactome network/ORA source data.")),

  SF4b = list(section(rel_path(stage, "SF4b_JenAge_Reactome_module_ranking.csv"), "Supplementary Fig. 4b JenAge human-blood Reactome module aging ranking.")),
  SF4c = list(section(rel_path(stage, "SF4c_JenAge_representative_Reactome_modules.csv"), "Supplementary Fig. 4c representative JenAge Reactome module aging effects.")),
  SF4d = list(section(rel_path(stage, "SF4d_JenAge_proteostasis_gene_aging_direction.csv"), "Supplementary Fig. 4d JenAge HSF/proteostasis gene aging-direction data.")),
  SF4e = list(
    section(
      rel_path(stage, "SF4e_participant_aware_SRdown_x_AgingDown_fisher_summary.csv"),
      "Supplementary Fig. 4e participant-aware SR-down x JenAge aging-down Venn/Fisher statistics.",
      function(dat) dat %>%
        filter(.data$dataset == "JenAge", .data$sr_threshold == "sr_nominalP05", .data$aging_threshold == "age_direction_down")
    ),
    section(rel_path(stage, "SF4e_participant_aware_JenAge_SRnominalP05_AgeDirection_overlap_genes.csv"), "Full SR-aging convergent gene list plotted in Supplementary Fig. 4e.")
  ),
  SF4f = list(
    section(rel_path(stage, "SF4f_participant_aware_Reactome_network_nodes.csv"), "Supplementary Fig. 4f Reactome network nodes."),
    section(rel_path(stage, "SF4f_participant_aware_Reactome_network_edges.csv"), "Supplementary Fig. 4f Reactome network edges."),
    section(rel_path(stage, "SF4f_SF4g_participant_aware_JenAge_SRnominalP05_AgeDirection_Reactome_ORA.csv"), "Supplementary Fig. 4f full Reactome ORA table.")
  ),
  SF4g = list(
    section(rel_path(stage, "SF4g_participant_aware_Reactome_family_summary.csv"), "Supplementary Fig. 4g Reactome family summary."),
    section(rel_path(stage, "SF4g_participant_aware_Reactome_family_gene_membership.csv"), "Supplementary Fig. 4g Reactome family gene membership.")
  ),

  SF5a = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Supplementary Fig. 5a lung Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Lung"))),
  SF5b = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Supplementary Fig. 5b aorta Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Aorta"))),
  SF5c = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Supplementary Fig. 5c heart Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Heart"))),
  SF5d = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Supplementary Fig. 5d muscle Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Muscle"))),
  SF5e = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Supplementary Fig. 5e kidney Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Kidney"))),
  SF5f = list(section(rel_path(gh, "Figure4/Tables/Figure_4D_enrichment_nonredundant.csv"), "Supplementary Fig. 5f skin Reactome network/enrichment source rows.", function(dat) filter_tissue(dat, "Skin"))),

  SF6a = list(
    section(rel_path(sacs, "SACS_triple_convergence_network_fdr_nominal_summary.csv"), "Supplementary Fig. 6a SACS definition overlap/Fisher summary."),
    section(rel_path(prj, "PRJNA757396_76_GSE113754_triple_convergence_chronic_not_significantly_up_genes.csv"), "Full SACS gene list: acute SD induction + chronic SF attenuation/not-up + aging decline.")
  ),
  SF6b = list(
    section(rel_path(sacs, "SACS_Reactome_network_FDR_only_nodes.csv"), "Supplementary Fig. 6b SACS Reactome network nodes (FDR < 0.05)."),
    section(rel_path(sacs, "SACS_Reactome_network_FDR_only_edges.csv"), "Supplementary Fig. 6b SACS Reactome network edges."),
    section(rel_path(sacs, "SACS_Reactome_network_FDR_only_summary.csv"), "Supplementary Fig. 6b SACS Reactome ORA summary.")
  ),
  SF6d = list(section(rel_path(mouse, "mouse_19_GSE128770_aging_attenuated_SD_responding_by_duration.csv"), "Supplementary Fig. 6d young-induced, old-induced, and aging-attenuated acute SD response gene counts by SD duration.")),
  SF6e = list(
    section(rel_path(mouse, "mouse_58_brain_SR_aging_overlap_aging_attenuated_SD_response_fisher.csv"), "Supplementary Fig. 6e aging-attenuated SD response x SACS Fisher statistics."),
    section(rel_path(mouse, "mouse_59_brain_SR_aging_overlap_aging_attenuated_SD_response_genes.csv"), "Supplementary Fig. 6e overlap gene list.")
  ),
  SF6f = list(section(rel_path(mouse, "mouse_22_GSE128770_aging_attenuated_SD_responding_HSF_HSP_genes_for_bar.csv"), "Supplementary Fig. 6f representative aging-attenuated SD-response gene trajectories.")),
  SF6g = list(
    section(rel_path(mouse, "mouse_28_GSE128770_aging_attenuated_SD_Reactome_network_nodes.csv"), "Supplementary Fig. 6g Reactome network nodes."),
    section(rel_path(mouse, "mouse_29_GSE128770_aging_attenuated_SD_Reactome_network_edges.csv"), "Supplementary Fig. 6g Reactome network edges."),
    section(rel_path(mouse, "mouse_26_GSE128770_aging_attenuated_SD_Reactome_enrichment.csv"), "Supplementary Fig. 6g Reactome ORA table.")
  )
)

missing_specs <- setdiff(panel_order, names(specs))
if (length(missing_specs)) stop("Missing workbook specs for sheets: ", paste(missing_specs, collapse = ", "), call. = FALSE)

wb <- createWorkbook(creator = "Round2 source-data exporter")
for (sheet in panel_order) add_sheet(wb, sheet, specs[[sheet]])

root_out <- file.path(project_root, "Supplementary_Data_1_round2.xlsx")
compat_out <- file.path(project_root, "Supplementary_Data_1_full_analysis_outputs.xlsx")
amendment_out <- file.path(dirname(project_root), "Supplementary_Data_1_round2.xlsx")

saveWorkbook(wb, root_out, overwrite = TRUE)
saveWorkbook(wb, compat_out, overwrite = TRUE)
if (dir.exists(dirname(amendment_out))) saveWorkbook(wb, amendment_out, overwrite = TRUE)
fix_xlsx_dimensions(root_out)
fix_xlsx_dimensions(compat_out)
if (file.exists(amendment_out)) fix_xlsx_dimensions(amendment_out)

audit <- data.frame(
  sheet = panel_order,
  section_count = vapply(panel_order, function(x) length(specs[[x]]), integer(1)),
  stringsAsFactors = FALSE
)
write_csv(audit, file.path(paths$tables, "round2_supplementary_data_sheet_audit.csv"))
message("Wrote round2 supplementary source-data workbook: ", amendment_out)
