# Correct Supplementary Fig. 4 source data so panels b-g use one healthy-blood
# aging dataset: JenAge blood RNA-seq (GSE103232 + GSE75337).

if (!exists("paths")) {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
  project_root_guess <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
  Sys.setenv(NPJAGING_PROJECT_ROOT = Sys.getenv("NPJAGING_PROJECT_ROOT", unset = project_root_guess))
  source(file.path(dirname(script_dir), "config.R"))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(msigdbr)
  library(openxlsx)
  library(purrr)
})

set.seed(20260811)

stage_dir <- file.path(paths$tables, "participant_aware_staging")
healthy_dir <- paths$human_healthy_blood
workbook_path <- file.path(project_root, "Supplementary_Data_1_full_analysis_outputs.xlsx")
backup_path <- file.path(project_root, "Supplementary_Data_1_full_analysis_outputs.pre_sf4_jenage_only_backup.xlsx")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

require_file(workbook_path, "Supplementary Data workbook")
require_file(file.path(healthy_dir, "JenAge_01_healthy_blood_age_log2RPKM_results.csv"), "JenAge aging gene-level results")
require_file(file.path(stage_dir, "SF4_participant_aware_human_blood_SR_aging_merged_gene_table.csv"), "participant-aware SR x aging merged table")

std_gene <- function(x) toupper(str_trim(as.character(x)))

clean_term <- function(x) {
  x %>%
    str_remove("^REACTOME_") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

short_p <- function(x) {
  x <- as.numeric(x)
  ifelse(is.na(x), NA_character_, ifelse(x <= 0, "<1e-300", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))))
}

jen_age <- read_csv(file.path(healthy_dir, "JenAge_01_healthy_blood_age_log2RPKM_results.csv"), show_col_types = FALSE) %>%
  transmute(
    gene = std_gene(gene),
    age_logFC_per_decade = as.numeric(jenage_blood_age_logFC_per_decade),
    age_t = as.numeric(jenage_blood_age_t),
    age_p = as.numeric(jenage_blood_age_p),
    age_FDR = as.numeric(jenage_blood_age_fdr),
    age_ave_log2_rpkm = as.numeric(jenage_blood_age_ave_log2_rpkm)
  ) %>%
  filter(gene != "", is.finite(age_logFC_per_decade)) %>%
  distinct(gene, .keep_all = TRUE)

reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  transmute(gs_name, gene = std_gene(gene_symbol)) %>%
  filter(gene %in% jen_age$gene) %>%
  distinct()

module_gene_tbl <- reactome %>%
  inner_join(jen_age, by = "gene")

module_base <- module_gene_tbl %>%
  group_by(gs_name) %>%
  summarise(
    term_label = clean_term(first(gs_name)),
    n_genes = n_distinct(gene),
    median_age_effect = median(age_logFC_per_decade, na.rm = TRUE),
    genes_tested = paste(sort(unique(gene)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 3)

universe <- jen_age %>% filter(is.finite(age_logFC_per_decade)) %>% pull(gene)
n_perm <- 5000
age_by_gene <- jen_age$age_logFC_per_decade[match(universe, jen_age$gene)]
null_by_size <- map_dfr(sort(unique(module_base$n_genes)), function(n) {
  tibble(
    n_genes = n,
    iteration = seq_len(n_perm),
    null_median_age_effect = replicate(n_perm, median(sample(age_by_gene, n, replace = FALSE), na.rm = TRUE))
  )
})

null_tbl <- module_base %>%
  select(gs_name, n_genes) %>%
  left_join(null_by_size, by = "n_genes")
null_summary <- null_by_size %>%
  group_by(n_genes) %>%
  summarise(
    null_q025 = quantile(null_median_age_effect, 0.025, na.rm = TRUE),
    null_q975 = quantile(null_median_age_effect, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
null_list <- split(null_by_size$null_median_age_effect, null_by_size$n_genes)
empirical_tail <- function(effect, n, tail) {
  vals <- null_list[[as.character(n)]]
  if (tail == "negative") {
    max(mean(vals <= effect), 1 / n_perm)
  } else {
    max(mean(vals >= effect), 1 / n_perm)
  }
}

sf4b <- module_base %>%
  left_join(null_summary, by = "n_genes") %>%
  rowwise() %>%
  mutate(
    empirical_p_more_negative = empirical_tail(median_age_effect, n_genes, "negative"),
    empirical_p_more_positive = empirical_tail(median_age_effect, n_genes, "positive"),
    direction = if_else(median_age_effect < 0, "negative", "positive"),
    empirical_p_directional = if_else(direction == "negative", empirical_p_more_negative, empirical_p_more_positive)
  ) %>%
  ungroup() %>%
  mutate(
    fdr_directional = p.adjust(empirical_p_directional, method = "BH"),
    hsf_proteostasis = str_detect(gs_name, regex("HSF1|HEAT_SHOCK|CHAPERONE|UNFOLDED_PROTEIN|PROTEIN_FOLDING|PROTEOSTASIS|HSP90", ignore_case = TRUE)),
    dataset = "JenAge",
    dataset_source = "JenAge blood RNA-seq; GSE103232 and GSE75337",
    rank = rank(median_age_effect, ties.method = "first")
  ) %>%
  arrange(median_age_effect) %>%
  mutate(rank = row_number()) %>%
  select(dataset, dataset_source, gs_name, term_label, n_genes, median_age_effect,
         null_q025, null_q975, empirical_p_directional, fdr_directional,
         hsf_proteostasis, genes_tested, rank)

program_terms <- tibble(
  program = c(
    "B-cell/Fc/complement module",
    "Signal-transduction module",
    "HSF1 heat-shock regulation",
    "Cellular heat-stress response",
    "Cell-death module",
    "Interleukin/inflammatory module"
  ),
  gs_name = c(
    "REACTOME_CD22_MEDIATED_BCR_REGULATION",
    "REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION",
    "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",
    "REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS",
    "REACTOME_PROGRAMMED_CELL_DEATH",
    "REACTOME_NON_CANONICAL_INFLAMMASOME_ACTIVATION"
  ),
  family = c(
    "Broad representative",
    "Broad representative",
    "HSF1/heat-shock arm",
    "HSF1/heat-shock arm",
    "Broad representative",
    "Broad representative"
  )
)

sf4c <- program_terms %>%
  left_join(sf4b, by = "gs_name") %>%
  mutate(
    row_label = paste0(program, " (n=", n_genes, ")"),
    sig_label = case_when(
      is.na(fdr_directional) ~ NA_character_,
      fdr_directional < 0.05 ~ "FDR<0.05",
      fdr_directional < 0.10 ~ "FDR<0.10",
      TRUE ~ "FDR>=0.10"
    )
  ) %>%
  select(dataset, dataset_source, program, row_label, family, sig_label, gs_name, term_label,
         n_genes, median_age_effect, null_q025, null_q975, empirical_p_directional,
         fdr_directional, hsf_proteostasis, genes_tested)

sf4d_genes <- c("HSPA1A", "HSPA1B", "HSPA6", "DNAJB5", "HSF1", "STIP1", "DNAJB1", "AHSA1",
                "FKBP4", "BAG3", "HSPH1", "DNAJA1", "HSP90AB1", "HSPE1", "HSP90AA1",
                "HSPA8", "HSPD1")
sf4d <- read_csv(file.path(stage_dir, "SF4_participant_aware_human_blood_SR_aging_merged_gene_table.csv"), show_col_types = FALSE) %>%
  filter(dataset == "JenAge", gene %in% sf4d_genes) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(
    dataset_source = "JenAge blood RNA-seq; GSE103232 and GSE75337",
    age_class = if_else(age_logFC_per_decade < 0, "Age-down", "Not age-down")
  ) %>%
  arrange(desc(age_logFC_per_decade), gene) %>%
  mutate(plot_order = row_number()) %>%
  select(dataset, dataset_source, gene, age_logFC_per_decade, age_p, age_FDR, age_class,
         participant_aware_phase_logFC, participant_aware_phase_p, participant_aware_phase_FDR,
         delta_mesor, mesor_p, mesor_FDR, plot_order)

write_csv(sf4b, file.path(stage_dir, "SF4b_JenAge_Reactome_module_ranking.csv"))
write_csv(sf4c, file.path(stage_dir, "SF4c_JenAge_representative_Reactome_modules.csv"))
write_csv(sf4d, file.path(stage_dir, "SF4d_JenAge_proteostasis_gene_aging_direction.csv"))
write_csv(null_summary, file.path(stage_dir, "SF4b_JenAge_Reactome_module_null_summary.csv"))

if (!file.exists(backup_path)) {
  file.copy(workbook_path, backup_path, overwrite = FALSE)
}

make_matrix <- function(source_label, description, dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  width <- max(1, ncol(dat))
  top <- matrix(NA, nrow = 2, ncol = width)
  top[1, 1] <- source_label
  top[2, 1] <- description
  headers <- matrix(names(dat), nrow = 1)
  rbind(top, headers, as.matrix(dat))
}

wb <- loadWorkbook(workbook_path)
original_order <- names(wb)
title_style <- createStyle(textDecoration = "bold", fgFill = "#D9EAF7")
note_style <- createStyle(fgFill = "#F3F6FA")
header_style <- createStyle(textDecoration = "bold", fgFill = "#EDEDED")

replace_sheet <- function(sheet, mat) {
  if (sheet %in% names(wb)) removeWorksheet(wb, sheet)
  addWorksheet(wb, sheet, gridLines = FALSE)
  writeData(wb, sheet, mat, colNames = FALSE, rowNames = FALSE, keepNA = FALSE)
  addStyle(wb, sheet, title_style, rows = 1, cols = seq_len(ncol(mat)), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet, note_style, rows = 2, cols = seq_len(ncol(mat)), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet, header_style, rows = 3, cols = seq_len(ncol(mat)), gridExpand = TRUE, stack = TRUE)
  freezePane(wb, sheet, firstActiveRow = 4)
  setColWidths(wb, sheet, cols = seq_len(ncol(mat)), widths = "auto")
}

replace_sheet(
  "SF4b",
  make_matrix(
    "Source: SF4b_JenAge_Reactome_module_ranking.csv",
    "JenAge blood-aging Reactome module ranking plotted in Supplementary Fig. 4b; source accessions GSE103232 and GSE75337.",
    sf4b
  )
)
replace_sheet(
  "SF4c",
  make_matrix(
    "Source: SF4c_JenAge_representative_Reactome_modules.csv",
    "Representative JenAge blood-aging Reactome module effects plotted in Supplementary Fig. 4c; source accessions GSE103232 and GSE75337.",
    sf4c
  )
)
replace_sheet(
  "SF4d",
  make_matrix(
    "Source: SF4d_JenAge_proteostasis_gene_aging_direction.csv",
    "Selected HSF/HSP and proteostasis genes in JenAge blood aging plotted in Supplementary Fig. 4d; source accessions GSE103232 and GSE75337.",
    sf4d
  )
)

worksheetOrder(wb) <- match(original_order, names(wb))
saveWorkbook(wb, workbook_path, overwrite = TRUE)

audit <- tibble(
  panel = c("SF4b", "SF4c", "SF4d", "SF4e", "SF4f", "SF4g"),
  healthy_blood_aging_dataset = "JenAge",
  accessions = "GSE103232; GSE75337",
  note = c(
    "Recomputed Reactome module medians from JenAge age effects.",
    "Recomputed selected Reactome module rows from JenAge age effects.",
    "Replaced GSE134080 representative genes with JenAge gene-level age effects.",
    "Already used JenAge participant-aware SR-aging overlap.",
    "Already used JenAge participant-aware Reactome ORA.",
    "Already used JenAge participant-aware Reactome family ORA."
  )
)
write_csv(audit, file.path(stage_dir, "SF4_JenAge_only_dataset_consistency_audit.csv"))

message("Corrected Supplementary Fig. 4 source sheets to JenAge-only: ", workbook_path)
