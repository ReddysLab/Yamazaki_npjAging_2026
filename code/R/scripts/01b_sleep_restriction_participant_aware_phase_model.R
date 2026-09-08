# Participant-aware phase/cosinor analysis for the Moller-Levet sleep-restriction dataset.
# This reruns the SR-control analysis from the public GSE39445 matrix while accounting
# for repeated sampling within each participant.

if (!exists("paths")) {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
  project_root_guess <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
  Sys.setenv(NPJAGING_PROJECT_ROOT = Sys.getenv("NPJAGING_PROJECT_ROOT", unset = project_root_guess))
  source(file.path(dirname(script_dir), "config.R"))
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(janitor)
  library(msigdbr)
  library(openxlsx)
})

set.seed(20260811)

out_dir <- file.path(paths$tables, "moller_levet_participant_aware")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

path_gse39445_matrix <- file.path(paths$source_data, "GSE39445", "GSE39445_series_matrix.txt.gz")
path_phase_map <- file.path(paths$source_data, "GSE39445", "curated_metadata", "GSE39445_sample_phase_mapping.csv")
require_file(path_gse39445_matrix, "GSE39445 series matrix")
require_file(path_phase_map, "GSE39445 sample phase mapping")

standardize_gene <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, " ///.*$", "")
  x <- str_replace(x, " //.*$", "")
  x <- str_trim(x)
  x <- na_if(x, "")
  toupper(x)
}

extract_gene_symbols_from_featuredata <- function(fd) {
  fd2 <- fd %>%
    as.data.frame() %>%
    rownames_to_column("probe_id") %>%
    clean_names()
  gene_col_candidates <- c("gene_symbol", "gene_symbols", "gene_symbol_ch1", "symbol", "gene", "gene_assignment")
  gene_col <- gene_col_candidates[gene_col_candidates %in% names(fd2)][1]
  if (is.na(gene_col) || length(gene_col) == 0) {
    fd2 <- fd2 %>% mutate(gene_symbol_use = probe_id)
  } else {
    fd2 <- fd2 %>% mutate(gene_symbol_use = .data[[gene_col]])
  }
  fd2 %>%
    mutate(gene_symbol_use = standardize_gene(gene_symbol_use)) %>%
    select(probe_id, gene_symbol_use)
}

collapse_expression_to_gene <- function(expr_mat, gene_map_tbl) {
  expr_mat %>%
    as.data.frame() %>%
    rownames_to_column("probe_id") %>%
    left_join(gene_map_tbl, by = "probe_id") %>%
    filter(!is.na(gene_symbol_use), gene_symbol_use != "") %>%
    group_by(gene_symbol_use) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("gene_symbol_use") %>%
    as.matrix()
}

clean_term <- function(x) {
  x %>%
    str_replace("^REACTOME_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence() %>%
    str_replace_all(regex("hsf1", ignore_case = TRUE), "HSF1") %>%
    str_replace_all(regex("bmal1", ignore_case = TRUE), "BMAL1") %>%
    str_replace_all(regex("clock", ignore_case = TRUE), "CLOCK") %>%
    str_replace_all(regex("n?pas2", ignore_case = TRUE), "NPAS2") %>%
    str_replace_all(regex("cry", ignore_case = TRUE), "CRY") %>%
    str_replace_all(regex("per", ignore_case = TRUE), "PER") %>%
    str_replace_all(regex("ros", ignore_case = TRUE), "ROS")
}

gset <- GEOquery::getGEO(filename = path_gse39445_matrix)
expr_gene <- collapse_expression_to_gene(Biobase::exprs(gset), extract_gene_symbols_from_featuredata(Biobase::fData(gset)))
pheno <- Biobase::pData(gset) %>% as_tibble() %>% clean_names()

sample_tbl <- pheno %>%
  transmute(
    sample_id = as.character(geo_accession),
    subject = as.character(subject_ch1),
    condition_raw = as.character(sleepprotocol_ch1),
    phase_num = suppressWarnings(as.numeric(circadianphase_ch1)),
    time_sample_taken = as.character(timesampletaken_ch1)
  ) %>%
  mutate(
    condition = case_when(
      condition_raw == "Sleep Extension" ~ "Control",
      condition_raw == "Sleep Restriction" ~ "SR",
      TRUE ~ NA_character_
    )
  )

if (all(is.na(sample_tbl$phase_num)) || any(is.na(sample_tbl$condition))) {
  phase_map <- read_csv(path_phase_map, show_col_types = FALSE) %>% clean_names()
  sample_col <- names(phase_map)[str_detect(names(phase_map), "sample|geo_accession|gsm")][1]
  phase_col <- names(phase_map)[str_detect(names(phase_map), "phase|clock|zt|time")][1]
  condition_col <- names(phase_map)[str_detect(names(phase_map), "group|condition|status|sr|control")][1]
  phase_map2 <- phase_map %>%
    transmute(
      sample_id = as.character(.data[[sample_col]]),
      phase_num_fallback = suppressWarnings(as.numeric(.data[[phase_col]])),
      condition_fallback = case_when(
        str_detect(str_to_lower(as.character(.data[[condition_col]])), "control|sleep extension") ~ "Control",
        str_detect(str_to_lower(as.character(.data[[condition_col]])), "sr|sleep restriction|restricted") ~ "SR",
        TRUE ~ NA_character_
      )
    )
  sample_tbl <- sample_tbl %>%
    left_join(phase_map2, by = "sample_id") %>%
    mutate(
      phase_num = coalesce(phase_num, phase_num_fallback),
      condition = coalesce(condition, condition_fallback)
    ) %>%
    select(sample_id, subject, condition_raw, condition, phase_num, time_sample_taken)
}

sample_tbl <- sample_tbl %>%
  filter(!is.na(sample_id), !is.na(subject), !is.na(condition), !is.na(phase_num)) %>%
  mutate(
    condition = factor(condition, levels = c("Control", "SR")),
    condition_sr = as.numeric(condition == "SR"),
    phase_mod = phase_num %% 24,
    phase_factor = factor(phase_num, levels = sort(unique(phase_num))),
    sin_phase = sin(2 * pi * phase_mod / 24),
    cos_phase = cos(2 * pi * phase_mod / 24),
    sr_sin = condition_sr * sin_phase,
    sr_cos = condition_sr * cos_phase
  )

common_samples <- intersect(colnames(expr_gene), sample_tbl$sample_id)
sample_tbl <- sample_tbl %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))
expr_gene <- expr_gene[, sample_tbl$sample_id, drop = FALSE]
stopifnot(identical(colnames(expr_gene), sample_tbl$sample_id))

subject_block <- factor(sample_tbl$subject)

# Phase-factor model: SR-control mean effect while adjusting for discrete circadian phase.
design_phase <- model.matrix(~ condition_sr + phase_factor, data = sample_tbl)
cor_phase <- duplicateCorrelation(expr_gene, design_phase, block = subject_block)
fit_phase <- lmFit(expr_gene, design_phase, block = subject_block, correlation = cor_phase$consensus)
fit_phase <- eBayes(fit_phase)
phase_tbl <- topTable(fit_phase, coef = "condition_sr", number = Inf, sort.by = "none") %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    participant_aware_phase_logFC = logFC,
    participant_aware_phase_t = t,
    participant_aware_phase_p = P.Value,
    participant_aware_phase_FDR = adj.P.Val
  )

# Cosinor-style model: condition mean plus SR-specific phase-waveform terms.
design_cos <- model.matrix(~ condition_sr + sin_phase + cos_phase + sr_sin + sr_cos, data = sample_tbl)
cor_cos <- duplicateCorrelation(expr_gene, design_cos, block = subject_block)
fit_cos <- lmFit(expr_gene, design_cos, block = subject_block, correlation = cor_cos$consensus)
fit_cos <- eBayes(fit_cos)

mesor_tbl <- topTable(fit_cos, coef = "condition_sr", number = Inf, sort.by = "none") %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    delta_mesor = logFC,
    mesor_t = t,
    mesor_p = P.Value,
    mesor_FDR = adj.P.Val
  )

wave_tbl <- topTable(fit_cos, coef = c("sr_sin", "sr_cos"), number = Inf, sort.by = "none") %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    waveform_F = F,
    waveform_p = P.Value,
    waveform_FDR = adj.P.Val
  )

coef_tbl <- as.data.frame(fit_cos$coefficients) %>%
  rownames_to_column("gene") %>%
  transmute(
    gene,
    sr_sin_coef = sr_sin,
    sr_cos_coef = sr_cos,
    rhythmic_SR_effect_amplitude = sqrt(sr_sin_coef^2 + sr_cos_coef^2),
    SR_effect_acrophase_h = (atan2(sr_sin_coef, sr_cos_coef) * 24 / (2 * pi)) %% 24
  )

gene_results <- phase_tbl %>%
  full_join(mesor_tbl, by = "gene") %>%
  full_join(wave_tbl, by = "gene") %>%
  full_join(coef_tbl, by = "gene") %>%
  arrange(mesor_p)

write_csv(sample_tbl, file.path(out_dir, "MollerLevet_participant_aware_sample_table.csv"))
write_csv(gene_results, file.path(out_dir, "MollerLevet_participant_aware_phase_cosinor_gene_results.csv"))

selected_terms <- c(
  "AGE-receptor signaling" = "REACTOME_ADVANCED_GLYCOSYLATION_ENDPRODUCT_RECEPTOR_SIGNALING",
  "ROS detoxification" = "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
  "Interleukin-1 signaling" = "REACTOME_INTERLEUKIN_1_SIGNALING",
  "Interleukin-1 family signaling" = "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",
  "BMAL1/CLOCK/NPAS2" = "REACTOME_BMAL1_CLOCK_NPAS2_ACTIVATES_CIRCADIAN_EXPRESSION",
  "CRY/PER degradation" = "REACTOME_DEGRADATION_OF_CRY_AND_PER_PROTEINS",
  "Circadian clock" = "REACTOME_CIRCADIAN_CLOCK",
  "HSF1 heat-shock regulation" = "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",
  "Heat-stress response" = "REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS",
  "Insulin processing" = "REACTOME_INSULIN_PROCESSING"
)

reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  transmute(module_id = gs_name, gene = standardize_gene(gene_symbol)) %>%
  filter(module_id %in% unname(selected_terms), gene %in% gene_results$gene) %>%
  distinct()

module_gene_tbl <- tibble(module_label = names(selected_terms), module_id = unname(selected_terms)) %>%
  left_join(reactome, by = "module_id") %>%
  filter(!is.na(gene))

circular_mean_h <- function(hours) {
  hours <- hours[is.finite(hours)]
  if (!length(hours)) return(NA_real_)
  (atan2(mean(sin(2 * pi * hours / 24)), mean(cos(2 * pi * hours / 24))) * 24 / (2 * pi)) %% 24
}

module_summary <- module_gene_tbl %>%
  left_join(gene_results, by = "gene") %>%
  group_by(module_label, module_id) %>%
  summarise(
    n_genes = n_distinct(gene),
    median_delta_mesor = median(delta_mesor, na.rm = TRUE),
    median_rhythmic_SR_effect_amplitude = median(rhythmic_SR_effect_amplitude, na.rm = TRUE),
    circular_mean_SR_effect_acrophase_h = circular_mean_h(SR_effect_acrophase_h),
    median_participant_aware_phase_logFC = median(participant_aware_phase_logFC, na.rm = TRUE),
    genes_tested = paste(sort(unique(gene)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 3)

universe <- gene_results %>%
  filter(is.finite(delta_mesor), is.finite(rhythmic_SR_effect_amplitude)) %>%
  pull(gene)

n_perm <- 5000
null_tbl <- lapply(seq_len(n_perm), function(iter) {
  bind_rows(lapply(seq_len(nrow(module_summary)), function(i) {
    sampled <- sample(universe, module_summary$n_genes[i], replace = FALSE)
    d <- gene_results[match(sampled, gene_results$gene), ]
    tibble(
      iteration = iter,
      module_label = module_summary$module_label[i],
      null_median_delta_mesor = median(d$delta_mesor, na.rm = TRUE),
      null_median_amplitude = median(d$rhythmic_SR_effect_amplitude, na.rm = TRUE)
    )
  }))
}) %>% bind_rows()

module_summary <- module_summary %>%
  left_join(
    null_tbl %>%
      group_by(module_label) %>%
      summarise(
        null_mesor_q025 = quantile(null_median_delta_mesor, 0.025, na.rm = TRUE),
        null_mesor_q975 = quantile(null_median_delta_mesor, 0.975, na.rm = TRUE),
        null_amp_q025 = quantile(null_median_amplitude, 0.025, na.rm = TRUE),
        null_amp_q975 = quantile(null_median_amplitude, 0.975, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "module_label"
  ) %>%
  rowwise() %>%
  mutate(
    empirical_p_more_negative = max(mean(null_tbl$null_median_delta_mesor[null_tbl$module_label == module_label] <= median_delta_mesor), 1 / n_perm),
    empirical_p_more_positive = max(mean(null_tbl$null_median_delta_mesor[null_tbl$module_label == module_label] >= median_delta_mesor), 1 / n_perm),
    empirical_p_amp_more_extreme = max(mean(null_tbl$null_median_amplitude[null_tbl$module_label == module_label] >= median_rhythmic_SR_effect_amplitude), 1 / n_perm),
    direction = if_else(median_delta_mesor < 0, "negative", "positive"),
    empirical_p_directional = if_else(direction == "negative", empirical_p_more_negative, empirical_p_more_positive)
  ) %>%
  ungroup() %>%
  mutate(
    empirical_FDR_directional = p.adjust(empirical_p_directional, method = "BH"),
    empirical_FDR_amplitude = p.adjust(empirical_p_amp_more_extreme, method = "BH")
  ) %>%
  arrange(median_delta_mesor)

write_csv(module_gene_tbl, file.path(out_dir, "MollerLevet_participant_aware_reactome_module_gene_membership.csv"))
write_csv(module_summary, file.path(out_dir, "MollerLevet_participant_aware_reactome_module_summary.csv"))
write_csv(null_tbl, file.path(out_dir, "MollerLevet_participant_aware_reactome_module_random_null.csv"))

hsf_hsp_genes <- c(
  "HSF1", "HSPA5", "HSPD1", "HSPH1", "HSPA8", "HSP90AB1", "HSP90AA1",
  "HSPA1A", "DNAJB1", "BAG3", "AHSA1", "STIP1", "HSPE1", "HSPB1", "CRYAB",
  "NFKB2", "PRDX2"
)
hsf_hsp_audit <- gene_results %>%
  filter(gene %in% hsf_hsp_genes) %>%
  arrange(match(gene, hsf_hsp_genes))
write_csv(hsf_hsp_audit, file.path(out_dir, "MollerLevet_participant_aware_HSF_HSP_gene_audit.csv"))

model_summary <- tibble(
  analysis = c("phase_factor", "cosinor"),
  model = c(
    "expression ~ condition_sr + phase_factor; limma duplicateCorrelation block=subject",
    "expression ~ condition_sr + sin(2*pi*phase/24) + cos(2*pi*phase/24) + condition_sr:sin + condition_sr:cos; limma duplicateCorrelation block=subject"
  ),
  n_genes = nrow(expr_gene),
  n_samples = ncol(expr_gene),
  n_subjects = n_distinct(sample_tbl$subject),
  duplicate_correlation = c(cor_phase$consensus, cor_cos$consensus)
)
write_csv(model_summary, file.path(out_dir, "MollerLevet_participant_aware_model_summary.csv"))

write.xlsx(
  list(
    model_summary = model_summary,
    gene_results = gene_results,
    HSF_HSP_gene_audit = hsf_hsp_audit,
    reactome_module_summary = module_summary,
    reactome_module_gene_membership = module_gene_tbl
  ),
  file.path(out_dir, "MollerLevet_participant_aware_phase_cosinor_results.xlsx"),
  overwrite = TRUE
)

message("Participant-aware Moller-Levet phase/cosinor rerun complete.")
print(model_summary)
print(module_summary %>% select(module_label, n_genes, median_delta_mesor, empirical_p_directional, empirical_FDR_directional, median_rhythmic_SR_effect_amplitude, empirical_p_amp_more_extreme))
print(hsf_hsp_audit %>% select(gene, delta_mesor, mesor_p, mesor_FDR, participant_aware_phase_logFC, participant_aware_phase_p, participant_aware_phase_FDR))
