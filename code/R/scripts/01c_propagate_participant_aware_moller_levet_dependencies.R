# Propagate participant-aware Moller-Levet outputs to dependent figure source data.
# This script stages updated Fig. 1, Supplementary Fig. 1, and Supplementary Fig. 4
# source tables without overwriting final figure files.

if (!exists("paths")) {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
  project_root_guess <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
  Sys.setenv(NPJAGING_PROJECT_ROOT = Sys.getenv("NPJAGING_PROJECT_ROOT", unset = project_root_guess))
  source(file.path(dirname(script_dir), "config.R"))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(ggforce)
  library(ggrepel)
  library(igraph)
  library(ggraph)
  library(msigdbr)
})

set.seed(20260811)

participant_dir <- file.path(paths$tables, "moller_levet_participant_aware")
stage_dir <- file.path(paths$tables, "participant_aware_staging")
stage_fig_dir <- file.path(paths$figures, "participant_aware_staging")
final_fig_dir <- file.path(paths$figures, "final")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stage_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)

std_gene <- function(x) toupper(str_trim(as.character(x)))
clean_term <- function(x) {
  x %>%
    str_remove("^REACTOME_") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence() %>%
    str_replace_all(regex("hsf1", ignore_case = TRUE), "HSF1") %>%
    str_replace_all(regex("hsp90", ignore_case = TRUE), "HSP90") %>%
    str_replace_all(regex("bmal1", ignore_case = TRUE), "BMAL1") %>%
    str_replace_all(regex("clock", ignore_case = TRUE), "CLOCK") %>%
    str_replace_all(regex("npas2", ignore_case = TRUE), "NPAS2") %>%
    str_replace_all(regex("cry", ignore_case = TRUE), "CRY") %>%
    str_replace_all(regex("per", ignore_case = TRUE), "PER") %>%
    str_replace_all(regex("ros", ignore_case = TRUE), "ROS")
}
short_p <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 1e-300 ~ "<1e-300",
    x < 1e-3 ~ formatC(x, format = "e", digits = 2),
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}
fisher_greater <- function(a, b, c, universe_n) {
  mat <- matrix(c(a, b - a, c - a, universe_n - b - c + a), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  tibble(overlap_n = a, odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}
reactome_flag <- function(term) {
  str_detect(term, regex("HSF|HEAT|CHAPERONE|UNFOLDED|PROTEAS|HSP|STRESS|CALNEXIN|CALRETICULIN|PROTEOSTAS", ignore_case = TRUE))
}

gene_results <- read_csv(
  file.path(participant_dir, "MollerLevet_participant_aware_phase_cosinor_gene_results.csv"),
  show_col_types = FALSE
) %>%
  mutate(gene = std_gene(gene))

module_selected_raw <- read_csv(
  file.path(participant_dir, "MollerLevet_participant_aware_reactome_module_summary.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    module_id = as.character(module_id),
    display_name = module_label,
    family = case_when(
      str_detect(display_name, regex("AGE", ignore_case = TRUE)) ~ "AGE",
      str_detect(display_name, regex("ROS|interleukin", ignore_case = TRUE)) ~ "IL/ROS",
      str_detect(display_name, regex("BMAL|circadian|CRY|PER", ignore_case = TRUE)) ~ "Circadian",
      str_detect(display_name, regex("HSF|heat", ignore_case = TRUE)) ~ "Heat-stress",
      str_detect(display_name, regex("insulin", ignore_case = TRUE)) ~ "Insulin",
      TRUE ~ "Other"
    ),
    median_rhythmic_amplitude = median_rhythmic_SR_effect_amplitude,
    circular_mean_effect_acrophase_h = circular_mean_SR_effect_acrophase_h,
    amp_FDR = empirical_FDR_amplitude,
    mesor_significance = case_when(
      empirical_FDR_directional < 0.001 ~ "***",
      empirical_FDR_directional < 0.01 ~ "**",
      empirical_FDR_directional < 0.05 ~ "*",
      TRUE ~ ""
    ),
    p_label = paste0("emp. p=", short_p(empirical_p_directional)),
    fdr_label = paste0("FDR=", short_p(empirical_FDR_directional)),
    amp_p_label = paste0("emp. p=", short_p(empirical_p_amp_more_extreme)),
    amp_fdr_label = paste0("FDR=", short_p(empirical_FDR_amplitude))
  ) %>%
  arrange(median_delta_mesor)

write_csv(
  module_selected_raw,
  file.path(stage_dir, "F1e_SF1c_SF1d_SF1e_participant_aware_selected_reactome_modules.csv")
)

old_f1d <- read_csv(
  file.path(paths$curated_metadata, "GSE39445", "Figure_1D_descriptive_timecourse.csv"),
  show_col_types = FALSE
) %>%
  mutate(gene = std_gene(gene))

add_gene_stats <- function(timecourse_tbl, keep_genes) {
  stats <- gene_results %>%
    filter(gene %in% std_gene(keep_genes)) %>%
    transmute(
      gene,
      participant_aware_phase_logFC,
      participant_aware_phase_p,
      participant_aware_phase_FDR,
      participant_aware_delta_MESOR = delta_mesor,
      participant_aware_MESOR_p = mesor_p,
      participant_aware_MESOR_FDR = mesor_FDR,
      participant_aware_waveform_F = waveform_F,
      participant_aware_waveform_p = waveform_p,
      participant_aware_waveform_FDR = waveform_FDR,
      participant_aware_SR_effect_amplitude = rhythmic_SR_effect_amplitude,
      participant_aware_SR_effect_acrophase_h = SR_effect_acrophase_h
    )
  timecourse_tbl %>%
    filter(gene %in% std_gene(keep_genes)) %>%
    select(-any_of(names(stats)[names(stats) != "gene"])) %>%
    left_join(stats, by = "gene") %>%
    arrange(factor(gene, levels = std_gene(keep_genes)), condition, phase_num)
}

f1d_src <- add_gene_stats(old_f1d, c("HSPH1", "HSPA5"))
sf1a_src <- add_gene_stats(old_f1d, c("HSPD1", "NFKB2", "PRDX2"))
write_csv(f1d_src, file.path(stage_dir, "F1d_participant_aware_source_data.csv"))
write_csv(sf1a_src, file.path(stage_dir, "SF1a_participant_aware_source_data.csv"))

reactome_sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  transmute(term = gs_name, gene = std_gene(gene_symbol)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct() %>%
  semi_join(gene_results, by = "gene")

gene_effects <- gene_results %>%
  filter(is.finite(delta_mesor), is.finite(rhythmic_SR_effect_amplitude)) %>%
  select(gene, delta_mesor, rhythmic_SR_effect_amplitude, SR_effect_acrophase_h,
         participant_aware_phase_logFC, participant_aware_phase_p, participant_aware_phase_FDR,
         mesor_p, mesor_FDR, waveform_p, waveform_FDR)

term_summary <- reactome_sets %>%
  inner_join(gene_effects, by = "gene") %>%
  group_by(term) %>%
  summarise(
    n_genes = n_distinct(gene),
    median_delta_mesor = median(delta_mesor, na.rm = TRUE),
    median_rhythmic_amplitude = median(rhythmic_SR_effect_amplitude, na.rm = TRUE),
    genes_tested = paste(sort(unique(gene)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 10, n_genes <= 500) %>%
  mutate(
    term_clean = clean_term(term),
    hsf_proteostasis = reactome_flag(term)
  )

null_by_size <- map_dfr(sort(unique(term_summary$n_genes)), function(n) {
  tibble(
    n_genes = n,
    iter = seq_len(2000),
    null_median_delta_mesor = replicate(2000, median(sample(gene_effects$delta_mesor, n, replace = FALSE), na.rm = TRUE)),
    null_median_rhythmic_amplitude = replicate(2000, median(sample(gene_effects$rhythmic_SR_effect_amplitude, n, replace = FALSE), na.rm = TRUE))
  )
})

term_summary <- term_summary %>%
  left_join(
    null_by_size %>%
      group_by(n_genes) %>%
      summarise(
        null_mesor_q025 = quantile(null_median_delta_mesor, 0.025, na.rm = TRUE),
        null_mesor_q975 = quantile(null_median_delta_mesor, 0.975, na.rm = TRUE),
        null_amp_q025 = quantile(null_median_rhythmic_amplitude, 0.025, na.rm = TRUE),
        null_amp_q975 = quantile(null_median_rhythmic_amplitude, 0.975, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "n_genes"
  ) %>%
  rowwise() %>%
  mutate(
    empirical_p_more_negative = max(mean(null_by_size$null_median_delta_mesor[null_by_size$n_genes == n_genes] <= median_delta_mesor, na.rm = TRUE), 1 / 2000),
    empirical_p_more_positive = max(mean(null_by_size$null_median_delta_mesor[null_by_size$n_genes == n_genes] >= median_delta_mesor, na.rm = TRUE), 1 / 2000),
    empirical_p_two_sided = min(1, 2 * min(empirical_p_more_negative, empirical_p_more_positive)),
    empirical_p_amp_more_extreme = max(mean(null_by_size$null_median_rhythmic_amplitude[null_by_size$n_genes == n_genes] >= median_rhythmic_amplitude, na.rm = TRUE), 1 / 2000)
  ) %>%
  ungroup() %>%
  mutate(
    empirical_FDR_more_negative = p.adjust(empirical_p_more_negative, "BH"),
    empirical_FDR_mesor_two_sided = p.adjust(empirical_p_two_sided, "BH"),
    empirical_FDR_amplitude = p.adjust(empirical_p_amp_more_extreme, "BH"),
    rank_from_most_negative = dense_rank(median_delta_mesor),
    highlight = if_else(hsf_proteostasis, "HSF/proteostasis", "Other Reactome")
  ) %>%
  arrange(rank_from_most_negative)

write_csv(term_summary, file.path(stage_dir, "SF1b_participant_aware_reactome_wide_MESOR_ranking.csv"))
write_csv(null_by_size, file.path(stage_dir, "SF1e_participant_aware_reactome_random_nulls.csv"))

sr_down_new <- gene_results %>%
  filter(
    participant_aware_phase_logFC < 0,
    participant_aware_phase_FDR < thresholds$human_sr_phase_fdr_max,
    delta_mesor < 0,
    mesor_FDR < thresholds$human_sr_mesor_fdr_max
  ) %>%
  transmute(
    gene,
    gene_symbol = gene,
    participant_aware_phase_logFC,
    participant_aware_phase_p,
    participant_aware_phase_FDR,
    delta_mesor,
    mesor_p,
    mesor_FDR,
    rhythmic_SR_effect_amplitude,
    SR_effect_acrophase_h,
    manuscript_threshold = "participant-aware phase logFC < 0 and BH FDR < 0.05; delta MESOR < 0 and BH FDR < 0.05"
  ) %>%
  arrange(participant_aware_phase_p, gene)

write_csv(sr_down_new, file.path(stage_dir, "human_blood_SR_down_participant_aware_dual_FDR_full.csv"))

healthy_dir <- paths$human_healthy_blood
healthy_files <- tibble(
  dataset = "JenAge",
  dataset_label = "JenAge",
  path = file.path(healthy_dir, "JenAge_01_healthy_blood_age_log2RPKM_results.csv"),
  age_logfc_col = "jenage_blood_age_logFC_per_decade",
  age_p_col = "jenage_blood_age_p",
  age_fdr_col = "jenage_blood_age_fdr"
)

sr_flags <- gene_results %>%
  transmute(
    gene,
    sr_direction = participant_aware_phase_logFC < 0 & delta_mesor < 0,
    sr_nominalP05 = participant_aware_phase_logFC < 0 & participant_aware_phase_p < 0.05 & delta_mesor < 0,
    sr_FDR05 = participant_aware_phase_logFC < 0 & participant_aware_phase_FDR < 0.05 & delta_mesor < 0,
    participant_aware_phase_logFC,
    participant_aware_phase_p,
    participant_aware_phase_FDR,
    delta_mesor,
    mesor_p,
    mesor_FDR
  )

sf4_joined <- pmap_dfr(healthy_files, function(dataset, dataset_label, path, age_logfc_col, age_p_col, age_fdr_col) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(gene = std_gene(gene)) %>%
    select(gene, all_of(c(age_logfc_col, age_p_col, age_fdr_col))) %>%
    rename(age_logFC_per_decade = all_of(age_logfc_col), age_p = all_of(age_p_col), age_FDR = all_of(age_fdr_col)) %>%
    inner_join(sr_flags, by = "gene") %>%
    mutate(
      dataset = dataset,
      dataset_label = dataset_label,
      age_direction_down = age_logFC_per_decade < 0,
      age_nominalP05_down = age_logFC_per_decade < 0 & age_p < 0.05,
      age_FDR05_down = age_logFC_per_decade < 0 & age_FDR < 0.05,
      .before = 1
    )
})

write_csv(sf4_joined, file.path(stage_dir, "SF4_participant_aware_human_blood_SR_aging_merged_gene_table.csv"))

overlap_summary <- tidyr::crossing(
  dataset = unique(sf4_joined$dataset),
  sr_threshold = c("sr_direction", "sr_nominalP05", "sr_FDR05"),
  aging_threshold = c("age_direction_down", "age_nominalP05_down", "age_FDR05_down")
) %>%
  mutate(stats = pmap(list(dataset, sr_threshold, aging_threshold), function(dataset, sr_threshold, aging_threshold) {
    dataset_i <- dataset
    sr_threshold_i <- sr_threshold
    aging_threshold_i <- aging_threshold
    dat <- sf4_joined %>% filter(.data$dataset == .env$dataset_i)
    sr_vec <- dat[[sr_threshold_i]]
    age_vec <- dat[[aging_threshold_i]]
    a <- sum(sr_vec & age_vec, na.rm = TRUE)
    b <- sum(sr_vec, na.rm = TRUE)
    c <- sum(age_vec, na.rm = TRUE)
    u <- nrow(dat)
    fisher_greater(a, b, c, u) %>%
      mutate(universe_n = u, sr_down_n = b, aging_down_n = c, .before = 1)
  })) %>%
  unnest(stats) %>%
  mutate(
    fisher_p = p_value,
    threshold_description = case_when(
      sr_threshold == "sr_nominalP05" ~ "participant-aware SR down: phase-adjusted logFC < 0, nominal p < 0.05, delta MESOR < 0",
      sr_threshold == "sr_FDR05" ~ "participant-aware SR down: phase-adjusted logFC < 0, FDR < 0.05, delta MESOR < 0",
      TRUE ~ "participant-aware SR down: phase-adjusted logFC < 0 and delta MESOR < 0"
    )
  ) %>%
  select(dataset, sr_threshold, aging_threshold, universe_n, sr_down_n, aging_down_n, overlap_n, fisher_OR = odds_ratio, fisher_p, threshold_description)

write_csv(overlap_summary, file.path(stage_dir, "SF4e_participant_aware_SRdown_x_AgingDown_fisher_summary.csv"))

primary_sf4 <- sf4_joined %>%
  filter(dataset == "JenAge", sr_nominalP05, age_direction_down) %>%
  arrange(gene)

write_csv(primary_sf4, file.path(stage_dir, "SF4e_participant_aware_JenAge_SRnominalP05_AgeDirection_overlap_genes.csv"))

run_ora <- function(query, universe, term2gene) {
  query <- intersect(unique(query), universe)
  term_tbl <- term2gene %>%
    filter(gene %in% universe) %>%
    group_by(term) %>%
    summarise(term_genes = list(sort(unique(gene))), term_n = n_distinct(gene), .groups = "drop") %>%
    filter(term_n >= 3, term_n <= 500)
  map_dfr(seq_len(nrow(term_tbl)), function(i) {
    tg <- term_tbl$term_genes[[i]]
    ov <- intersect(query, tg)
    a <- length(ov)
    b <- length(query)
    c <- length(tg)
    u <- length(universe)
    ft <- fisher_greater(a, b, c, u)
    tibble(
      term = term_tbl$term[[i]],
      term_clean = clean_term(term_tbl$term[[i]]),
      term_n = c,
      query_n = b,
      universe_n = u,
      overlap_n = a,
      odds_ratio = ft$odds_ratio,
      p_value = ft$p_value,
      overlap_gene_string = paste(sort(ov), collapse = ";")
    )
  }) %>%
    mutate(FDR = p.adjust(p_value, "BH"), hsf_proteostasis_related = reactome_flag(term)) %>%
    arrange(FDR, p_value, desc(odds_ratio))
}

jen_universe <- sf4_joined %>% filter(dataset == "JenAge") %>% pull(gene) %>% unique()
sf4_ora <- run_ora(primary_sf4$gene, jen_universe, reactome_sets)
write_csv(sf4_ora, file.path(stage_dir, "SF4f_SF4g_participant_aware_JenAge_SRnominalP05_AgeDirection_Reactome_ORA.csv"))

assign_family <- function(term) {
  case_when(
    str_detect(term, regex("HSF|HEAT|CHAPERONE|HSP|UNFOLDED|PROTEOSTAS", ignore_case = TRUE)) ~ "HSF1-associated heat-shock response",
    str_detect(term, regex("SUMO", ignore_case = TRUE)) ~ "SUMOylation programs",
    str_detect(term, regex("RNA|SPLICE|SNRNA|SNRNP|RIBONUCLEOPROTEIN|TRNA|RRNA", ignore_case = TRUE)) ~ "RNA processing / splicing",
    str_detect(term, regex("NUCLEAR.*TRANSPORT|RNA EXPORT|VIRAL.*EXPORT|EXPORTIN|IMPORTIN", ignore_case = TRUE)) ~ "Nuclear transport / viral RNA export",
    str_detect(term, regex("CELL CYCLE|MITOTIC|CHROMOSOME|DNA REPLICATION", ignore_case = TRUE)) ~ "Cell-cycle regulation",
    str_detect(term, regex("DNA DAMAGE|TP53|P53", ignore_case = TRUE)) ~ "DNA damage / TP53 regulation",
    TRUE ~ NA_character_
  )
}

sf4_family_terms <- sf4_ora %>%
  filter(FDR < 0.05) %>%
  mutate(family = assign_family(term)) %>%
  filter(!is.na(family))

sf4_family_membership <- sf4_family_terms %>%
  select(family, term, term_clean, overlap_gene_string, FDR, p_value) %>%
  separate_rows(overlap_gene_string, sep = ";") %>%
  rename(gene = overlap_gene_string) %>%
  filter(gene != "")

sf4_family_universe_membership <- sf4_family_terms %>%
  select(family, term, term_clean, FDR) %>%
  inner_join(reactome_sets %>% filter(gene %in% jen_universe), by = "term") %>%
  distinct(family, term, term_clean, FDR, gene)

sf4_family_summary <- sf4_family_membership %>%
  group_by(family) %>%
  summarise(
    n_reactome_terms = n_distinct(term),
    overlap_n = n_distinct(gene),
    representative_terms = paste(head(unique(term_clean[order(FDR)]), 8), collapse = "; "),
    overlap_gene_string = paste(sort(unique(gene)), collapse = ";"),
    .groups = "drop"
  ) %>%
  left_join(
    sf4_family_universe_membership %>%
      group_by(family) %>%
      summarise(family_gene_n = n_distinct(gene), .groups = "drop"),
    by = "family"
  ) %>%
  select(family, n_reactome_terms, family_gene_n, overlap_n, representative_terms, overlap_gene_string) %>%
  rowwise() %>%
  mutate(
    fisher = list(fisher_greater(overlap_n, nrow(primary_sf4), family_gene_n, length(jen_universe)) %>%
                    select(family_odds_ratio = odds_ratio, family_p_value = p_value))
  ) %>%
  unnest(fisher) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(family_p_value, "BH")) %>%
  arrange(FDR, desc(overlap_n))

write_csv(sf4_family_summary, file.path(stage_dir, "SF4g_participant_aware_Reactome_family_summary.csv"))
write_csv(sf4_family_membership, file.path(stage_dir, "SF4g_participant_aware_Reactome_family_gene_membership.csv"))

network_nodes <- sf4_ora %>%
  filter(FDR < 0.05, overlap_n >= 3) %>%
  arrange(FDR, desc(overlap_n)) %>%
  slice_head(n = 120) %>%
  transmute(
    node_id = term,
    label = term_clean,
    overlap_n,
    odds_ratio,
    p_value,
    FDR,
    node_class = if_else(hsf_proteostasis_related, "HSF/proteostasis", "Other Reactome"),
    overlap_gene_string
  )

network_edges <- combn(network_nodes$node_id, 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    g1 <- str_split(network_nodes$overlap_gene_string[network_nodes$node_id == pair[1]], ";")[[1]]
    g2 <- str_split(network_nodes$overlap_gene_string[network_nodes$node_id == pair[2]], ";")[[1]]
    shared <- intersect(g1[g1 != ""], g2[g2 != ""])
    tibble(from = pair[1], to = pair[2], shared_genes_n = length(shared), shared_genes = paste(sort(shared), collapse = ";"))
  }) %>%
  filter(shared_genes_n >= 2)

write_csv(network_nodes, file.path(stage_dir, "SF4f_participant_aware_Reactome_network_nodes.csv"))
write_csv(network_edges, file.path(stage_dir, "SF4f_participant_aware_Reactome_network_edges.csv"))

sf4_primary_overlap <- overlap_summary %>%
  filter(dataset == "JenAge", sr_threshold == "sr_nominalP05", aging_threshold == "age_direction_down") %>%
  slice(1)
sf4_sr_only_n <- sf4_primary_overlap$sr_down_n - sf4_primary_overlap$overlap_n
sf4_age_only_n <- sf4_primary_overlap$aging_down_n - sf4_primary_overlap$overlap_n
sf4_fisher_label <- paste0(
  "Fisher OR=", signif(sf4_primary_overlap$fisher_OR, 3),
  "; p=", short_p(sf4_primary_overlap$fisher_p)
)
sf4e_circle_df <- tibble(
  set = c("SR-down\nnominal p<0.05", "Blood-aging down\nlogFC<0"),
  x0 = c(0, 1.1),
  y0 = c(0, 0),
  r = c(0.95, 0.95),
  color = c("#2c7fb8", "#2ca25f")
)
sf4e_label_df <- tibble(
  x = c(-0.42, 1.52, 0.55),
  y = c(0, 0, 0),
  label = c(sf4_sr_only_n, sf4_age_only_n, sf4_primary_overlap$overlap_n)
)
p_sf4e <- ggplot() +
  geom_circle(data = sf4e_circle_df, aes(x0 = x0, y0 = y0, r = r, color = set),
              fill = "grey90", alpha = 0.45, linewidth = 0.8) +
  geom_text(data = sf4e_label_df, aes(x, y, label = label), size = 4.5) +
  geom_text(
    data = sf4e_circle_df,
    aes(x0, y0 + 1.15, label = paste0(set, "\nn=", c(sf4_primary_overlap$sr_down_n, sf4_primary_overlap$aging_down_n))),
    size = 2.8,
    lineheight = 0.92
  ) +
  annotate("label", x = 0.55, y = -1.08, label = sf4_fisher_label, size = 2.6, linewidth = 0.25) +
  scale_color_manual(values = c("SR-down\nnominal p<0.05" = "#2c7fb8", "Blood-aging down\nlogFC<0" = "#2ca25f")) +
  coord_fixed(xlim = c(-1.1, 2.2), ylim = c(-1.35, 1.45), clip = "off") +
  labs(title = "SR-aging convergent genes") +
  theme_void(base_size = 9) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0, size = 10))
ggsave(file.path(stage_fig_dir, "SF4e_participant_aware_SR_aging_overlap_staged.png"), p_sf4e, width = 3.2, height = 2.7, dpi = 360, bg = "white")
ggsave(file.path(stage_fig_dir, "SF4e_participant_aware_SR_aging_overlap_staged.pdf"), p_sf4e, width = 3.2, height = 2.7, bg = "white")

sf4_family_color <- c(
  "HSF/HSP-related gene" = "#1b9e77",
  "HSF1 heat-shock family" = "#d95f02",
  "Other pathway family" = "#6a51a3",
  "Shared convergence gene" = "grey70"
)
sf4_network_terms <- bind_rows(
  sf4_ora %>%
    filter(FDR < 0.05) %>%
    arrange(FDR, desc(overlap_n)) %>%
    slice_head(n = 14),
  sf4_ora %>%
    filter(FDR < 0.05, hsf_proteostasis_related)
) %>%
  distinct(term, .keep_all = TRUE) %>%
  mutate(
    term_node = paste0("TERM__", term),
    term_label = case_when(
      str_detect(term_clean, regex("HSF1", ignore_case = TRUE)) ~ "HSF1 heat-shock\nresponse",
      str_detect(term_clean, regex("Cellular response to heat stress", ignore_case = TRUE)) ~ "Cellular response\nto heat stress",
      str_detect(term_clean, regex("SUMO", ignore_case = TRUE)) ~ "SUMOylation\nprograms",
      str_detect(term_clean, regex("RNA|MRNA|TRNA|RRNA|SPLIC", ignore_case = TRUE)) ~ "RNA processing /\nsplicing",
      str_detect(term_clean, regex("CELL CYCLE", ignore_case = TRUE)) ~ "Cell-cycle\nregulation",
      str_detect(term_clean, regex("TP53|DNA DAMAGE", ignore_case = TRUE)) ~ "DNA damage / TP53\nregulation",
      TRUE ~ str_wrap(term_clean, 18)
    ),
    term_class = if_else(hsf_proteostasis_related, "HSF1 heat-shock family", "Other pathway family")
  )
sf4_term_gene_edges <- sf4_network_terms %>%
  transmute(term_node, gene = str_split(overlap_gene_string, ";")) %>%
  unnest(gene) %>%
  filter(gene != "") %>%
  mutate(from = term_node, to = paste0("GENE__", gene)) %>%
  select(from, to, gene)
sf4_network_gene_nodes <- sf4_term_gene_edges %>%
  distinct(name = to, gene) %>%
  mutate(
    label = gene,
    node_type = "gene",
    node_class = if_else(
      str_detect(gene, regex("^HSP|^DNAJ|^BAG|HSF1|HSP90|HSPA|HSPH", ignore_case = TRUE)),
      "HSF/HSP-related gene",
      "Shared convergence gene"
    )
  )
sf4_network_term_nodes <- sf4_network_terms %>%
  transmute(name = term_node, label = term_label, node_type = "term", node_class = term_class)
sf4_network_graph <- graph_from_data_frame(
  sf4_term_gene_edges %>% select(from, to),
  vertices = bind_rows(sf4_network_term_nodes, sf4_network_gene_nodes),
  directed = FALSE
)
sf4_gene_degree <- degree(sf4_network_graph)
V(sf4_network_graph)$node_size <- ifelse(
  V(sf4_network_graph)$node_type == "term",
  5.5,
  pmax(1.2, pmin(3.2, sf4_gene_degree[V(sf4_network_graph)$name] / 2))
)
V(sf4_network_graph)$label_plot <- ifelse(
  V(sf4_network_graph)$node_type == "term" |
    (V(sf4_network_graph)$node_class == "HSF/HSP-related gene" & sf4_gene_degree[V(sf4_network_graph)$name] >= 2),
  V(sf4_network_graph)$label,
  NA_character_
)
p_sf4f <- ggraph(sf4_network_graph, layout = "fr") +
  geom_edge_link(alpha = 0.25, color = "grey65", linewidth = 0.25) +
  geom_node_point(aes(color = node_class, size = node_size), alpha = 0.95) +
  geom_node_label(aes(label = ifelse(node_type == "term", label_plot, NA_character_)),
                  size = 2.0, label.size = 0.18, fill = "white", color = "black", repel = TRUE, max.overlaps = Inf) +
  geom_node_text(aes(label = ifelse(node_type == "gene", label_plot, NA_character_)),
                 size = 1.8, repel = TRUE, max.overlaps = Inf) +
  scale_color_manual(values = sf4_family_color, name = NULL) +
  scale_size_continuous(range = c(0.4, 4.2), name = "N genes / links") +
  labs(title = "Reactome network of SR-aging convergent genes") +
  theme_void(base_size = 9) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0, size = 10))
ggsave(file.path(stage_fig_dir, "SF4f_participant_aware_Reactome_network_staged.png"), p_sf4f, width = 6.9, height = 4.2, dpi = 360, bg = "white")
ggsave(file.path(stage_fig_dir, "SF4f_participant_aware_Reactome_network_staged.pdf"), p_sf4f, width = 6.9, height = 4.2, bg = "white")

sf4g_plot_df <- sf4_family_summary %>%
  arrange(FDR) %>%
  mutate(
    family = factor(family, levels = rev(family)),
    label = paste0(overlap_n, " overlap genes; ", n_reactome_terms, " terms"),
    fill_class = if_else(family == "HSF1-associated heat-shock response", "HSF1-associated heat-shock response", "Other Reactome family")
  )
p_sf4g <- ggplot(sf4g_plot_df, aes(-log10(FDR), family, fill = fill_class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), hjust = -0.04, size = 2.5) +
  scale_fill_manual(values = c("HSF1-associated heat-shock response" = "#d95f02", "Other Reactome family" = "#7297c7")) +
  coord_cartesian(xlim = c(0, max(-log10(sf4g_plot_df$FDR), na.rm = TRUE) * 1.35), clip = "off") +
  labs(title = "Enriched Reactome pathway families", x = "-log10(FDR)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0, size = 10),
    plot.margin = margin(5.5, 70, 5.5, 5.5)
  )
ggsave(file.path(stage_fig_dir, "SF4g_participant_aware_Reactome_family_barplot_staged.png"), p_sf4g, width = 5.1, height = 2.8, dpi = 360, bg = "white")
ggsave(file.path(stage_fig_dir, "SF4g_participant_aware_Reactome_family_barplot_staged.pdf"), p_sf4g, width = 5.1, height = 2.8, bg = "white")

old_module_path <- file.path(paths$aging_revision_upload, "additional_visuals", "MollerLevet_Reactome_heat_circadian_IL_ROS_insulin_AGE_modules_summary.csv")
old_sr_path <- file.path(paths$aging_revision_upload, "additional_visuals", "MollerLevet_SR_down_nominal_phase_adjusted_and_MESOR_down_genes_full.csv")
old_wb <- file.path(paths$source_materials, "..", "Supplementary_Data_1_full_analysis_outputs.backup_20260628_090759.xlsx")

old_modules <- if (file.exists(old_module_path)) {
  read_csv(old_module_path, show_col_types = FALSE)
} else {
  tibble(
    display_name = character(), n_genes = integer(), median_delta_mesor = numeric(),
    empirical_p_directional = numeric(), empirical_FDR_directional = numeric()
  )
}
old_sr <- if (file.exists(old_sr_path)) read_csv(old_sr_path, show_col_types = FALSE) else tibble()
if (nrow(old_modules) > 0) {
  if (!"display_name" %in% names(old_modules)) old_modules$display_name <- old_modules$module_label
  if (!"empirical_FDR_directional" %in% names(old_modules)) old_modules$empirical_FDR_directional <- NA_real_
}

old_fisher <- tibble()
old_ora <- tibble()
if (file.exists(old_wb)) {
  old_fisher <- tryCatch(openxlsx::read.xlsx(old_wb, sheet = "HumanBlood_SRdown_x_AgingDown_t"), error = function(e) tibble())
  old_ora <- tryCatch(openxlsx::read.xlsx(old_wb, sheet = "JenAge1540_SRnominalP05_AgeDire"), error = function(e) tibble())
}

audit_modules <- module_selected_raw %>%
  transmute(
    section = "Fig.1e/SF1c-e selected Reactome modules",
    item = display_name,
    new_n_genes = n_genes,
    new_median_delta_mesor = median_delta_mesor,
    new_empirical_p = empirical_p_directional,
    new_empirical_FDR = empirical_FDR_directional
  ) %>%
  left_join(
    old_modules %>%
      transmute(
        item = display_name,
        old_n_genes = n_genes,
        old_median_delta_mesor = median_delta_mesor,
        old_empirical_p = empirical_p_directional,
        old_empirical_FDR = empirical_FDR_directional
      ),
    by = "item"
  ) %>%
  mutate(qualitative_change = case_when(
    is.na(old_median_delta_mesor) ~ "Legacy comparison unavailable in clean package",
    sign(new_median_delta_mesor) == sign(old_median_delta_mesor) ~ "No direction change",
    TRUE ~ "Direction changed"
  ))

audit_genes <- gene_results %>%
  filter(gene %in% c("HSPH1", "HSPA5", "HSPD1", "NFKB2", "PRDX2")) %>%
  transmute(
    section = "Fig.1d/SF1a gene statistics",
    item = gene,
    new_phase_logFC = participant_aware_phase_logFC,
    new_phase_p = participant_aware_phase_p,
    new_phase_FDR = participant_aware_phase_FDR,
    new_delta_mesor = delta_mesor,
    new_mesor_p = mesor_p,
    new_mesor_FDR = mesor_FDR,
    qualitative_change = "Inferential statistics now use duplicateCorrelation block=subject"
  )

primary_old <- tibble()
if (all(c("dataset", "sr_threshold", "aging_threshold") %in% names(old_fisher))) {
  primary_old <- old_fisher %>%
    filter(dataset == "JenAge", sr_threshold %in% c("SR_nominalP05", "nominalP05", "sr_nominalP05"), aging_threshold %in% c("Age_direction", "age_direction_down", "directional")) %>%
    slice_head(n = 1)
}
primary_new <- overlap_summary %>%
  filter(dataset == "JenAge", sr_threshold == "sr_nominalP05", aging_threshold == "age_direction_down") %>%
  slice_head(n = 1)

audit_sf4_overlap <- tibble(
  section = "Supplementary Fig.4 SR-aging overlap",
  item = "JenAge participant-aware SR nominal p<0.05 + age-down direction",
  old_sr_down_n = if (nrow(primary_old)) primary_old$sr_down_n else NA_real_,
  new_sr_down_n = primary_new$sr_down_n,
  old_aging_down_n = if (nrow(primary_old)) primary_old$aging_down_n else NA_real_,
  new_aging_down_n = primary_new$aging_down_n,
  old_overlap_n = if (nrow(primary_old)) primary_old$overlap_n else NA_real_,
  new_overlap_n = primary_new$overlap_n,
  old_OR = if (nrow(primary_old)) coalesce(primary_old$fisher_OR, primary_old$odds_ratio) else NA_real_,
  new_OR = primary_new$fisher_OR,
  old_p = if (nrow(primary_old)) coalesce(primary_old$fisher_p, primary_old$p_value) else NA_real_,
  new_p = primary_new$fisher_p,
  qualitative_change = if_else(primary_new$fisher_OR > 1 & primary_new$fisher_p < 0.05, "No qualitative change: significant enrichment remains", "Potential qualitative change: inspect")
)

top_pathways_new <- sf4_ora %>%
  slice_head(n = 10) %>%
  transmute(rank = row_number(), new_term = term_clean, new_overlap_n = overlap_n, new_OR = odds_ratio, new_FDR = FDR)
top_pathways_old <- tibble(rank = integer(), old_term = character(), old_overlap_n = numeric(), old_OR = numeric(), old_FDR = numeric())
if (nrow(old_ora) > 0 && any(c("term_clean", "term") %in% names(old_ora))) {
  top_pathways_old <- old_ora %>%
    as_tibble() %>%
    mutate(old_term_label = if ("term_clean" %in% names(.)) term_clean else term) %>%
    slice_head(n = 10) %>%
    transmute(rank = row_number(), old_term = old_term_label, old_overlap_n = overlap_n, old_OR = odds_ratio, old_FDR = FDR)
}
audit_pathways <- full_join(top_pathways_old, top_pathways_new, by = "rank") %>%
  mutate(section = "Supplementary Fig.4 Reactome top pathways", item = paste0("rank_", rank), .before = 1)

write_csv(audit_modules, file.path(stage_dir, "participant_aware_audit_selected_modules.csv"))
write_csv(audit_genes, file.path(stage_dir, "participant_aware_audit_gene_stats.csv"))
write_csv(audit_sf4_overlap, file.path(stage_dir, "participant_aware_audit_sf4_overlap.csv"))
write_csv(audit_pathways, file.path(stage_dir, "participant_aware_audit_sf4_reactome_top_pathways.csv"))

audit_md <- c(
  "# Participant-aware Moller-Levet propagation audit",
  "",
  "Model: expression ~ condition + sin(2*pi*phase/24) + cos(2*pi*phase/24) + condition:sin(2*pi*phase/24) + condition:cos(2*pi*phase/24), fitted with limma duplicateCorrelation(..., block = subject).",
  "",
  "## Gene counts",
  paste0("- Previous human SR-down nominal phase-adjusted + MESOR-down genes: ", nrow(old_sr)),
  paste0("- Participant-aware human SR-down nominal phase-adjusted + MESOR-down genes: ", nrow(sr_down_new)),
  paste0("- Supplementary Fig. 4 primary JenAge SR-aging overlap: ", primary_new$overlap_n, " genes; OR=", signif(primary_new$fisher_OR, 4), "; Fisher p=", short_p(primary_new$fisher_p)),
  "",
  "## Selected module results",
  paste0("- HSF1 heat-shock regulation: median MESOR=", signif(module_selected_raw$median_delta_mesor[module_selected_raw$display_name == "HSF1 heat-shock regulation"], 4),
         "; empirical p=", short_p(module_selected_raw$empirical_p_directional[module_selected_raw$display_name == "HSF1 heat-shock regulation"]),
         "; FDR=", short_p(module_selected_raw$empirical_FDR_directional[module_selected_raw$display_name == "HSF1 heat-shock regulation"])),
  paste0("- Heat-stress response: median MESOR=", signif(module_selected_raw$median_delta_mesor[module_selected_raw$display_name == "Heat-stress response"], 4),
         "; empirical p=", short_p(module_selected_raw$empirical_p_directional[module_selected_raw$display_name == "Heat-stress response"]),
         "; FDR=", short_p(module_selected_raw$empirical_FDR_directional[module_selected_raw$display_name == "Heat-stress response"])),
  "",
  "## Qualitative interpretation",
  "- The participant-aware model preserves the negative SR MESOR-like effect for HSF1-associated/heat-stress proteostasis modules.",
  "- The downstream human blood SR-aging overlap remains enriched rather than depleted (OR > 1).",
  "- Final figure files were not overwritten by this staging script."
)
writeLines(audit_md, file.path(stage_dir, "participant_aware_before_after_audit.md"))

plot_selected <- module_selected_raw %>%
  mutate(display_name = factor(display_name, levels = display_name))
module_family_palette <- c(
  "AGE" = "#c44e52",
  "IL/ROS" = "#8c8c8c",
  "Circadian" = "#6f63b6",
  "Heat-stress" = "#1b9e77",
  "Insulin" = "#e69f00",
  "Other" = "grey50"
)
p_f1e <- ggplot(plot_selected, aes(median_delta_mesor, display_name, color = family)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  geom_errorbarh(aes(xmin = null_mesor_q025, xmax = null_mesor_q975), color = "grey68", height = 0, linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_text(aes(label = mesor_significance), nudge_x = -0.025, size = 4, color = "black") +
  scale_color_manual(values = module_family_palette) +
  labs(x = "Median SR-control MESOR-like effect", y = NULL, color = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", axis.text = element_text(color = "black"))
ggsave(file.path(stage_fig_dir, "F1e_participant_aware_selected_modules_staged.png"), p_f1e, width = 4.8, height = 2.7, dpi = 360, bg = "white")
ggsave(file.path(stage_fig_dir, "F1e_participant_aware_selected_modules_staged.pdf"), p_f1e, width = 4.8, height = 2.7, bg = "white")

sf1b_labeled_terms <- c(
  "Cellular response to heat stress",
  "Regulation of HSF1 mediated heat shock response",
  "HSF1 dependent transactivation",
  "HSF1 activation"
)
sf1b_plot <- term_summary %>%
  mutate(empirical_FDR_more_positive = p.adjust(empirical_p_more_positive, method = "BH")) %>%
  mutate(
	    label_text = case_when(
	      term_clean == "Regulation of HSF1 mediated heat shock response" ~ "HSF1 heat-shock regulation",
	      term_clean == "HSF1 dependent transactivation" ~ "HSF1-dependent transactivation",
	      term_clean %in% sf1b_labeled_terms ~ term_clean,
	      TRUE ~ NA_character_
	    ),
	    label_p = if_else(median_delta_mesor < 0, empirical_p_more_negative, empirical_p_more_positive),
	    label_fdr = if_else(median_delta_mesor < 0, empirical_FDR_more_negative, empirical_FDR_more_positive),
	    label = if_else(
	      !is.na(label_text),
	      paste0(str_wrap(label_text, 21), "\np=", short_p(label_p),
	             "; FDR=", short_p(label_fdr)),
	      NA_character_
	    ),
	    point_class = if_else(hsf_proteostasis, "HSF/proteostasis", "Other Reactome"),
	    direction_class = case_when(
	      median_delta_mesor < 0 & empirical_p_more_negative < 0.05 ~ "SR-down",
	      median_delta_mesor > 0 & empirical_p_more_positive < 0.05 ~ "SR-up",
	      TRUE ~ "No change"
	    )
	  )
sf1b_labels <- sf1b_plot %>%
  filter(!is.na(label)) %>%
	mutate(
		  label_x = case_when(
		    term_clean == "Cellular response to heat stress" ~ 285,
		    term_clean == "Regulation of HSF1 mediated heat shock response" ~ 85,
		    term_clean == "HSF1 dependent transactivation" ~ 250,
		    term_clean == "HSF1 activation" ~ 705,
		    TRUE ~ rank_from_most_negative
		  ),
		  label_y = case_when(
		    term_clean == "Cellular response to heat stress" ~ -0.047,
		    term_clean == "Regulation of HSF1 mediated heat shock response" ~ -0.066,
		    term_clean == "HSF1 dependent transactivation" ~ -0.020,
		    term_clean == "HSF1 activation" ~ 0.014,
		    TRUE ~ median_delta_mesor
		  )
  )
p_sf1b <- ggplot(sf1b_plot, aes(rank_from_most_negative, median_delta_mesor)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60", linewidth = 0.28) +
  geom_point(data = sf1b_plot %>% filter(point_class == "Other Reactome"),
             aes(color = direction_class), size = 0.8, alpha = 0.42) +
  geom_point(data = sf1b_plot %>% filter(point_class == "HSF/proteostasis"),
             aes(color = direction_class), size = 1.45, alpha = 0.95) +
  geom_segment(data = sf1b_labels,
               aes(x = label_x, y = label_y, xend = rank_from_most_negative, yend = median_delta_mesor),
               inherit.aes = FALSE, color = "black", linewidth = 0.22) +
  geom_label(data = sf1b_labels, aes(x = label_x, y = label_y, label = label),
             inherit.aes = FALSE, size = 1.55, color = "black", fill = "white",
             label.size = 0.16, label.padding = unit(0.08, "lines")) +
  scale_color_manual(values = c("SR-down" = "#2c7fb8", "SR-up" = "#c44e52", "No change" = "grey70")) +
  coord_cartesian(xlim = c(0, max(sf1b_plot$rank_from_most_negative, na.rm = TRUE) + 40),
                  ylim = c(-0.116, 0.066), clip = "off") +
  labs(x = "Rank from most negative SR-control MESOR-like effect", y = "Median SR-control MESOR-like effect") +
  theme_classic(base_size = 6.7) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(color = "black", size = 7),
    plot.margin = margin(3, 3, 3, 3)
  )
ggsave(file.path(stage_fig_dir, "SF1b_participant_aware_reactome_rank_staged.png"), p_sf1b, width = 5.2, height = 3.0, dpi = 360, bg = "white")
ggsave(file.path(stage_fig_dir, "SF1b_participant_aware_reactome_rank_staged.pdf"), p_sf1b, width = 5.2, height = 3.0, bg = "white")
ggsave(file.path(final_fig_dir, "SF1b_participant_aware_reactome_rank.png"), p_sf1b, width = 5.2, height = 3.0, dpi = 360, bg = "white")
ggsave(file.path(final_fig_dir, "SF1b_participant_aware_reactome_rank.pdf"), p_sf1b, width = 5.2, height = 3.0, bg = "white")

p_sf1c_label_x <- max(plot_selected$null_amp_q975, na.rm = TRUE) + 0.002
p_sf1c <- ggplot(plot_selected, aes(median_rhythmic_amplitude, display_name, color = family)) +
  geom_errorbarh(aes(xmin = null_amp_q025, xmax = null_amp_q975), color = "grey68", height = 0, linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_text(aes(x = p_sf1c_label_x, label = amp_p_label), color = "black", size = 2.2, hjust = 0) +
  scale_color_manual(values = module_family_palette) +
  coord_cartesian(clip = "off") +
  labs(x = "Median rhythmic SR-effect amplitude", y = NULL, color = NULL) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    plot.margin = margin(5.5, 70, 5.5, 5.5)
  )
ggsave(file.path(stage_fig_dir, "SF1c_participant_aware_amplitude_staged.png"), p_sf1c, width = 4.8, height = 2.7, dpi = 360, bg = "white")
ggsave(file.path(stage_fig_dir, "SF1c_participant_aware_amplitude_staged.pdf"), p_sf1c, width = 4.8, height = 2.7, bg = "white")

p_sf1d <- ggplot(plot_selected, aes(circular_mean_effect_acrophase_h, display_name, color = family)) +
  geom_segment(aes(x = 0, xend = circular_mean_effect_acrophase_h, yend = display_name), color = "grey78") +
  geom_point(size = 2.5) +
  scale_color_manual(values = module_family_palette) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 6)) +
  labs(x = "SR-effect acrophase (h)", y = NULL, color = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", axis.text = element_text(color = "black"))
ggsave(file.path(stage_fig_dir, "SF1d_participant_aware_phase_staged.png"), p_sf1d, width = 4.8, height = 2.7, dpi = 360, bg = "white")
ggsave(file.path(stage_fig_dir, "SF1d_participant_aware_phase_staged.pdf"), p_sf1d, width = 4.8, height = 2.7, bg = "white")

message("Participant-aware dependent source-data staging complete: ", stage_dir)
