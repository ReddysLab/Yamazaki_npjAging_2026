#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggforce)
})

if (!exists("paths", inherits = TRUE)) {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
  source(file.path(dirname(script_dir), "config.R"))
}

stage_dir <- file.path(paths$tables, "figure3_participant_aware")
export_dir <- file.path(paths$figures, "figure3_participant_aware")
gh_fig3 <- file.path(paths$tables, "final_source_data", "Figure3")
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gh_fig3, recursive = TRUE, showWarnings = FALSE)

h5ad <- Sys.getenv("NPJAGING_SCHAUM_H5AD", unset = file.path(paths$source_data, "GSE132040", "tabula_muris_senis.h5ad"))
pb_mat_path <- file.path(paths$curated_metadata, "GSE132040", "Figure3D_neuronal_binning_gene_slopes.csv")
pb_def_path <- file.path(paths$curated_metadata, "GSE132040", "Figure3D_neuronal_binning_group_definitions.csv")
sr_down_path <- file.path(paths$tables, "moller_levet_participant_aware", "MollerLevet_participant_aware_phase_cosinor_gene_results.csv")
sr_all_path <- sr_down_path
map_path <- file.path(paths$source_data, "GSE132040", "mouse_gene_map.csv")

if (!file.exists(h5ad)) {
  stop("Missing Schaum et al. GSE132040 h5ad file. Set NPJAGING_SCHAUM_H5AD or place it at source_data/GSE132040/tabula_muris_senis.h5ad.", call. = FALSE)
}

clean_label <- function(x) {
  x |>
    gsub("_", " ", x = _) |>
    gsub("\\s+", " ", x = _) |>
    trimws()
}

read_obs_factor <- function(field) {
  cats <- h5read(h5ad, paste0("/obs/", field, "/categories"))
  codes <- h5read(h5ad, paste0("/obs/", field, "/codes"))
  ifelse(codes >= 0, cats[codes + 1], NA_character_)
}

mean_sparse_rows <- function(cell_idx, label) {
  indptr <- h5read(h5ad, "/X/indptr")
  n_genes <- length(h5read(h5ad, "/var/gene_symbols"))
  sum_vec <- numeric(n_genes)
  for (i in cell_idx) {
    start0 <- indptr[i] + 1
    end0 <- indptr[i + 1]
    if (end0 >= start0) {
      len <- end0 - start0 + 1
      gene_idx <- h5read(h5ad, "/X/indices", index = list(seq.int(start0, end0))) + 1L
      vals <- h5read(h5ad, "/X/data", index = list(seq.int(start0, end0)))
      sum_vec[gene_idx] <- sum_vec[gene_idx] + vals
    }
  }
  sum_vec / length(cell_idx)
}

fit_slopes <- function(expr_by_age, ages, genes, group_name) {
  res <- t(apply(expr_by_age, 1, function(y) {
    fit <- lm(as.numeric(y) ~ ages)
    c(slope = unname(coef(fit)[2]), p_value = summary(fit)$coefficients[2, 4])
  }))
  tibble(
    gene_mouse = as.character(genes),
    cell_group = group_name,
    aging_slope = as.numeric(res[, "slope"]),
    aging_p = as.numeric(res[, "p_value"])
  )
}

representative_order <- c("Hsf1", "Hspa1a", "Hsp90aa1", "Hsp90ab1", "Hspa5", "Hspd1",
                          "Hsph1", "Dnajb1", "Dnajb6", "Bag3", "Psmb5", "Eif2ak3",
                          "Atf4", "Atf6", "Xbp1", "Ddit3")

bin_slopes <- read_csv(pb_mat_path, show_col_types = FALSE)
fig3d_long <- bin_slopes |>
  filter(gene %in% representative_order) |>
  mutate(
    panel_cell_type = recode(
      group,
      "Atlas neuron label" = "Excitatory",
      "GABA/inhibitory-like neuronal labels" = "Inhibitory",
      "Broad neurons" = "All neurons",
      "Astrocyte" = "Astrocyte",
      "Oligodendrocyte" = "Oligodend.",
      "Microglia" = "Microglia",
      "Endothelial" = "Endothelial",
      .default = NA_character_
    ),
    plotted_value = slope_expr_vs_age,
    panel = "Fig. 3d",
    plotted_value_definition = "Slope from pseudobulk expression ~ age, using the cell-type grouping shown in the figure."
  ) |>
  filter(!is.na(panel_cell_type)) |>
  mutate(
    gene = factor(gene, levels = representative_order),
    panel_cell_type = factor(panel_cell_type, levels = c("Excitatory", "Inhibitory", "All neurons", "Astrocyte", "Oligodend.", "Microglia", "Endothelial"))
  ) |>
  arrange(gene, panel_cell_type) |>
  mutate(gene = as.character(gene), panel_cell_type = as.character(panel_cell_type)) |>
  select(panel, gene, panel_cell_type, plotted_value, source_group = group, n_age_points, min_cells_per_age, gene_class, plotted_value_definition)

fig3d_wide <- fig3d_long |>
  select(gene, panel_cell_type, plotted_value) |>
  pivot_wider(names_from = panel_cell_type, values_from = plotted_value) |>
  mutate(gene = factor(gene, levels = representative_order)) |>
  arrange(gene) |>
  mutate(gene = as.character(gene))

fig3d_defs <- tibble(
  panel = "Fig. 3d",
  panel_label = c("Excitatory", "Inhibitory", "All neurons", "Astrocyte", "Oligodend.", "Microglia", "Endothelial"),
  source_group = c("Atlas neuron label", "GABA/inhibitory-like neuronal labels", "Broad neurons", "Astrocyte", "Oligodendrocyte", "Microglia", "Endothelial"),
  operational_definition = c(
    "Atlas cell_type label 'neuron'; used as the excitatory/non-interneuron neuronal group in the displayed heatmap.",
    "Pooled inhibitory/GABAergic-like neuronal labels, including interneuron and medium spiny neuron labels where available.",
    "All neuronal labels pooled for the broad neuron group.",
    "Atlas astrocyte label.",
    "Atlas oligodendrocyte label.",
    "Atlas microglial cell label.",
    "Atlas endothelial cell label."
  )
)

write_csv(fig3d_long, file.path(stage_dir, "Fig3d_cell_type_gene_slopes_long.csv"))
write_csv(fig3d_wide, file.path(stage_dir, "Fig3d_cell_type_gene_slopes_wide.csv"))
write_csv(fig3d_defs, file.path(stage_dir, "Fig3d_cell_type_group_definitions.csv"))
write_csv(fig3d_wide, file.path(gh_fig3, "Figure_3D_source_data.csv"))

tissue <- read_obs_factor("tissue")
cell_type <- read_obs_factor("cell_type")
age <- read_obs_factor("age")
genes <- as.character(h5read(h5ad, "/var/gene_symbols"))

obs <- tibble(row_index = seq_along(tissue), tissue, cell_type, age) |>
  mutate(age_num = as.numeric(gsub("m", "", age)))

groups <- list(
  `Excitatory neurons` = c("neuron"),
  `Inhibitory neurons` = c("interneuron", "medium spiny neuron")
)

all_slopes <- lapply(names(groups), function(group_name) {
  cts <- groups[[group_name]]
  idx_by_age <- lapply(sort(unique(obs$age_num[obs$tissue == "brain" & obs$cell_type %in% cts])), function(a) {
    obs$row_index[obs$tissue == "brain" & obs$cell_type %in% cts & obs$age_num == a]
  })
  ages <- sort(unique(obs$age_num[obs$tissue == "brain" & obs$cell_type %in% cts]))
  mat <- sapply(seq_along(idx_by_age), function(i) mean_sparse_rows(idx_by_age[[i]], group_name))
  colnames(mat) <- paste0(ages, "m")
  fit_slopes(mat, ages, genes, group_name) |>
    mutate(
      source_cell_type_labels = paste(cts, collapse = "; "),
      n_cells_3m = length(idx_by_age[[which(ages == 3)]]),
      n_cells_18m = length(idx_by_age[[which(ages == 18)]]),
      n_cells_24m = length(idx_by_age[[which(ages == 24)]])
    )
}) |>
  bind_rows()

gene_map <- read_csv(map_path, show_col_types = FALSE) |>
  mutate(gene_upper = toupper(gene_upper)) |>
  distinct(gene_mouse, gene_upper)

all_slopes_mapped <- all_slopes |>
  left_join(gene_map, by = "gene_mouse") |>
  mutate(gene_upper = if_else(is.na(gene_upper), toupper(gene_mouse), gene_upper))

sr_down <- read_csv(sr_down_path, show_col_types = FALSE) |>
  transmute(gene_upper = toupper(gene), participant_aware_phase_logFC, participant_aware_phase_p, participant_aware_phase_FDR, delta_mesor, mesor_p, mesor_FDR) |>
  filter(participant_aware_phase_logFC < 0, participant_aware_phase_FDR < 0.05,
         delta_mesor < 0, mesor_FDR < 0.05) |>
  distinct(gene_upper, .keep_all = TRUE)

sr_all <- read_csv(sr_all_path, show_col_types = FALSE) |>
  transmute(gene_upper = toupper(gene)) |>
  distinct(gene_upper)

fisher_for_group <- function(group_name) {
  neuron_tbl <- all_slopes_mapped |>
    filter(cell_group == group_name) |>
    distinct(gene_upper, .keep_all = TRUE) |>
    mutate(neuron_aging_down = aging_slope < 0 & aging_p < 0.05)
  universe <- intersect(sr_all$gene_upper, neuron_tbl$gene_upper)
  sr_set <- intersect(sr_down$gene_upper, universe)
  neuron_set <- neuron_tbl |>
    filter(gene_upper %in% universe, neuron_aging_down) |>
    pull(gene_upper)
  overlap <- intersect(sr_set, neuron_set)
  m <- matrix(c(
    length(overlap),
    length(setdiff(sr_set, neuron_set)),
    length(setdiff(neuron_set, sr_set)),
    length(setdiff(universe, union(sr_set, neuron_set)))
  ), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m)
  tibble(
    panel = "Fig. 3e",
    comparison = group_name,
    sr_set_label = "SR-down human blood",
    neuron_set_label = paste0(sub(" neurons", "", group_name), " neuron aging-down"),
    universe_n = length(universe),
    sr_down_n = length(sr_set),
    neuron_aging_down_n = length(neuron_set),
    overlap_n = length(overlap),
    sr_only_n = length(setdiff(sr_set, neuron_set)),
    neuron_only_n = length(setdiff(neuron_set, sr_set)),
    neither_n = length(setdiff(universe, union(sr_set, neuron_set))),
    fisher_OR = unname(ft$estimate),
    fisher_p = ft$p.value,
    sr_threshold = "Participant-aware Moller-Levet SR-down: phase-adjusted logFC < 0 and BH FDR < 0.05; delta MESOR < 0 and BH FDR < 0.05.",
    neuron_threshold = "Schaum et al. brain neuronal pseudobulk aging-down: expression ~ age slope < 0 and nominal p < 0.05.",
    background_universe = "Genes present in both the participant-aware human sleep-restriction model and the mapped Schaum et al. neuronal aging pseudobulk table."
  )
}

fig3e_summary <- bind_rows(lapply(names(groups), fisher_for_group))

fig3e_membership <- all_slopes_mapped |>
  filter(cell_group %in% names(groups)) |>
  distinct(cell_group, gene_upper, .keep_all = TRUE) |>
  select(cell_group, source_cell_type_labels, gene_mouse, gene_upper, aging_slope, aging_p, starts_with("n_cells")) |>
  left_join(sr_down |> mutate(participant_aware_SR_down = TRUE), by = "gene_upper") |>
  mutate(
    participant_aware_SR_down = if_else(is.na(participant_aware_SR_down), FALSE, participant_aware_SR_down),
    neuron_aging_down = aging_slope < 0 & aging_p < 0.05,
    in_fig3e_overlap = participant_aware_SR_down & neuron_aging_down
  ) |>
  filter(gene_upper %in% sr_all$gene_upper) |>
  arrange(cell_group, desc(in_fig3e_overlap), gene_upper)

fig3e_overlap_genes <- fig3e_membership |>
  filter(in_fig3e_overlap) |>
  select(cell_group, source_cell_type_labels, gene_mouse, gene_upper, aging_slope, aging_p, participant_aware_phase_logFC, participant_aware_phase_p, participant_aware_phase_FDR, delta_mesor, mesor_p, mesor_FDR)

write_csv(all_slopes_mapped, file.path(stage_dir, "Fig3e_Schaum_brain_neuron_subtype_aging_slopes_all_genes.csv"))
write_csv(fig3e_summary, file.path(stage_dir, "Fig3e_participant_aware_SRdown_neuronal_agingdown_fisher_summary.csv"))
write_csv(fig3e_membership, file.path(stage_dir, "Fig3e_participant_aware_gene_membership_full.csv"))
write_csv(fig3e_overlap_genes, file.path(stage_dir, "Fig3e_participant_aware_overlap_genes.csv"))
write_csv(fig3e_summary, file.path(gh_fig3, "Figure_3E_participant_aware_venn_fisher_summary.csv"))
write_csv(fig3e_membership, file.path(gh_fig3, "Figure_3E_participant_aware_gene_membership_full.csv"))
write_csv(fig3e_overlap_genes, file.path(gh_fig3, "Figure_3E_participant_aware_overlap_genes.csv"))

fmt_p <- function(p) {
  ifelse(p < 1e-4, "p<1e-4", paste0("p=", signif(p, 2)))
}

venn_df <- fig3e_summary |>
  mutate(
    x_left = if_else(comparison == "Excitatory neurons", -1.15, 2.15),
    x_right = if_else(comparison == "Excitatory neurons", -0.25, 3.05),
    label_x = if_else(comparison == "Excitatory neurons", -0.70, 2.60),
    heading_x = if_else(comparison == "Excitatory neurons", -0.70, 2.60),
    panel_title = comparison
  )

label_df <- bind_rows(
  venn_df |> transmute(comparison, x = x_left - 0.18, y = 1.12, label = paste0("SR-down\nhuman blood\nn=", sr_down_n)),
  venn_df |> transmute(comparison, x = x_right + 0.18, y = 1.12, label = paste0(sub(" neurons", "", comparison), "\naging-down\nslope < 0; p<0.05\nn=", neuron_aging_down_n)),
  venn_df |> transmute(comparison, x = x_left - 0.24, y = 0.0, label = as.character(sr_only_n)),
  venn_df |> transmute(comparison, x = label_x, y = 0.0, label = as.character(overlap_n)),
  venn_df |> transmute(comparison, x = x_right + 0.24, y = 0.0, label = as.character(neuron_only_n))
)

stat_df <- venn_df |>
  transmute(
    comparison,
    x = label_x,
    y = -0.92,
    label = paste0("Fisher OR=", sprintf("%.2f", fisher_OR), "; ", fmt_p(fisher_p))
  )

p <- ggplot() +
  geom_circle(data = venn_df, aes(x0 = x_left, y0 = 0, r = 0.78), fill = "#dceef9", color = "#2b83c6", linewidth = 0.7, alpha = 0.75) +
  geom_circle(data = venn_df, aes(x0 = x_right, y0 = 0, r = 0.78), fill = "#e0f2df", color = "#2ca25f", linewidth = 0.7, alpha = 0.75) +
  geom_text(data = label_df, aes(x, y, label = label), size = 3.1, lineheight = 0.9) +
  geom_label(data = stat_df, aes(x, y, label = label), size = 2.35, label.size = 0.25, fill = "white", label.padding = unit(0.10, "lines")) +
  geom_text(data = venn_df, aes(heading_x, 1.68, label = panel_title), fontface = "bold", size = 3.25) +
  coord_fixed(xlim = c(-2.25, 4.15), ylim = c(-1.18, 1.90), clip = "off") +
  theme_void(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(6, 6, 6, 6, "pt")
  )

ggsave(file.path(export_dir, "F3e_participant_aware_SRdown_neuronal_agingdown_venn.pdf"), p, width = 6.2, height = 2.5)
ggsave(file.path(export_dir, "F3e_participant_aware_SRdown_neuronal_agingdown_venn.png"), p, width = 6.2, height = 2.5, dpi = 600)

message("Fig. 3d/e participant-aware outputs written to: ", stage_dir)
print(fig3e_summary)
