# Export Supplementary Fig. 4b-g in the original compact style, with
# participant-aware Moller-Levet-dependent data propagated to panels e-g.

if (!exists("paths")) {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
  project_root_guess <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
  Sys.setenv(NPJAGING_PROJECT_ROOT = Sys.getenv("NPJAGING_PROJECT_ROOT", unset = project_root_guess))
  source(file.path(dirname(script_dir), "config.R"))
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(igraph)
  library(ggraph)
  library(stringr)
  library(tibble)
})

workbook_path <- file.path(project_root, "Supplementary_Data_1_full_analysis_outputs.xlsx")
stage_fig_dir <- file.path(paths$figures, "participant_aware_staging")
final_fig_dir <- file.path(project_root, "final_figure_exports")
stage_dir <- file.path(paths$tables, "participant_aware_staging")
dir.create(stage_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_sf <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1)
    )
}

read_panel_sheet <- function(sheet) {
  raw <- read.xlsx(workbook_path, sheet = sheet, colNames = FALSE)
  header <- as.character(unlist(raw[3, ]))
  dat <- raw[-seq_len(3), , drop = FALSE]
  names(dat) <- header
  dat
}

short_p <- function(x) {
  x <- as.numeric(x)
  ifelse(is.na(x), "NA", ifelse(x <= 0, "<1e-300", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))))
}

save_panel <- function(plot, name, width, height) {
  for (dir_i in c(stage_fig_dir, final_fig_dir)) {
    ggsave(file.path(dir_i, paste0(name, ".png")), plot, width = width, height = height, dpi = 450, bg = "white")
    ggsave(file.path(dir_i, paste0(name, ".pdf")), plot, width = width, height = height, bg = "white")
  }
}

# Supplementary Fig. 4b
sf4b <- read_panel_sheet("SF4b") %>%
  mutate(
    n_genes = as.numeric(n_genes),
    median_age_effect = as.numeric(median_age_effect),
    empirical_p_directional = as.numeric(empirical_p_directional),
    fdr_directional = as.numeric(fdr_directional),
    hsf_proteostasis = as.logical(hsf_proteostasis)
  ) %>%
  filter(is.finite(median_age_effect)) %>%
  arrange(median_age_effect) %>%
  mutate(
    rank = row_number(),
	    direction_class = case_when(
	      median_age_effect < 0 ~ "Age-down",
	      median_age_effect > 0 ~ "Age-up",
	      TRUE ~ "Other"
	    ),
	    highlight = case_when(
	      fdr_directional < 0.05 & hsf_proteostasis ~ "HSF/proteostasis",
	      fdr_directional < 0.05 & direction_class == "Age-down" ~ "Age-down",
	      fdr_directional < 0.05 & direction_class == "Age-up" ~ "Age-up",
	      TRUE ~ "Other"
	    )
	  )

sf4b_labels <- sf4b %>%
  filter(str_detect(term_label, regex("Cd22 mediated bcr regulation|Fceri mediated mapk activation|Regulation of hsf1 mediated heat shock response|HSF1 heat-shock regulation|Cellular response to heat stress|^Programmed cell death$|Non canonical inflammasome activation", ignore_case = TRUE))) %>%
  distinct(term_label, .keep_all = TRUE) %>%
  mutate(
    clean_label = case_when(
      str_detect(term_label, regex("Cd22 mediated", ignore_case = TRUE)) ~ "B-cell/Fc/complement module",
      str_detect(term_label, regex("Fceri mediated", ignore_case = TRUE)) ~ "Signal-transduction module",
      str_detect(term_label, regex("Regulation of hsf1 mediated heat shock response", ignore_case = TRUE)) ~ "HSF1 heat-shock regulation",
      TRUE ~ term_label
    ),
    label = paste0(str_wrap(clean_label, 24), "\np", if_else(empirical_p_directional <= 0, "<1e-300", paste0("=", short_p(empirical_p_directional))),
                   "; FDR", if_else(fdr_directional <= 0, "<1e-300", paste0("=", short_p(fdr_directional)))),
    label_x = case_when(
      str_detect(clean_label, "B-cell") ~ 280,
      str_detect(clean_label, "Signal") ~ 520,
      clean_label == "HSF1 heat-shock regulation" ~ 135,
      clean_label == "Cellular response to heat stress" ~ 300,
      clean_label == "Programmed cell death" ~ 860,
      str_detect(clean_label, "inflammasome") ~ 1140,
      TRUE ~ rank
    ),
    label_y = case_when(
      str_detect(clean_label, "B-cell") ~ -0.128,
      str_detect(clean_label, "Signal") ~ -0.073,
      clean_label == "HSF1 heat-shock regulation" ~ -0.044,
      clean_label == "Cellular response to heat stress" ~ 0.026,
      clean_label == "Programmed cell death" ~ 0.031,
      str_detect(clean_label, "inflammasome") ~ 0.058,
      TRUE ~ median_age_effect
    )
  )

p_sf4b <- ggplot(sf4b, aes(rank, median_age_effect)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey65", linewidth = 0.35) +
  geom_point(data = sf4b %>% filter(highlight == "Other"), color = "grey78", size = 0.9, alpha = 0.65) +
  geom_point(data = sf4b %>% filter(highlight != "Other"), aes(color = highlight), size = 1.7, alpha = 0.95) +
  geom_segment(data = sf4b_labels, aes(x = label_x, y = label_y, xend = rank, yend = median_age_effect),
               inherit.aes = FALSE, color = "black", linewidth = 0.25) +
  geom_label(data = sf4b_labels, aes(x = label_x, y = label_y, label = label),
             inherit.aes = FALSE, size = 2.15, color = "black", fill = "white",
             label.size = 0.18, label.padding = unit(0.09, "lines")) +
  scale_color_manual(values = c("Age-down" = "#2c7fb8", "Age-up" = "#c44e52", "HSF/proteostasis" = "#1b9e77")) +
  labs(title = "Aging direction of Reactome-defined modules",
       x = "Reactome modules ranked by median blood-aging effect",
       y = "Median age log2FC per decade") +
  coord_cartesian(ylim = c(-0.145, 0.065), clip = "off") +
  theme_sf(8) +
  theme(legend.position = "none", plot.margin = margin(5.5, 18, 5.5, 5.5))

# Supplementary Fig. 4c
sf4c <- read_panel_sheet("SF4c") %>%
  mutate(
    n_genes = as.numeric(n_genes),
    median_age_effect = as.numeric(median_age_effect),
    null_q025 = as.numeric(null_q025),
    null_q975 = as.numeric(null_q975),
    empirical_p_directional = as.numeric(empirical_p_directional),
    fdr_directional = as.numeric(fdr_directional),
    row_label = factor(row_label, levels = rev(c(
      "Interleukin/inflammatory module (n=14)",
      "Cell-death module (n=184)",
      "Cellular heat-stress response (n=92)",
      "HSF1 heat-shock regulation (n=76)",
      "Signal-transduction module (n=83)",
      "B-cell/Fc/complement module (n=58)"
    ))),
    point_class = case_when(
      family == "HSF1/heat-shock arm" ~ "HSF/proteostasis",
      median_age_effect < 0 ~ "Age-down",
      TRUE ~ "Age-up"
    )
  )

p_sf4c <- ggplot(sf4c, aes(median_age_effect, row_label, color = point_class)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey65", linewidth = 0.25) +
  geom_segment(aes(x = null_q025, xend = null_q975, yend = row_label),
               color = "grey70", linewidth = 1.15) +
  geom_point(aes(x = x, y = row_label), color = "grey70", size = 1.35, inherit.aes = FALSE,
             data = sf4c %>% mutate(x = 0)) +
  geom_point(size = 1.95) +
  scale_color_manual(values = c("Age-down" = "#2166ac", "Age-up" = "#2166ac", "HSF/proteostasis" = "#1b9e77")) +
  coord_cartesian(xlim = c(-0.13, 0.055), clip = "off") +
  labs(x = "Median age log2FC per decade", y = NULL) +
  theme_classic(base_size = 6.5) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = "black", size = 5.8),
    axis.text.x = element_text(color = "black", size = 6),
    axis.title.x = element_text(color = "black", size = 7),
    axis.title.y = element_blank(),
    plot.margin = margin(2.5, 2.5, 2.5, 2.5)
  )

# Supplementary Fig. 4d
old_gene_order <- c("HSPA1A", "HSPA1B", "HSPA6", "DNAJB5", "HSF1", "STIP1", "DNAJB1",
                    "AHSA1", "FKBP4", "BAG3", "HSPH1", "DNAJA1", "HSP90AB1",
                    "HSPE1", "HSP90AA1", "HSPA8", "HSPD1")
sf4d <- read_panel_sheet("SF4d") %>%
  mutate(
    age_logFC_decade = as.numeric(age_logFC_per_decade),
    age_p = as.numeric(age_p),
    plot_order = as.numeric(plot_order)
  ) %>%
  filter(dataset == "JenAge", gene %in% old_gene_order) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(age_logFC_decade), gene) %>%
  select(gene, age_logFC_decade, age_p)
plot_order <- sf4d$gene
sf4d <- sf4d %>%
  filter(gene %in% plot_order) %>%
  mutate(
    gene = factor(gene, levels = plot_order),
    age_class = if_else(age_logFC_decade < 0, "Age-down", "Not age-down"),
    sig_bin = cut(-log10(pmax(age_p, 1e-300)), breaks = c(-Inf, 1, 2, 3, Inf),
                  labels = c("1", "2", "3", ">3"))
  )

p_sf4d <- ggplot(sf4d, aes(gene, age_logFC_decade, fill = age_class, alpha = sig_bin)) +
  geom_col(width = 0.72, color = "grey35", linewidth = 0.15) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.25) +
  scale_fill_manual(values = c("Age-down" = "#2b8cbe", "Not age-down" = "#9bd3a9"), drop = FALSE) +
  scale_alpha_manual(values = c("1" = 0.45, "2" = 0.65, "3" = 0.85, ">3" = 1), drop = FALSE) +
  labs(title = "Aging direction of proteostasis genes", x = NULL, y = "Aging log2FC per decade",
       fill = NULL, alpha = "-log10(p)") +
  theme_sf(8) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 6.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# Supplementary Fig. 4e
sf4_overlap <- read.csv(file.path(stage_dir, "SF4e_participant_aware_SRdown_x_AgingDown_fisher_summary.csv")) %>%
  filter(dataset == "JenAge", sr_threshold == "sr_nominalP05", aging_threshold == "age_direction_down") %>%
  slice(1)
sr_only_n <- sf4_overlap$sr_down_n - sf4_overlap$overlap_n
age_only_n <- sf4_overlap$aging_down_n - sf4_overlap$overlap_n

p_sf4e <- ggplot() +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "#2c7fb8", fill = "#d9edf7", alpha = 0.75, linewidth = 0.8) +
  geom_circle(aes(x0 = 1.15, y0 = 0, r = 1), color = "#2ca25f", fill = "#dff0df", alpha = 0.75, linewidth = 0.8) +
  annotate("text", x = -0.46, y = 0, label = sr_only_n, size = 4.2) +
  annotate("text", x = 0.57, y = 0, label = sf4_overlap$overlap_n, size = 4.2) +
  annotate("text", x = 1.60, y = 0, label = age_only_n, size = 4.2) +
  annotate("text", x = -0.30, y = 1.20, label = paste0("SR-down\nnominal p<0.05\nn=", sf4_overlap$sr_down_n),
           size = 2.25, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 1.42, y = 1.20, label = paste0("Blood-aging down\nlogFC<0\nn=", sf4_overlap$aging_down_n),
           size = 2.25, fontface = "bold", lineheight = 0.9) +
  annotate("label", x = 0.57, y = -1.38,
           label = paste0("Fisher OR=", signif(sf4_overlap$fisher_OR, 3), "; p=", short_p(sf4_overlap$fisher_p)),
           size = 2.3, label.size = 0.22) +
  coord_fixed(xlim = c(-1.15, 2.3), ylim = c(-1.55, 1.45), clip = "off") +
  labs(title = "SR-aging convergent genes") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 9),
        plot.margin = margin(4, 4, 4, 4))

# Supplementary Fig. 4f
sf4_ora <- read.csv(file.path(stage_dir, "SF4f_SF4g_participant_aware_JenAge_SRnominalP05_AgeDirection_Reactome_ORA.csv"))
network_terms <- bind_rows(
  sf4_ora %>% filter(FDR < 0.05) %>% arrange(FDR, desc(overlap_n)) %>% slice_head(n = 9),
  sf4_ora %>% filter(FDR < 0.05, hsf_proteostasis_related)
) %>%
  distinct(term, .keep_all = TRUE) %>%
  mutate(
    term_node = paste0("TERM__", term),
	    term_label = case_when(
	      str_detect(term_clean, regex("HSF1", ignore_case = TRUE)) ~ "HSF1-associated\nheat-shock response",
	      str_detect(term_clean, regex("Cellular response to heat stress", ignore_case = TRUE)) ~ "Cellular response\nto heat stress",
	      str_detect(term_clean, regex("SUMO", ignore_case = TRUE)) ~ "SUMOylation programs",
	      str_detect(term_clean, regex("Processing of capped intron", ignore_case = TRUE)) ~ "mRNA processing /\nsplicing",
	      str_detect(term_clean, regex("^SnRNP Assembly$", ignore_case = TRUE)) ~ "SnRNP assembly",
	      str_detect(term_clean, regex("RNA|MRNA|TRNA|RRNA|SPLIC", ignore_case = TRUE)) ~ "RNA processing /\nsplicing",
	      str_detect(term_clean, regex("CELL CYCLE", ignore_case = TRUE)) ~ "Cell-cycle regulation",
	      str_detect(term_clean, regex("TP53|DNA DAMAGE", ignore_case = TRUE)) ~ "DNA damage / TP53\nregulation",
	      TRUE ~ str_wrap(term_clean, 18)
	    ),
	    suppress_term_label = str_detect(term_clean, regex(
	      "^M Phase$|^Mitotic Prometaphase$|^Mitotic Spindle Checkpoint$",
	      ignore_case = TRUE
	    )) | (duplicated(term_label) & !str_detect(term_clean, regex("Processing of capped intron", ignore_case = TRUE))),
	    term_class = if_else(hsf_proteostasis_related, "HSF1 heat-shock family", "Other pathway family")
	  )
edge_tbl <- network_terms %>%
  transmute(from = term_node, gene = str_split(overlap_gene_string, ";")) %>%
  unnest(gene) %>%
  filter(gene != "") %>%
  mutate(to = paste0("GENE__", gene))
hsf_family_genes <- network_terms %>%
  filter(hsf_proteostasis_related) %>%
  transmute(gene = str_split(overlap_gene_string, ";")) %>%
  unnest(gene) %>%
  filter(gene != "") %>%
  distinct(gene) %>%
  pull(gene)
hsf_family_label_genes <- hsf_family_genes[str_detect(
  hsf_family_genes,
  regex("^HSP|^DNAJ|^BAG", ignore_case = TRUE)
)]
hsf_family_label_genes <- setdiff(hsf_family_label_genes, c("HSPA4", "HSPA14", "BAG2"))
term_nodes <- network_terms %>%
  transmute(name = term_node, label = if_else(suppress_term_label, NA_character_, term_label),
            node_type = "term", node_class = term_class)
gene_nodes <- edge_tbl %>%
  distinct(name = to, gene) %>%
  mutate(
	    label = gene,
	    node_type = "gene",
	    node_class = case_when(
	      gene %in% hsf_family_genes ~ "HSF1 heat-shock family gene",
	      TRUE ~ "Shared convergence gene"
	    )
	  )
graph <- graph_from_data_frame(edge_tbl %>% select(from, to), vertices = bind_rows(term_nodes, gene_nodes), directed = FALSE)
deg <- degree(graph)
V(graph)$node_size <- ifelse(V(graph)$node_type == "term", 4.6, ifelse(V(graph)$node_class == "HSF1 heat-shock family gene", 1.9, 1.0))
V(graph)$label_plot <- ifelse(
  V(graph)$node_type == "term" |
    (V(graph)$node_class == "HSF1 heat-shock family gene" &
       V(graph)$label %in% hsf_family_label_genes &
       deg[V(graph)$name] >= 1),
  V(graph)$label,
  NA_character_
)

set.seed(41)
sf4f_layout <- create_layout(graph, layout = "fr")
sf4f_gene_labels <- as_tibble(sf4f_layout) %>%
  filter(node_class == "HSF1 heat-shock family gene", !is.na(label_plot)) %>%
  mutate(
    x_label = if_else(label == "BAG2", x - 1.05, x),
    y_label = if_else(label == "BAG2", y - 0.48, y),
    nudge_x = if_else(label == "BAG2", -0.55, 0),
    nudge_y = if_else(label == "BAG2", -0.10, 0)
  )
p_sf4f <- ggraph(sf4f_layout) +
  geom_edge_link(color = "grey80", alpha = 0.55, linewidth = 0.2) +
  geom_node_point(aes(color = node_class, size = node_size), alpha = 0.95) +
  geom_node_label(aes(label = ifelse(node_type == "term", label_plot, NA_character_)),
                  color = "black", size = 1.85, label.size = 0.16, fill = "white",
                  repel = TRUE, box.padding = unit(0.32, "lines"),
                  point.padding = unit(0.25, "lines"), max.overlaps = Inf,
                  show.legend = FALSE) +
  geom_text_repel(data = sf4f_gene_labels,
                  aes(x = x_label, y = y_label, label = label_plot),
                  inherit.aes = FALSE, color = "black", size = 1.45,
                  fontface = "italic",
                  box.padding = unit(0.28, "lines"),
                  point.padding = unit(0.22, "lines"), segment.color = NA,
                  max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c("HSF1 heat-shock family gene" = "#1b9e77",
                                "HSF1 heat-shock family" = "#d95f02",
                                "Other pathway family" = "#6a51a3",
                                "Shared convergence gene" = "grey70"),
	                     breaks = c("HSF1 heat-shock family gene", "HSF1 heat-shock family", "Other pathway family", "Shared convergence gene")) +
  scale_size_identity(guide = "none") +
  labs(title = "Reactome network of SR-aging convergent genes", color = NULL) +
  theme_void(base_size = 8) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0, size = 9),
        legend.text = element_text(size = 6.5), plot.margin = margin(4, 4, 4, 4))

# Supplementary Fig. 4g
sf4g <- read.csv(file.path(stage_dir, "SF4g_participant_aware_Reactome_family_summary.csv")) %>%
  arrange(FDR) %>%
  mutate(
    family = factor(family, levels = rev(family)),
    label = paste0(overlap_n, " overlap genes; ", n_reactome_terms, " terms"),
    fill_class = if_else(family == "HSF1-associated heat-shock response", "HSF1-associated heat-shock response", "Other Reactome family")
  )

p_sf4g <- ggplot(sf4g, aes(-log10(FDR), family, fill = fill_class)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = label), hjust = -0.04, size = 2.35) +
  scale_fill_manual(values = c("HSF1-associated heat-shock response" = "#d95f02", "Other Reactome family" = "#7297c7")) +
  coord_cartesian(xlim = c(0, max(-log10(sf4g$FDR), na.rm = TRUE) * 1.38), clip = "off") +
  labs(title = "Enriched Reactome pathway families", x = "-log10(FDR)", y = NULL) +
  theme_sf(8) +
  theme(legend.position = "none", axis.text.y = element_text(size = 6.5),
        plot.margin = margin(5.5, 80, 5.5, 5.5))

save_panel(p_sf4b, "SF4b_final_style_human_blood_Reactome_module_ranking", 6.3, 3.05)
save_panel(p_sf4c, "SF4c_final_style_human_blood_representative_Reactome_modules", 3.45, 1.55)
save_panel(p_sf4d, "SF4d_final_style_human_blood_proteostasis_gene_aging_direction", 5.1, 2.65)
save_panel(p_sf4e, "SF4e_final_style_participant_aware_SR_aging_overlap", 2.85, 2.55)
save_panel(p_sf4f, "SF4f_final_style_participant_aware_Reactome_network", 7.1, 3.9)
save_panel(p_sf4g, "SF4g_final_style_participant_aware_Reactome_family_barplot", 4.85, 2.45)

message("Exported final-style Supplementary Fig. 4b-g panels to: ", final_fig_dir)
