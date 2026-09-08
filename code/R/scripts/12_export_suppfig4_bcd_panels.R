# Export Supplementary Fig. 4b-d panels from the submission source-data workbook.

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
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

workbook_path <- file.path(project_root, "Supplementary_Data_1_full_analysis_outputs.xlsx")
stage_fig_dir <- file.path(paths$figures, "participant_aware_staging")
final_fig_dir <- file.path(project_root, "final_figure_exports")
dir.create(stage_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)

read_panel_sheet <- function(sheet) {
  raw <- read.xlsx(workbook_path, sheet = sheet, colNames = FALSE)
  header <- as.character(unlist(raw[3, ]))
  dat <- raw[-seq_len(3), , drop = FALSE]
  names(dat) <- header
  dat
}

short_p <- function(x) {
  x <- as.numeric(x)
  ifelse(
    is.na(x),
    "NA",
    ifelse(x <= 0, "<1e-300", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x)))
  )
}

sf4b <- read_panel_sheet("SF4b") %>%
  mutate(
    n_genes = as.numeric(n_genes),
    median_age_effect = as.numeric(median_age_effect),
    empirical_p_directional = as.numeric(empirical_p_directional),
    fdr_directional = as.numeric(fdr_directional),
    hsf_proteostasis = as.logical(hsf_proteostasis)
  ) %>%
  arrange(median_age_effect) %>%
  mutate(
    rank = row_number(),
    point_class = case_when(
      hsf_proteostasis ~ "HSF/proteostasis",
      median_age_effect < 0 ~ "Age-down",
      TRUE ~ "Age-up"
    ),
    label = if_else(term_label %in% c(
      "Cd22 mediated bcr regulation",
      "HSF1 heat-shock regulation",
      "Cellular response to heat stress",
      "Non canonical inflammasome activation",
      "Programmed cell death"
    ), paste0(str_wrap(term_label, 22), "\np=", short_p(empirical_p_directional), "; FDR=", short_p(fdr_directional)), NA_character_)
  )

p_sf4b <- ggplot(sf4b, aes(rank, median_age_effect)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_point(aes(color = point_class, size = -log10(pmax(empirical_p_directional, 1e-300))),
             alpha = 0.78) +
  ggrepel::geom_text_repel(aes(label = label), seed = 20260811, size = 2.2,
                           max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = c("HSF/proteostasis" = "#1b9e77", "Age-down" = "#2c7fb8", "Age-up" = "#c44e52")) +
  scale_size_continuous(range = c(0.3, 3.0), name = "-log10(p)") +
  labs(title = "Aging direction of Reactome-defined modules",
       x = "Reactome modules ranked by median blood-aging effect",
       y = "Median age log2FC per decade", color = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0, size = 10))

sf4c <- read_panel_sheet("SF4c") %>%
  mutate(
    n_genes = as.numeric(n_genes),
    median_age_effect = as.numeric(median_age_effect),
    empirical_p_directional = as.numeric(empirical_p_directional),
    fdr_directional = as.numeric(fdr_directional),
    row_label = factor(row_label, levels = rev(row_label)),
    point_class = case_when(
      family == "HSF1/heat-shock arm" ~ "HSF/proteostasis",
      median_age_effect < 0 ~ "Age-down",
      TRUE ~ "Age-up"
    )
  )

p_sf4c <- ggplot(sf4c, aes(median_age_effect, row_label, color = point_class)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_segment(aes(x = 0, xend = median_age_effect, yend = row_label), color = "grey72", linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("HSF/proteostasis" = "#1b9e77", "Age-down" = "#2c7fb8", "Age-up" = "#c44e52")) +
  labs(title = "Representative Reactome modules", x = "Median age log2FC per decade", y = NULL, color = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0, size = 10))

sf4d <- read_panel_sheet("SF4d") %>%
  mutate(
    sr_logFC = as.numeric(sr_logFC),
    sr_p = as.numeric(sr_p),
    sr_fdr = as.numeric(sr_fdr),
    age_logFC_decade = as.numeric(age_logFC_decade),
    age_p = as.numeric(age_p),
    age_fdr = as.numeric(age_fdr)
  ) %>%
  filter(sr_threshold == "SR_nominalP05", aging_threshold == "Age_direction") %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(age_logFC_decade)) %>%
  mutate(
    gene = factor(gene, levels = gene),
    age_class = if_else(age_logFC_decade < 0, "Age-down", "Not age-down"),
    sig_bin = cut(-log10(pmax(age_p, 1e-300)), breaks = c(-Inf, 1, 2, 3, Inf),
                  labels = c("<1", "1-2", "2-3", ">=3"))
  )

p_sf4d <- ggplot(sf4d, aes(gene, age_logFC_decade, fill = age_class, alpha = sig_bin)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.25) +
  scale_fill_manual(values = c("Age-down" = "#2c7fb8", "Not age-down" = "#8fd19e")) +
  scale_alpha_manual(values = c("<1" = 0.35, "1-2" = 0.55, "2-3" = 0.75, ">=3" = 1.0), drop = FALSE) +
  labs(title = "Aging direction of proteostasis genes",
       x = NULL, y = "Aging log2FC per decade", fill = NULL, alpha = "-log10(p)") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0, size = 10))

exports <- list(
  SF4b = list(plot = p_sf4b, width = 5.0, height = 3.0, name = "SF4b_human_blood_Reactome_module_ranking"),
  SF4c = list(plot = p_sf4c, width = 4.2, height = 2.4, name = "SF4c_human_blood_representative_Reactome_modules"),
  SF4d = list(plot = p_sf4d, width = 6.2, height = 3.0, name = "SF4d_human_blood_proteostasis_gene_aging_direction")
)

for (panel in names(exports)) {
  item <- exports[[panel]]
  for (dir_i in c(stage_fig_dir, final_fig_dir)) {
    ggsave(file.path(dir_i, paste0(item$name, ".png")), item$plot, width = item$width, height = item$height, dpi = 360, bg = "white")
    ggsave(file.path(dir_i, paste0(item$name, ".pdf")), item$plot, width = item$width, height = item$height, bg = "white")
  }
}

message("Exported Supplementary Fig. 4b-d panels to: ", final_fig_dir)
