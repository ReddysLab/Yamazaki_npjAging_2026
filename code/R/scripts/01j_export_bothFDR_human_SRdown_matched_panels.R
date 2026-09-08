#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Mm.eg.db)
  library(msigdbr)
})

if (!exists("paths", inherits = TRUE)) {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
  source(file.path(dirname(script_dir), "config.R"))
}

tables <- paths$tables
out <- paths$final_dual_fdr
unlink(out, recursive = TRUE, force = TRUE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

threshold_text <- paste(
  "participant-aware phase-adjusted logFC < 0; phase-model FDR < 0.05;",
  "delta MESOR < 0; MESOR FDR < 0.05"
)

clean_term <- function(x) {
  x |>
    str_remove("^REACTOME_") |>
    str_replace_all("_", " ") |>
    str_squish() |>
    str_to_sentence()
}

hsf_flag <- function(x) {
  str_detect(x, regex("HSF1|HEAT SHOCK|HEAT STRESS|HSP90|CHAPERONE|UNFOLDED PROTEIN|PROTEOSTASIS", ignore_case = TRUE))
}

short_label <- function(x) {
  x |>
    str_replace(regex("Regulation of HSF1.?mediated heat shock response", ignore_case = TRUE), "HSF1 heat-shock regulation") |>
    str_replace(regex("Cellular response to heat stress", ignore_case = TRUE), "Cellular heat-stress response") |>
    str_replace(regex("Processing of capped intron containing pre mrna", ignore_case = TRUE), "RNA processing / splicing") |>
    str_replace(regex("Transport of mature transcript to cytoplasm", ignore_case = TRUE), "Nuclear RNA export") |>
    str_replace(regex("Sumoylation of sumoylation proteins", ignore_case = TRUE), "SUMOylation programs") |>
    str_replace(regex("Sumoylation of dna replication proteins", ignore_case = TRUE), "SUMOylation of DNA-replication proteins") |>
    str_replace(regex("^Hsp90 chaperone cycle.*", ignore_case = TRUE), "HSP90 chaperone cycle") |>
    str_replace(regex("^Hsf1.?dependent transactivation$", ignore_case = TRUE), "HSF1-dependent transactivation") |>
    str_replace(regex("Chaperone mediated autophagy", ignore_case = TRUE), "Chaperone-mediated autophagy") |>
    str_replace(regex("Late endosomal microautophagy", ignore_case = TRUE), "Late endosomal microautophagy") |>
    str_replace(regex("Dna damage recognition in gg ner", ignore_case = TRUE), "DNA-damage recognition (GG-NER)") |>
    str_replace(regex("Global genome nucleotide excision repair gg ner", ignore_case = TRUE), "Global-genome NER") |>
    str_replace(regex("Fanconi anemia pathway", ignore_case = TRUE), "Fanconi anemia pathway") |>
    str_replace(regex("Mismatch repair", ignore_case = TRUE), "Mismatch repair") |>
    str_replace(regex("Mapk1 mapk3 signaling", ignore_case = TRUE), "MAPK1/MAPK3 signaling") |>
    str_replace(regex("Rhoa gtpase cycle", ignore_case = TRUE), "RHOA GTPase cycle") |>
    str_replace(regex("Auf1 hnnp d0 binds and destabilizes mrna", ignore_case = TRUE), "AUF1/hnRNP D0 mRNA destabilization") |>
    str_replace(regex("^N.?glycan trimming.*calnexin.*calreticulin.*$", ignore_case = TRUE), "N-glycan trimming / calnexin-calreticulin") |>
    str_replace(regex("^Antigen presentation.*peptide loading.*class i mhc$", ignore_case = TRUE), "MHC-I antigen presentation / peptide loading") |>
    str_replace(regex("^Metabolism of rna$", ignore_case = TRUE), "RNA metabolism") |>
    str_replace(regex("^Post-translational protein modification$", ignore_case = TRUE), "Post-translational protein modification") |>
    str_replace(regex("^Regulation of pd-l1\\(cd274\\) post-translational modification$", ignore_case = TRUE), "PD-L1 post-translational regulation") |>
    str_replace(regex("^Metabolism of water-soluble vitamins and cofactors$", ignore_case = TRUE), "Water-soluble vitamin / cofactor metabolism") |>
    str_replace(regex("^Oxidative stress induced senescence$", ignore_case = TRUE), "Oxidative stress-induced senescence") |>
    str_replace(regex("^Class i mhc mediated antigen processing & presentation$", ignore_case = TRUE), "Class I MHC antigen processing / presentation") |>
    str_replace(regex("^Activation of nmda receptors and postsynaptic events$", ignore_case = TRUE), "NMDA-receptor postsynaptic signaling") |>
    str_replace(regex("^Major pathway of rrna processing in the nucleolus and cytosol$", ignore_case = TRUE), "rRNA processing") |>
    str_replace(regex("^Sumoylation of ubiquitinylation proteins$", ignore_case = TRUE), "SUMOylation of ubiquitinylation proteins") |>
    str_replace(regex("Regulation of endogenous retroelements", ignore_case = TRUE), "Endogenous retroelement regulation") |>
    str_wrap(width = 24)
}

term_family <- function(x) {
  case_when(
    hsf_flag(x) ~ "HSF/proteostasis",
    str_detect(x, regex("RNA|MRNA|TRNA|RRNA|SPLIC|RIBOSOM|TRANSLAT|NUCLEAR (IMPORT|EXPORT|PORE)|SUMO", ignore_case = TRUE)) ~ "RNA processing / nuclear transport",
    str_detect(x, regex("CELL CYCLE|MITOT|CHROMAT|DNA|TP53|RHO.?GTPASE", ignore_case = TRUE)) ~ "Cell cycle / chromatin",
    str_detect(x, regex("NEUR|NTRK|NGF|SYNAP|SIGNAL|KINASE|GTPASE", ignore_case = TRUE)) ~ "Signaling / neuronal",
    str_detect(x, regex("ER |ENDOPLASMIC|CALNEXIN|CALRETICULIN|GLYCOSYL|VESIC|GOLGI|CARGO|PROTEIN TRANSPORT", ignore_case = TRUE)) ~ "ER / protein transport",
    str_detect(x, regex("MITOCHON|RESPIR|METABOL|PYRUVATE|OXID", ignore_case = TRUE)) ~ "Metabolism / mitochondria",
    TRUE ~ "Other"
  )
}

fisher_summary <- function(set_a, set_b, universe) {
  set_a <- intersect(unique(set_a), universe)
  set_b <- intersect(unique(set_b), universe)
  mat <- matrix(c(
    length(intersect(set_a, set_b)),
    length(setdiff(set_a, set_b)),
    length(setdiff(set_b, set_a)),
    length(setdiff(universe, union(set_a, set_b)))
  ), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat, alternative = "greater")
  tibble(
    universe_n = length(universe),
    set_a_n = length(set_a),
    set_b_n = length(set_b),
    overlap_n = length(intersect(set_a, set_b)),
    set_a_only_n = length(setdiff(set_a, set_b)),
    set_b_only_n = length(setdiff(set_b, set_a)),
    fisher_OR = unname(ft$estimate),
    fisher_p = ft$p.value
  )
}

plot_venn <- function(stats, left_label, right_label, caption, stem, width = 3.3, height = 2.7) {
  p_label <- if (stats$fisher_p > 0.999) {
    sprintf("%.4f", stats$fisher_p)
  } else if (stats$fisher_p >= 0.001) {
    sprintf("%.3f", stats$fisher_p)
  } else {
    format(stats$fisher_p, digits = 2, scientific = TRUE)
  }
  p <- ggplot() +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "#2c7fb8", fill = "#d9edf7", alpha = 0.75, linewidth = 0.8) +
    geom_circle(aes(x0 = 1.15, y0 = 0, r = 1), color = "#2ca25f", fill = "#dff0df", alpha = 0.75, linewidth = 0.8) +
    annotate("text", x = -0.46, y = 0, label = stats$set_a_only_n, size = 4.2) +
    annotate("text", x = 0.57, y = 0, label = stats$overlap_n, size = 4.2) +
    annotate("text", x = 1.60, y = 0, label = stats$set_b_only_n, size = 4.2) +
    annotate("text", x = -0.30, y = 1.20, label = paste0(left_label, "\nn=", stats$set_a_n), size = 2.35, fontface = "bold", lineheight = 0.9) +
    annotate("text", x = 1.42, y = 1.20, label = paste0(right_label, "\nn=", stats$set_b_n), size = 2.35, fontface = "bold", lineheight = 0.9) +
    annotate("label", x = 0.57, y = -1.30,
             label = paste0("Fisher OR=", sprintf("%.2f", stats$fisher_OR), "; p=", p_label),
             size = 2.35, linewidth = 0.2) +
    coord_fixed(xlim = c(-1.2, 2.35), ylim = c(-1.50, 1.45), clip = "off") +
    labs(caption = str_wrap(caption, width = 64)) +
    theme_void(base_size = 8) +
    theme(plot.caption = element_text(size = 6.5, color = "grey35", hjust = 0), plot.margin = margin(4, 4, 4, 4))
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height, dpi = 450, bg = "white")
}

run_human_ora <- function(query, universe) {
  sets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") |>
    transmute(term = gs_name, gene = toupper(gene_symbol)) |>
    distinct()
  enr <- enricher(unique(query), universe = unique(universe), TERM2GENE = sets, pvalueCutoff = 1, qvalueCutoff = 1)
  if (is.null(enr)) return(tibble())
  as_tibble(as.data.frame(enr)) |>
    mutate(pathway = clean_term(Description), is_hsf = hsf_flag(pathway)) |>
    arrange(p.adjust, pvalue)
}

run_mouse_ora <- function(query, universe, tissue) {
  q <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = unique(query), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first") |>
    unname() |>
    na.omit() |>
    unique()
  u <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = unique(universe), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first") |>
    unname() |>
    na.omit() |>
    unique()
  enr <- enrichPathway(q, universe = u, organism = "mouse", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  if (is.null(enr)) return(tibble())
  as_tibble(as.data.frame(enr)) |>
    mutate(tissue = tissue, pathway = clean_term(Description), is_hsf = hsf_flag(pathway), .before = 1) |>
    arrange(p.adjust, pvalue)
}

select_nodes <- function(tbl, max_nodes = 12) {
  if (!nrow(tbl)) return(tbl)
  tbl |>
    filter(p.adjust < 0.05) |>
    mutate(family = term_family(pathway)) |>
    filter(!str_detect(pathway, regex("glucokinase regulatory|viral|infection|HCMV|HIV|VPR|REV|dependant mature mrna|mitotic prophase|glucose metabolism", ignore_case = TRUE))) |>
    group_by(family) |>
    arrange(p.adjust, pvalue, .by_group = TRUE) |>
    slice_head(n = 2) |>
    ungroup() |>
    arrange(p.adjust, pvalue) |>
    mutate(display_label = short_label(pathway)) |>
    distinct(display_label, .keep_all = TRUE) |>
    slice_head(n = max_nodes) |>
    mutate(node_id = row_number())
}

select_named_nodes <- function(tbl, patterns) {
  map_dfr(seq_along(patterns), function(i) {
    tbl |>
      filter(pvalue < 0.05, str_detect(pathway, regex(patterns[[i]], ignore_case = TRUE))) |>
      arrange(p.adjust, pvalue) |>
      slice_head(n = 1) |>
      mutate(display_order = i)
  }) |>
    distinct(ID, .keep_all = TRUE) |>
    arrange(display_order) |>
    mutate(family = term_family(pathway), display_label = short_label(pathway), node_id = row_number())
}

select_bothFDR_nodes <- function(tbl, max_nodes = 10) {
  tbl |>
    filter(pvalue < 0.05, Count >= 2) |>
    mutate(family = term_family(pathway)) |>
    filter(!str_detect(pathway, regex("glucokinase regulatory|viral|infection|HCMV|HIV|VPR|REV|NS1|mitotic|M phase", ignore_case = TRUE))) |>
    group_by(family) |>
    arrange(pvalue, p.adjust, .by_group = TRUE) |>
    slice_head(n = 2) |>
    ungroup() |>
    arrange(pvalue, p.adjust) |>
    mutate(display_label = short_label(pathway)) |>
    distinct(display_label, .keep_all = TRUE) |>
    slice_head(n = max_nodes) |>
    mutate(node_id = row_number())
}

network_edges <- function(nodes) {
  if (nrow(nodes) < 2) return(tibble(from = integer(), to = integer(), shared_n = integer()))
  genes <- strsplit(nodes$geneID, "/", fixed = TRUE)
  map_dfr(seq_len(nrow(nodes) - 1L), function(i) {
    map_dfr(seq.int(i + 1L, nrow(nodes)), function(j) {
      n <- length(intersect(genes[[i]], genes[[j]]))
      if (!n) return(tibble())
      tibble(from = i, to = j, shared_n = n)
    })
  })
}

place_nodes <- function(nodes) {
  family_order <- c("HSF/proteostasis", "ER / protein transport", "RNA processing / nuclear transport", "Cell cycle / chromatin", "Signaling / neuronal", "Metabolism / mitochondria", "Other")
  centers <- tibble(
    family = family_order,
    cx = c(-0.72, -0.25, 0.55, 0.72, 0.20, -0.42, 0.00),
    cy = c(-0.42, 0.48, 0.32, -0.30, -0.70, -0.70, 0.00)
  )
  nodes |>
    mutate(family = factor(family, levels = family_order)) |>
    group_by(family) |>
    arrange(p.adjust, pvalue, .by_group = TRUE) |>
    mutate(i = row_number(), n = n(), angle = if_else(n == 1, 0, 2 * pi * (i - 1) / n), radius = if_else(n == 1, 0, 0.11 + 0.04 * i)) |>
    ungroup() |>
    left_join(centers, by = "family") |>
    mutate(x = cx + radius * cos(angle), y = cy + radius * sin(angle))
}

plot_network <- function(nodes, caption, stem, label_n = 6, width = 7.3, height = 4.9, label_size = 2.9) {
  if (!nrow(nodes)) return(invisible(NULL))
  if (!"force_label" %in% names(nodes)) nodes$force_label <- FALSE
  xy <- place_nodes(nodes) |>
    mutate(label = if_else(is_hsf | force_label | rank(p.adjust, ties.method = "first") <= label_n, short_label(pathway), ""))
  e <- network_edges(xy)
  if (nrow(e)) {
    cut <- max(1L, as.integer(quantile(e$shared_n, 0.55, names = FALSE, type = 1)))
    e <- e |>
      filter(shared_n >= cut) |>
      arrange(desc(shared_n)) |>
      slice_head(n = 20) |>
      left_join(xy |> dplyr::select(node_id, x, y), by = c("from" = "node_id")) |>
      dplyr::rename(x1 = x, y1 = y) |>
      left_join(xy |> dplyr::select(node_id, x, y), by = c("to" = "node_id")) |>
      dplyr::rename(x2 = x, y2 = y)
  } else {
    e <- tibble(x1 = numeric(), y1 = numeric(), x2 = numeric(), y2 = numeric(), shared_n = numeric())
  }
  p <- ggplot() +
    geom_segment(data = e, aes(x1, y1, xend = x2, yend = y2, linewidth = shared_n), color = "grey72", alpha = 0.55) +
    geom_point(data = xy, aes(x, y, size = Count, fill = is_hsf), shape = 21, color = "black", stroke = 0.45) +
    geom_label_repel(data = xy |> filter(label != ""), aes(x, y, label = label), size = label_size, lineheight = 0.9,
                     fill = "white", label.size = 0.18, box.padding = 0.28, point.padding = 0.20,
                     min.segment.length = 0, seed = 41, max.overlaps = Inf) +
    scale_fill_manual(values = c(`TRUE` = "#e66101", `FALSE` = "#b8c7d6"), labels = c(`TRUE` = "HSF/proteostasis", `FALSE` = "Other Reactome")) +
    scale_size_continuous(range = c(3, 9), breaks = scales::breaks_pretty(4)) +
    scale_linewidth_continuous(range = c(0.25, 1.3), guide = "none") +
    coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.05, 0.95), clip = "off") +
    labs(fill = NULL, size = "Overlap genes", caption = str_wrap(caption, width = 88)) +
    theme_void(base_size = 10) +
    theme(legend.position = "right", plot.caption = element_text(size = 7, color = "grey35", hjust = 0, lineheight = 1.05), plot.margin = margin(8, 20, 18, 16))
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(out, paste0(stem, ".png")), p, width = width, height = height, dpi = 450, bg = "white")
  write_csv(xy, file.path(out, paste0(stem, "_nodes.csv")))
  write_csv(e, file.path(out, paste0(stem, "_edges.csv")))
}

sr <- read_csv(file.path(tables, "moller_levet_participant_aware", "MollerLevet_participant_aware_phase_cosinor_gene_results.csv"), show_col_types = FALSE) |>
  mutate(gene = toupper(gene))
fdr_sr <- sr |>
  filter(
    participant_aware_phase_logFC < 0,
    participant_aware_phase_FDR < 0.05,
    delta_mesor < 0,
    mesor_FDR < 0.05
  ) |>
  distinct(gene, .keep_all = TRUE)
fdr_genes <- fdr_sr$gene
write_csv(fdr_sr, file.path(out, "SRdown_bothFDR_full_gene_list.csv"))

## Fig. 3e-f and Supplementary Fig. 3b
fig3 <- read_csv(file.path(tables, "figure3_participant_aware", "Fig3e_participant_aware_gene_membership_full.csv"), show_col_types = FALSE)
fig3_stats <- map_dfr(c("Excitatory neurons", "Inhibitory neurons"), function(group) {
  dat <- fig3 |> filter(cell_group == group)
  universe <- unique(dat$gene_upper)
  sr_set <- intersect(fdr_genes, universe)
  aging_set <- dat |> filter(aging_slope < 0, aging_p < 0.05) |> pull(gene_upper) |> unique()
  fisher_summary(sr_set, aging_set, universe) |> mutate(group = group, .before = 1)
})
write_csv(fig3_stats, file.path(out, "F3e_bothFDR_fisher_summary.csv"))
fig3_membership <- fig3 |>
  transmute(
    group = cell_group,
    gene = gene_upper,
    aging_slope,
    aging_p,
    sr_down = gene_upper %in% fdr_genes,
    neuronal_aging_down = aging_slope < 0 & aging_p < 0.05,
    convergence = sr_down & neuronal_aging_down
  ) |>
  distinct(group, gene, .keep_all = TRUE)
write_csv(fig3_membership, file.path(out, "F3e_SF3b_bothFDR_gene_membership.csv"))

venn_left <- fig3_stats |> filter(group == "Excitatory neurons")
venn_right <- fig3_stats |> filter(group == "Inhibitory neurons")
venn_df <- bind_rows(
  tibble(group = "Excitatory neurons", x_left = -1.15, x_right = -0.25, stats = list(venn_left)),
  tibble(group = "Inhibitory neurons", x_left = 2.15, x_right = 3.05, stats = list(venn_right))
)
circle_df <- map_dfr(seq_len(nrow(venn_df)), function(i) {
  theta <- seq(0, 2 * pi, length.out = 240)
  bind_rows(
    tibble(group = venn_df$group[i], circle = "SR", x = venn_df$x_left[i] + 0.78 * cos(theta), y = 0.78 * sin(theta)),
    tibble(group = venn_df$group[i], circle = "Aging", x = venn_df$x_right[i] + 0.78 * cos(theta), y = 0.78 * sin(theta))
  )
})
p3e <- ggplot() +
  geom_polygon(data = circle_df, aes(x, y, group = interaction(group, circle), fill = circle, color = circle), alpha = 0.22, linewidth = 0.7) +
  annotate("text", x = c(-1.39, -0.70, -0.01, 1.91, 2.60, 3.29), y = 0,
           label = c(venn_left$set_a_only_n, venn_left$overlap_n, venn_left$set_b_only_n, venn_right$set_a_only_n, venn_right$overlap_n, venn_right$set_b_only_n), size = 3.7) +
  annotate("text", x = c(-1.33, -0.07, 1.97, 3.17), y = 1.12,
           label = c(paste0("SR-down\nhuman blood\nn=", venn_left$set_a_n), paste0("Excitatory\naging-down\nslope<0; p<0.05\nn=", venn_left$set_b_n), paste0("SR-down\nhuman blood\nn=", venn_right$set_a_n), paste0("Inhibitory\naging-down\nslope<0; p<0.05\nn=", venn_right$set_b_n)), size = 2.35, lineheight = 0.9) +
  annotate("text", x = c(-0.70, 2.60), y = 1.68, label = c("Excitatory neurons", "Inhibitory neurons"), fontface = "bold", size = 3.2) +
  annotate("label", x = c(-0.70, 2.60), y = -0.94,
           label = c(paste0("Fisher OR=", sprintf("%.2f", venn_left$fisher_OR), "; p=", format(venn_left$fisher_p, digits = 2, scientific = TRUE)), paste0("Fisher OR=", sprintf("%.2f", venn_right$fisher_OR), "; p=", format(venn_right$fisher_p, digits = 2, scientific = TRUE))), size = 2.2, linewidth = 0.18) +
  scale_fill_manual(values = c(SR = "#dceef9", Aging = "#e0f2df"), guide = "none") +
  scale_color_manual(values = c(SR = "#2b83c6", Aging = "#2ca25f"), guide = "none") +
  coord_fixed(xlim = c(-2.25, 4.15), ylim = c(-1.18, 1.90), clip = "off") +
  labs(caption = threshold_text) +
  theme_void(base_size = 8) +
  theme(plot.caption = element_text(size = 6.5, color = "grey35", hjust = 0), plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(out, "F3e_bothFDR_SRdown_neuronal_agingdown_venn.pdf"), p3e, width = 6.2, height = 3.0, bg = "white")
ggsave(file.path(out, "F3e_bothFDR_SRdown_neuronal_agingdown_venn.png"), p3e, width = 6.2, height = 3.0, dpi = 450, bg = "white")

for (group in c("Excitatory neurons", "Inhibitory neurons")) {
  dat <- fig3 |> filter(cell_group == group)
  universe <- unique(dat$gene_upper)
  query <- dat |> filter(gene_upper %in% fdr_genes, aging_slope < 0, aging_p < 0.05) |> pull(gene_upper) |> unique()
  ora <- run_human_ora(query, universe)
  stem <- if (group == "Excitatory neurons") "F3f_bothFDR_clean_Reactome_network" else "SF3b_bothFDR_clean_Reactome_network"
  write_csv(ora, file.path(out, paste0(stem, "_full_ORA.csv")))
  display_note <- "representative Reactome terms with nominal p < 0.05 and at least two overlapping genes; BH FDR is retained in source data"
  nodes <- select_bothFDR_nodes(ora) |>
    filter(!is_hsf | p.adjust < 0.05) |>
    mutate(
      family = term_family(pathway),
      display_label = short_label(pathway),
      node_id = row_number(),
      force_label = FALSE,
      display_status = "nominal Reactome p < 0.05"
    )
  plot_network(
    nodes,
    paste0(threshold_text, "; neuronal aging slope < 0 and p < 0.05; ", display_note),
    stem,
    label_n = if_else(group == "Excitatory neurons", nrow(nodes), 6),
    width = if_else(group == "Excitatory neurons", 8.3, 7.3),
    height = if_else(group == "Excitatory neurons", 5.3, 4.9),
    label_size = if_else(group == "Excitatory neurons", 2.7, 2.9)
  )
}

## Supplementary Fig. 4e-g
jen <- read_csv(file.path(paths$source_data, "GSE134080_or_JenAge", "JenAge_01_healthy_blood_age_log2RPKM_results.csv"), show_col_types = FALSE) |>
  transmute(gene = toupper(gene), age_logFC = jenage_blood_age_logFC_per_decade, age_p = jenage_blood_age_p, age_FDR = jenage_blood_age_fdr) |>
  distinct(gene, .keep_all = TRUE)
sf4_universe <- intersect(sr$gene, jen$gene)
sf4_sr <- intersect(fdr_genes, sf4_universe)
sf4_age <- jen |> filter(gene %in% sf4_universe, age_logFC < 0) |> pull(gene)
sf4_stats <- fisher_summary(sf4_sr, sf4_age, sf4_universe)
write_csv(sf4_stats |> mutate(sr_threshold = threshold_text, aging_threshold = "JenAge blood age logFC per decade < 0"), file.path(out, "SF4e_bothFDR_fisher_summary.csv"))
sf4_overlap <- intersect(sf4_sr, sf4_age)
write_csv(jen |> filter(gene %in% sf4_overlap) |> left_join(fdr_sr, by = "gene"), file.path(out, "SF4e_bothFDR_overlap_genes.csv"))
plot_venn(sf4_stats, "SR-down\nboth FDR<0.05", "Blood-aging down\nlogFC<0", threshold_text, "SF4e_bothFDR_SR_aging_overlap")

sf4_ora <- run_human_ora(sf4_overlap, sf4_universe)
write_csv(sf4_ora, file.path(out, "SF4f_SF4g_bothFDR_full_Reactome_ORA.csv"))
plot_network(
  select_bothFDR_nodes(sf4_ora, max_nodes = 16),
  paste0(threshold_text, "; JenAge blood-aging logFC < 0; plotted terms have nominal Reactome p < 0.05 and at least two overlapping genes"),
  "SF4f_bothFDR_clean_Reactome_network"
)

assign_family <- function(x) {
  case_when(
    hsf_flag(x) ~ "HSF1-associated heat-shock response",
    str_detect(x, regex("SUMO", ignore_case = TRUE)) ~ "SUMOylation programs",
    str_detect(x, regex("RNA|SPLIC|SNRNA|SNRNP|RIBONUCLEOPROTEIN|TRNA|RRNA", ignore_case = TRUE)) ~ "RNA processing / splicing",
    str_detect(x, regex("NUCLEAR.*TRANSPORT|RNA EXPORT|EXPORTIN|IMPORTIN", ignore_case = TRUE)) ~ "Nuclear transport / RNA export",
    str_detect(x, regex("CELL CYCLE|MITOTIC|CHROMOSOME|DNA REPLICATION", ignore_case = TRUE)) ~ "Cell-cycle regulation",
    str_detect(x, regex("DNA DAMAGE|TP53|P53", ignore_case = TRUE)) ~ "DNA damage / TP53 regulation",
    TRUE ~ NA_character_
  )
}
sf4_family <- sf4_ora |>
  filter(p.adjust < 0.05) |>
  mutate(family = assign_family(pathway)) |>
  filter(!is.na(family)) |>
  group_by(family) |>
  summarise(
    n_terms = n(),
    overlap_genes = paste(unique(unlist(strsplit(geneID, "/", fixed = TRUE))), collapse = ";"),
    overlap_n = length(unique(unlist(strsplit(geneID, "/", fixed = TRUE)))),
    family_FDR = min(p.adjust),
    .groups = "drop"
  ) |>
  arrange(family_FDR)
write_csv(sf4_family, file.path(out, "SF4g_bothFDR_Reactome_family_summary.csv"))
if (nrow(sf4_family)) {
  p4g <- sf4_family |>
    mutate(family = factor(family, levels = rev(family)), label = paste0(overlap_n, " genes; ", n_terms, " terms"), hsf = family == "HSF1-associated heat-shock response") |>
    ggplot(aes(-log10(family_FDR), family, fill = hsf)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label), hjust = -0.04, size = 2.6) +
    scale_fill_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "#7297c7"), guide = "none") +
    coord_cartesian(xlim = c(0, max(-log10(sf4_family$family_FDR)) * 1.45), clip = "off") +
    labs(x = "-log10(FDR)", y = NULL, caption = str_wrap(paste0(threshold_text, "; Reactome pathway families require BH FDR < 0.05"), width = 100)) +
    theme_classic(base_size = 9) +
    theme(plot.caption = element_text(size = 6.5, color = "grey35", hjust = 0), plot.margin = margin(5, 70, 5, 5))
  ggsave(file.path(out, "SF4g_bothFDR_Reactome_family_barplot.pdf"), p4g, width = 5.2, height = 2.9, bg = "white")
  ggsave(file.path(out, "SF4g_bothFDR_Reactome_family_barplot.png"), p4g, width = 5.2, height = 2.9, dpi = 450, bg = "white")
}

## Fig. 4c/e/f/g and Supplementary Fig. 5a-f
fig4_dir <- file.path(paths$source_data, "Takasugi_2024")
traj <- read_csv(file.path(fig4_dir, "Figure_4_trajectory_features.csv"), show_col_types = FALSE) |>
  mutate(tissue = as.character(tissue), gene = toupper(feature_upper), feature_mouse = feature)
tissues <- c("Aorta", "Brain", "Heart", "Kidney", "Liver", "Lung", "Muscle", "Skin")
multiomics_membership_path <- file.path(tables, "multiomics_aging_definition", "F4_SF5_participant_aware_tissue_membership.csv")
multiomics_summary_path <- file.path(tables, "multiomics_aging_definition", "F4c_SF5_participant_aware_overlap_summary.csv")
stopifnot(file.exists(multiomics_membership_path), file.exists(multiomics_summary_path))

gene_map <- traj |>
  distinct(tissue, gene, feature_mouse)
multiomics_membership <- read_csv(multiomics_membership_path, show_col_types = FALSE) |>
  dplyr::select(-dplyr::any_of(c("participant_aware_SR_down", "convergence"))) |>
  mutate(
    participant_aware_SR_down = gene %in% fdr_genes,
    convergence = aging_decline & participant_aware_SR_down
  ) |>
  left_join(gene_map, by = c("tissue", "gene"))
multiomics_summary <- multiomics_membership |>
  group_by(tissue) |>
  group_modify(~{
    universe <- unique(.x$gene)
    stats <- fisher_summary(
      .x$gene[.x$participant_aware_SR_down],
      .x$gene[.x$aging_decline],
      universe
    )
    transmute(
      stats,
      universe_n,
      sr_n = set_a_n,
      aging_n = set_b_n,
      overlap_n,
      sr_only_n = set_a_only_n,
      aging_only_n = set_b_only_n,
      fisher_or = fisher_OR,
      fisher_p
    )
  }) |>
  ungroup() |>
  mutate(
    sr_threshold = threshold_text,
    aging_threshold = "replicate-level log2(value+1) RNA age slope < 0 and at least one WTL/LSF protein age slope < 0",
    background_universe = "genes tested in the original GSE39445 phase-adjusted analysis and measured in tissue RNA plus protein"
  )
decline <- multiomics_membership |>
  filter(convergence) |>
  distinct(tissue, gene, feature_mouse, rna_slope, protein_layers_down)
write_csv(decline, file.path(out, "F4_SF5_bothFDR_tissue_convergence_genes.csv"))
write_csv(multiomics_membership, file.path(out, "F4_SF5_bothFDR_tissue_membership_full.csv"))

brain_row <- multiomics_summary |> filter(tissue == "Brain")
stopifnot(nrow(brain_row) == 1L, brain_row$aging_n == 1634L)
brain_stats <- tibble(
  set_a_n = brain_row$sr_n,
  set_b_n = brain_row$aging_n,
  overlap_n = brain_row$overlap_n,
  set_a_only_n = brain_row$sr_only_n,
  set_b_only_n = brain_row$aging_only_n,
  universe_n = brain_row$universe_n,
  fisher_OR = brain_row$fisher_or,
  fisher_p = brain_row$fisher_p
)
write_csv(brain_row, file.path(out, "F4c_bothFDR_fisher_summary.csv"))
plot_venn(
  brain_stats,
  "SR-down\nhuman blood",
  "RNA+protein\naging decline",
  "Participant-aware SR: phase-adjusted logFC < 0 and FDR < 0.05; delta MESOR < 0 and MESOR FDR < 0.05; aging: RNA slope < 0 plus WTL/LSF slope < 0",
  "F4c_bothFDR_SR_multiomics_overlap",
  width = 3.4,
  height = 2.7
)

reactome <- map_dfr(tissues, function(tissue) {
  query <- decline |> filter(.data$tissue == .env$tissue) |> pull(feature_mouse) |> unique()
  universe <- multiomics_membership |> filter(.data$tissue == .env$tissue) |> pull(feature_mouse) |> unique()
  run_mouse_ora(query, universe, tissue)
})
write_csv(reactome, file.path(out, "F4_SF5_bothFDR_full_Reactome_ORA.csv"))
write_csv(
  reactome |>
    group_by(tissue) |>
    summarise(
      fdr_significant_terms_n = sum(p.adjust < 0.05, na.rm = TRUE),
      fdr_significant_hsf_terms_n = sum(p.adjust < 0.05 & is_hsf, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(out, "F4_SF5_bothFDR_Reactome_network_significance_audit.csv")
)

tissue_context_patterns <- list(
  Brain = c("calnexin|calreticulin|n glycan|endoplasmic|er to golgi|protein transport", "neuronal|ntrk|ngf|synap"),
  Liver = c("rna processing|splicing|snrnp|trna|rrna", "nuclear.*transport|rna export|sumoylation"),
  Lung = c("rho gtpase|foxo|developmental", "transcription|rna processing"),
  Aorta = c("nuclear.*transport|rna export|sumoylation", "chromatin|histone"),
  Heart = c("respiratory|mitochondrial|complex ii|atp", "rrna|rna polymerase"),
  Muscle = c("nuclear envelope|signaling|met", "respiratory|pyruvate"),
  Kidney = c("complex ii|pyruvate|respiratory", "ras signaling|nitric oxide|nmdar"),
  Skin = c("golgi|er |cargo|protein transport", "antigen presentation|interleukin")
)

select_tissue_nodes <- function(ora, tissue, max_nodes = 16) {
  select_bothFDR_nodes(ora, max_nodes = max_nodes) |>
    mutate(
      force_label = FALSE,
      display_status = "nominal Reactome p < 0.05; overlap genes >= 2"
    )
}

for (tissue in tissues) {
  ora <- reactome |> filter(.data$tissue == .env$tissue)
  nodes <- select_tissue_nodes(ora, tissue)
  reactome_display_cutoff <- "plotted terms have nominal Reactome p < 0.05 and at least two overlapping genes; no terms are forced"
  stem <- case_when(
    tissue == "Brain" ~ "F4e_bothFDR_clean_Reactome_network",
    tissue == "Liver" ~ "F4f_bothFDR_clean_Reactome_network",
    TRUE ~ paste0("SF5_", tolower(tissue), "_bothFDR_clean_Reactome_network")
  )
  plot_network(
    nodes,
    paste0(
      threshold_text,
      "; aging decline = RNA slope < 0 plus at least one WTL/LSF protein slope < 0; ",
      reactome_display_cutoff
    ),
    stem,
    label_n = if_else(tissue %in% c("Brain", "Liver"), nrow(nodes), 4L),
    width = if_else(tissue %in% c("Brain", "Liver"), 8.3, 7.3),
    height = if_else(tissue %in% c("Brain", "Liver"), 5.3, 4.9),
    label_size = if_else(tissue %in% c("Brain", "Liver"), 2.7, 2.9)
  )
}

focus <- tibble(
  pathway = c(
    "CELLULAR RESPONSE TO HEAT STRESS",
    "HSF1 DEPENDENT TRANSACTIVATION",
    "HSF1-ASSOCIATED HEAT-SHOCK RESPONSE"
  ),
  pattern = c(
    "^Cellular response to heat stress$",
    "^Hsf1-dependent transactivation$",
    "^Regulation of hsf1-mediated heat shock response$"
  )
)
heat <- expand_grid(tissue = tissues, pathway = focus$pathway) |>
  mutate(stats = map2(tissue, pathway, function(tissue, pathway) {
    pat <- focus$pattern[match(pathway, focus$pathway)]
    hit <- reactome |>
      filter(.data$tissue == .env$tissue, str_detect(.data$pathway, regex(pat, ignore_case = TRUE))) |>
      arrange(pvalue, p.adjust) |>
      slice_head(n = 1)
    if (!nrow(hit)) return(tibble(pvalue = 1, FDR = 1, overlap_genes_n = 0L))
    hit |> transmute(pvalue, FDR = p.adjust, overlap_genes_n = Count)
  })) |>
  unnest(stats) |>
  mutate(logP = -log10(pvalue), logFDR = -log10(FDR))
write_csv(heat, file.path(out, "F4g_bothFDR_heatmap_source.csv"))
pheat <- heat |>
  mutate(tissue = factor(tissue, levels = tissues), pathway = factor(pathway, levels = rev(focus$pathway)), value = pmin(logP, 5)) |>
  ggplot(aes(tissue, pathway, fill = value)) +
  geom_tile(color = "#BDBDBD", linewidth = 0.35) +
  scale_fill_gradient(low = "white", high = "#3B7DBA", limits = c(0, 5), breaks = 0:5, name = expression(-log[10](italic(P)))) +
  labs(x = NULL, y = NULL, caption = str_wrap(paste0(threshold_text, "; fill = -log10 nominal Reactome p; raw p, BH FDR, and overlap counts are reported in source data"), width = 74)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.text.y = element_text(face = "bold", size = 7), axis.line = element_blank(), axis.ticks = element_blank(), plot.caption = element_text(size = 6.5, color = "grey35", hjust = 0, lineheight = 1.05), plot.margin = margin(5, 12, 12, 5))
ggsave(file.path(out, "F4g_bothFDR_recurrent_Reactome_heatmap.pdf"), pheat, width = 5.6, height = 3.0, bg = "white")
ggsave(file.path(out, "F4g_bothFDR_recurrent_Reactome_heatmap.png"), pheat, width = 5.6, height = 3.0, dpi = 450, bg = "white")

audit <- bind_rows(
  fig3_stats |> transmute(panel = paste0("F3e_", group), set_a_n, set_b_n, overlap_n, fisher_OR, fisher_p),
  sf4_stats |> transmute(panel = "SF4e", set_a_n, set_b_n, overlap_n, fisher_OR, fisher_p),
  brain_stats |> transmute(panel = "F4c", set_a_n, set_b_n, overlap_n, fisher_OR, fisher_p)
) |>
  mutate(sr_down_total_n = length(fdr_genes), sr_threshold = threshold_text)
write_csv(audit, file.path(out, "bothFDR_panel_audit_summary.csv"))

nominal_genes_for_comparison <- sr |>
  filter(
    participant_aware_phase_logFC < 0,
    participant_aware_phase_p < 0.05,
    delta_mesor < 0
  ) |>
  pull(gene) |>
  unique()
write_csv(
  tibble(
    comparison = "participant-aware human blood SR-down",
    nominal_p_definition_n = length(nominal_genes_for_comparison),
    both_FDR_definition_n = length(fdr_genes),
    shared_n = length(intersect(nominal_genes_for_comparison, fdr_genes)),
    nominal_only_n = length(setdiff(nominal_genes_for_comparison, fdr_genes)),
    both_FDR_only_n = length(setdiff(fdr_genes, nominal_genes_for_comparison))
  ),
  file.path(out, "SRdown_nominal_vs_bothFDR_gene_set_comparison.csv")
)

nominal_audit_path <- file.path(paths$tables, "sensitivity", "nominal_panel_audit_summary.csv")
if (file.exists(nominal_audit_path)) {
  nominal_audit <- read_csv(nominal_audit_path, show_col_types = FALSE) |>
    mutate(threshold_version = "nominal phase p < 0.05 plus negative delta MESOR", .before = 1)
  write_csv(
    bind_rows(
      nominal_audit,
      audit |> mutate(threshold_version = "phase FDR < 0.05 plus negative delta MESOR with MESOR FDR < 0.05", .before = 1)
    ),
    file.path(out, "panel_statistics_nominal_vs_bothFDR_comparison.csv")
  )
}

writeLines(c(
  "Human blood SR-down both-FDR panel review set",
  "",
  "SR-down definition:",
  "participant-aware phase-adjusted logFC < 0; phase-model FDR < 0.05; delta MESOR < 0; MESOR FDR < 0.05.",
  "The participant-aware model uses limma duplicateCorrelation with participant identity as the repeated-measures block.",
  "Fisher overlap tests are one-sided enrichment tests (alternative = greater) within each panel-specific shared tested-gene universe.",
  "Reactome pathway P values are adjusted by the Benjamini-Hochberg procedure where FDR is reported.",
  "HSF1/proteostasis terms are never inserted contextually; they appear only when they meet the stated panel-specific enrichment rule.",
  "All Reactome network panels plot terms with nominal Reactome p<0.05 and at least two overlapping genes; global Reactome BH-FDR values remain reported in source data.",
  "",
  "Saved plot groups:",
  "Fig. 3e-f; Fig. 4c,e-g; Supplementary Fig. 3b; Supplementary Fig. 4e-g; Supplementary Fig. 5a-f.",
  "PDF and high-resolution PNG versions are provided with panel-specific CSV source tables.",
  "",
  "Important Fig. 4c result:",
  paste0(
    "Within the ", brain_row$universe_n, "-gene brain RNA+protein universe, SR-down n=", brain_row$sr_n,
    ", aging-decline n=", brain_row$aging_n, ", overlap n=", brain_row$overlap_n,
    ", one-sided Fisher OR=", sprintf("%.3f", brain_row$fisher_or),
    ", p=", format(brain_row$fisher_p, digits = 4), "."
  ),
  ifelse(brain_row$fisher_p < 0.05 & brain_row$fisher_or > 1, "This comparison shows enrichment.", "This comparison does not show enrichment.")
), file.path(out, "README.txt"))

message("Both-FDR SR-down panel review set written to: ", out)
