suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(msigdbr)
  library(limma)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(tibble)
})

fig1_dir <- file.path(paths$github_ready, "Figure1")
fig1_tables <- file.path(fig1_dir, "Tables")
stage_dir <- file.path(paths$tables, "participant_aware_staging")
stage_fig_dir <- file.path(paths$figures, "participant_aware_staging")
final_fig_dir <- file.path(paths$figures, "final")
for (d in c(fig1_tables, file.path(fig1_dir, "Plots"), stage_dir, stage_fig_dir, final_fig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

gene_results <- read_csv(
  file.path(paths$tables, "moller_levet_participant_aware", "MollerLevet_participant_aware_phase_cosinor_gene_results.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    gene_symbol = gene,
    log_fc = participant_aware_phase_logFC,
    t = participant_aware_phase_t,
    p_value = participant_aware_phase_p,
    adj_p_val = participant_aware_phase_FDR,
    delta_mesor = delta_mesor,
    mesor_p = mesor_p,
    mesor_FDR = mesor_FDR,
    rhythmic_SR_effect_amplitude = rhythmic_SR_effect_amplitude,
    SR_effect_acrophase_h = SR_effect_acrophase_h
  ) %>%
  filter(!is.na(gene_symbol), !is.na(t)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

write_csv(gene_results, file.path(stage_dir, "F1_gene_level_participant_aware_phase_cosinor_results.csv"))
write_csv(gene_results, file.path(fig1_tables, "Figure_1E_1F_gene_level_phase_adjusted_limma.csv"))
write_csv(gene_results, file.path(fig1_dir, "Figure_1_SR_gene_level_results.csv"))

msig <- suppressWarnings(msigdbr(species = "Homo sapiens"))
get_msig <- function(collection, subcollection = NULL) {
  if ("gs_collection" %in% names(msig)) {
    out <- msig %>% filter(gs_collection == collection)
    if (!is.null(subcollection)) out <- out %>% filter(gs_subcollection == subcollection)
  } else {
    out <- msig %>% filter(category == collection)
    if (!is.null(subcollection)) out <- out %>% filter(subcategory == subcollection)
  }
  out %>% select(gs_name, gene_symbol) %>% distinct() %>% filter(gene_symbol %in% gene_results$gene_symbol)
}

run_camera <- function(gene_sets) {
  stat <- gene_results$t
  names(stat) <- gene_results$gene_symbol
  idx <- split(match(gene_sets$gene_symbol, names(stat)), gene_sets$gs_name)
  idx <- idx[vapply(idx, function(i) length(i) >= 10 && length(i) <= 500, logical(1))]
  cameraPR(stat, idx, use.ranks = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("GeneSet") %>%
    transmute(
      NGenes = NGenes,
      Direction = Direction,
      PValue = PValue,
      FDR = FDR,
      GeneSet = GeneSet,
      logP = -log10(pmax(PValue, .Machine$double.xmin)),
      GeneSet_clean = str_replace_all(str_remove(GeneSet, "^REACTOME_|^HALLMARK_"), "_", " ")
    ) %>%
    arrange(PValue)
}

hallmark_all <- run_camera(get_msig("H"))
write_csv(hallmark_all, file.path(fig1_tables, "Figure_1F_CAMERA_Hallmark_all.csv"))
write_csv(hallmark_all, file.path(fig1_tables, "Figure_1E_CAMERA_Hallmark_all.csv"))
write_csv(hallmark_all, file.path(fig1_dir, "Hallmark_CAMERA_all.csv"))
write_csv(hallmark_all, file.path(fig1_dir, "Figure_1E_CAMERA_hallmark_full.csv"))
write_csv(hallmark_all, file.path(stage_dir, "F1f_participant_aware_Hallmark_CAMERA_all.csv"))

hallmark_plot <- hallmark_all %>%
  slice_head(n = 7) %>%
  mutate(label = paste0(GeneSet_clean, " (n=", NGenes, ")"), label = fct_reorder(label, logP))
write_csv(hallmark_plot, file.path(fig1_tables, "Figure_1F_source_data.csv"))
write_csv(hallmark_plot, file.path(fig1_dir, "Figure_1F_source_data.csv"))
write_csv(hallmark_plot, file.path(stage_dir, "F1f_participant_aware_Hallmark_plot_source_data.csv"))

p_hallmark <- ggplot(hallmark_plot, aes(x = logP, y = label, fill = Direction)) +
  geom_col(width = 0.76, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(Up = "#c9252d", Down = "#3f7fb3"), breaks = c("Up", "Down")) +
  labs(title = "Top Hallmark pathways", x = expression(-log[10]~"(p-value)"), y = NULL, fill = "Direction") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )
for (pfile in c(file.path(fig1_dir, "Plots", "Figure_1F.png"), file.path(fig1_dir, "Figure_1F.png"), file.path(stage_fig_dir, "F1f_participant_aware_Hallmark_barplot.png"), file.path(final_fig_dir, "F1f_participant_aware_Hallmark_barplot.png"))) {
  ggsave(pfile, p_hallmark, width = 5.2, height = 3.6, dpi = 300, bg = "white")
}
for (pfile in c(file.path(stage_fig_dir, "F1f_participant_aware_Hallmark_barplot.pdf"), file.path(final_fig_dir, "F1f_participant_aware_Hallmark_barplot.pdf"))) {
  ggsave(pfile, p_hallmark, width = 5.2, height = 3.6, bg = "white")
}

reactome_sets <- get_msig("C2", "CP:REACTOME")
reactome_all <- run_camera(reactome_sets)
write_csv(reactome_all, file.path(fig1_tables, "Figure_1F_CAMERA_Reactome_all.csv"))
write_csv(reactome_all %>% transmute(pathway = GeneSet, pathway_clean = GeneSet_clean, n_genes = NGenes, direction = Direction, p_value = PValue, fdr = FDR, logp = logP), file.path(fig1_tables, "Figure_1F_CAMERA_reactome_full.csv"))
write_csv(reactome_all, file.path(fig1_dir, "Reactome_CAMERA_all.csv"))
write_csv(reactome_all, file.path(fig1_dir, "Figure_1F_Reactome_ORA_full.csv"))
write_csv(reactome_all, file.path(stage_dir, "F1g_participant_aware_Reactome_CAMERA_all.csv"))

term_genes <- split(reactome_sets$gene_symbol, reactome_sets$gs_name)
sig <- reactome_all %>% filter(FDR < 0.05) %>% arrange(FDR, PValue)
kept <- character()
nonred <- logical(nrow(sig))
if (nrow(sig) > 0) {
  for (i in seq_len(nrow(sig))) {
    gs <- sig$GeneSet[i]
    g <- unique(term_genes[[gs]])
    redundant <- FALSE
    if (length(kept)) {
      for (kgs in kept) {
        kg <- unique(term_genes[[kgs]])
        if (length(intersect(g, kg)) / min(length(g), length(kg)) >= 0.75) {
          redundant <- TRUE
          break
        }
      }
    }
    nonred[i] <- !redundant
    if (!redundant) kept <- c(kept, gs)
  }
}
reactome_nonredundant <- sig %>% mutate(nonredundant = nonred)
write_csv(reactome_nonredundant, file.path(fig1_tables, "Figure_1F_CAMERA_reactome_nonredundant.csv"))
write_csv(reactome_nonredundant, file.path(fig1_dir, "Reactome_CAMERA_nonredundant_flagged.csv"))
write_csv(reactome_nonredundant %>% filter(nonredundant), file.path(fig1_dir, "Reactome_significant_nonredundant_table.csv"))

network_terms <- reactome_nonredundant %>% filter(nonredundant) %>% arrange(FDR, PValue) %>% slice_head(n = 28)
if (nrow(network_terms) < 8) network_terms <- sig %>% arrange(FDR, PValue) %>% slice_head(n = 28) %>% mutate(nonredundant = TRUE)
hsf1_term <- reactome_all %>%
  filter(GeneSet == "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE") %>%
  mutate(nonredundant = TRUE)
if (nrow(hsf1_term) == 1 && !hsf1_term$GeneSet %in% network_terms$GeneSet) {
  network_terms <- bind_rows(network_terms, hsf1_term) %>% arrange(FDR, PValue)
}
node_df <- network_terms %>%
  mutate(
    node_id = as.character(row_number()),
    label = if_else(GeneSet == "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE", "HSF1 HEAT SHOCK RESPONSE", GeneSet_clean),
    degree = 0L,
    community = NA_integer_,
    class_proteostasis = str_detect(str_to_upper(GeneSet), "HSF|HEAT|CHAPERONE|UNFOLDED|PROTEIN_FOLD|PROTEOSTAS|HSP|SUMO"),
    rank_global = row_number(),
    show_label = rank_global <= 12 | class_proteostasis | GeneSet == "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",
    label_plot = if_else(show_label, label, "")
  ) %>%
  mutate(
    label_plot = if_else(
      label %in% c(
        "DEFECTIVE PYROPTOSIS",
        "TRANSPORT OF MATURE TRANSCRIPT TO CYTOPLASM",
        "MITOTIC PROPHASE",
        "RRNA MODIFICATION IN THE NUCLEUS AND CYTOSOL"
      ),
      "",
      label_plot
    )
  )
edges <- list()
if (nrow(node_df) >= 2) {
  for (i in seq_len(nrow(node_df) - 1)) {
    for (j in (i + 1):nrow(node_df)) {
      g1 <- unique(term_genes[[node_df$GeneSet[i]]])
      g2 <- unique(term_genes[[node_df$GeneSet[j]]])
      ov <- length(intersect(g1, g2)) / min(length(g1), length(g2))
      if (!is.na(ov) && ov >= 0.25) edges[[length(edges) + 1]] <- data.frame(from = node_df$node_id[i], to = node_df$node_id[j], overlap_coef = ov)
    }
  }
}
edge_df <- if (length(edges)) bind_rows(edges) else data.frame(from = character(), to = character(), overlap_coef = numeric())
if (nrow(edge_df) > 0) {
  graph_tmp <- graph_from_data_frame(edge_df, vertices = node_df %>% select(name = node_id, everything()), directed = FALSE)
  node_df$degree <- as.integer(degree(graph_tmp)[node_df$node_id])
  node_df$community <- as.integer(membership(cluster_louvain(graph_tmp))[node_df$node_id])
}
node_source <- node_df %>% select(GeneSet, label, logP, Direction, NGenes, degree, community, class_proteostasis, rank_global, show_label, label_plot, PValue, FDR)
for (pfile in c(file.path(fig1_tables, "Figure_1F_nodes_source_data.csv"), file.path(fig1_dir, "Figure_1F_nodes_source_data.csv"), file.path(fig1_dir, "Figure_1F_nodes.csv"), file.path(stage_dir, "F1g_participant_aware_network_nodes.csv"))) write_csv(node_source, pfile)
for (pfile in c(file.path(fig1_tables, "Figure_1F_edges_source_data.csv"), file.path(fig1_dir, "Figure_1F_edges_source_data.csv"), file.path(fig1_dir, "Figure_1F_edges.csv"), file.path(stage_dir, "F1g_participant_aware_network_edges.csv"))) write_csv(edge_df, pfile)

vertices <- node_df %>% select(name = node_id, label, label_plot, Direction, NGenes)
graph <- if (nrow(edge_df) > 0) graph_from_data_frame(edge_df, vertices = vertices, directed = FALSE) else make_empty_graph(n = nrow(vertices), directed = FALSE) %>% set_vertex_attr("name", value = vertices$name)
p_reactome <- ggraph(graph, layout = "fr")
if (ecount(graph) > 0) p_reactome <- p_reactome + geom_edge_link(aes(width = overlap_coef), color = "grey72", alpha = 0.55, show.legend = FALSE)
p_reactome <- p_reactome +
  geom_node_point(aes(size = NGenes, fill = Direction), shape = 21, color = "black", stroke = 0.35, alpha = 0.95) +
  geom_node_label(aes(label = label_plot), repel = TRUE, size = 4.5, label.size = 0.22, fill = "white", label.padding = unit(0.12, "lines"), na.rm = TRUE) +
  scale_fill_manual(values = c(Up = "#c9252d", Down = "#3f7fb3"), breaks = c("Up", "Down")) +
  scale_size_continuous(range = c(3, 12), name = "N genes") +
  labs(title = "Reactome pathway overlap network", fill = "Direction") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(size = 18, hjust = 0, margin = margin(b = 8)), legend.position = "right", plot.margin = margin(5, 5, 5, 5))
for (pfile in c(file.path(stage_fig_dir, "F1g_participant_aware_Reactome_network.png"), file.path(final_fig_dir, "F1g_participant_aware_Reactome_network.png"))) ggsave(pfile, p_reactome, width = 7, height = 7, dpi = 300, bg = "white")
for (pfile in c(file.path(stage_fig_dir, "F1g_participant_aware_Reactome_network.pdf"), file.path(final_fig_dir, "F1g_participant_aware_Reactome_network.pdf"))) ggsave(pfile, p_reactome, width = 7, height = 7, bg = "white")
