# Full figure-to-script/source manifest for final package consistency.
# This is a lightweight audit: it does not rewrite panels, but confirms that the
# final R runner documents source-data/script coverage for every figure in the
# resubmission deck and that the expected source folders/files are present.

count_files <- function(rel_dir, pattern) {
  full <- file.path(paths$source_materials, rel_dir)
  if (!dir.exists(full)) return(0L)
  length(list.files(full, pattern = pattern, recursive = TRUE, full.names = TRUE))
}

required_present <- function(rel_paths) {
  if (length(rel_paths) == 0) return(TRUE)
  all(file.exists(file.path(paths$source_materials, rel_paths)))
}

manifest <- data.frame(
  figure = c(
    "Fig. 1", "Fig. 2", "Fig. 3", "Fig. 4", "Fig. 5", "Fig. 6",
    "Supplementary Fig. 1", "Supplementary Fig. 2", "Supplementary Fig. 3",
    "Supplementary Fig. 4", "Supplementary Fig. 5", "Supplementary Fig. 6"
  ),
  title = c(
    "HSF1-associated proteostasis pathways are attenuated under chronic sleep restriction",
    "Human prefrontal cortex transcriptome exhibits directionally organized aging trajectories",
    "HSF1-associated proteostasis decline emerges as a recurrent molecular axis linking chronic sleep insufficiency and neuronal aging",
    "Multi-omics convergence of aging and sleep restriction reveals systemic HSF1-associated proteostasis decline",
    "Spatial transcriptomics reveals regional and cellular vulnerability of the Sleep-Aging Convergence Signature during aging and Alzheimer’s disease progression",
    "Model linking chronic sleep disruption to proteostasis decline and reciprocal sleep-proteostasis interactions",
    "Phase-resolved and cosinor-style analyses of HSF1-associated genes in the human sleep-restriction dataset",
    "HSF1-associated acute sleep-deprivation responses exhibit preferential chronic attenuation in mouse brain",
    "Neuronal subtype-specific aging trajectories and inhibitory-neuron convergence network",
    "Human sleep-restriction and healthy-blood aging signatures share transcriptional features enriched for HSF1-associated proteostasis pathways",
    "Multi-omics convergence of aging and sleep restriction reveals systemic HSF1-associated proteostasis decline",
    "Aging attenuates acute SD-responsive transcriptional programs in mouse brain"
  ),
  panel_labels_in_deck = c(
    "a,b,c,d,e,f,g,h", "a,b,c,d", "a,b,c,d,e,f", "a,b,c,d", "a,b,c,d,e,f,g,h", "a,b",
    "a,b,c,e", "a,b,c,d,e,f,g,h", "a,b", "b,c,d,e,f,g", "a,d", "b,c,d,e"
  ),
  primary_r_modules = c(
    "01_sleep_restriction_phase_aware.R; 09_render_figures.R",
    "03_human_brain_aging.R; 09_render_figures.R",
    "03_human_brain_aging.R; 04_mouse_brain_sacs_definition.R; 09_render_figures.R",
    "06_multiomics_concordance.R; 09_render_figures.R",
    "07_spatial_seaad_sacs.R; 09_render_figures.R",
    "conceptual model in figure deck; source conclusions audited by upstream modules",
    "01_sleep_restriction_phase_aware.R; 10_supplementary_figures_public_source_check.R",
    "02_acute_sd_chronic_sf_attenuation.R; 10_supplementary_figures_public_source_check.R",
    "03_human_brain_aging.R; 10_supplementary_figures_public_source_check.R",
    "01_sleep_restriction_phase_aware.R; 03_human_brain_aging.R; 10_supplementary_figures_public_source_check.R",
    "06_multiomics_concordance.R; 10_supplementary_figures_public_source_check.R",
    "05_age_attenuated_sd_response.R; 10_supplementary_figures_public_source_check.R"
  ),
  source_scope = c(
    "GitHub_ready/Figure1",
    "GitHub_ready/Figure2",
    "GitHub_ready/Figure3",
    "GitHub_ready/Figure4",
    "GitHub_ready/Figure5",
    "Figure_aging-sleep_revised.pptx conceptual panel",
    "supplementary_figure_public_source_manifest.csv rows for Supplementary Fig. 1",
    "supplementary_figure_public_source_manifest.csv rows for Supplementary Fig. 2",
    "supplementary_figure_public_source_manifest.csv rows for Supplementary Fig. 3",
    "supplementary_figure_public_source_manifest.csv rows for Supplementary Fig. 4",
    "supplementary_figure_public_source_manifest.csv rows for Supplementary Fig. 5",
    "supplementary_figure_public_source_manifest.csv rows for Supplementary Fig. 6"
  ),
  stringsAsFactors = FALSE
)

main_dirs <- c(
  "GitHub_ready/Figure1", "GitHub_ready/Figure2", "GitHub_ready/Figure3",
  "GitHub_ready/Figure4", "GitHub_ready/Figure5", NA,
  rep(NA, 6)
)

manifest$source_table_count <- vapply(main_dirs, function(d) {
  if (is.na(d)) return(NA_integer_)
  count_files(d, "\\.(csv|xlsx|tsv)$")
}, integer(1))

manifest$plot_file_count <- vapply(main_dirs, function(d) {
  if (is.na(d)) return(NA_integer_)
  count_files(d, "\\.(png|pdf|tiff|tif|svg)$")
}, integer(1))

supp_manifest_path <- file.path(paths$tables, "supplementary_figure_public_source_manifest.csv")
supp_rows <- if (file.exists(supp_manifest_path)) read.csv(supp_manifest_path, stringsAsFactors = FALSE) else data.frame()
manifest$all_required_sources_present <- TRUE

for (i in seq_len(nrow(manifest))) {
  fig <- manifest$figure[i]
  if (startsWith(fig, "Fig.") && fig != "Fig. 6") {
    manifest$all_required_sources_present[i] <- isTRUE(manifest$source_table_count[i] > 0 && manifest$plot_file_count[i] > 0)
  } else if (startsWith(fig, "Supplementary")) {
    sub <- supp_rows[supp_rows$supplementary_figure == fig, , drop = FALSE]
    manifest$all_required_sources_present[i] <- nrow(sub) > 0 && all(sub$present)
  }
}

write.csv(manifest, file.path(paths$tables, "full_figure_source_script_manifest.csv"), row.names = FALSE)

if (!all(manifest$all_required_sources_present)) {
  write.csv(manifest[!manifest$all_required_sources_present, ],
            file.path(paths$tables, "full_figure_source_script_manifest_missing.csv"),
            row.names = FALSE)
  stop("One or more figure source/script manifest checks failed.", call. = FALSE)
}
