args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
source(file.path(script_dir, "config.R"))

dir.create(paths$processed, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)

check_packages()

modules <- c(
  "00_inventory_inputs.R",
  "01b_sleep_restriction_participant_aware_phase_model.R",
  "01c_propagate_participant_aware_moller_levet_dependencies.R",
  "01d_update_participant_aware_fig1f_fig1g.R",
  "01e_update_participant_aware_fig3de.R",
  "04a_multiomics_rna_protein_aging_decline.R",
  "01j_export_bothFDR_human_SRdown_matched_panels.R",
  "08_export_round2_supplementary_data.R",
  "09_export_sf6e_exact_fisher_venn.R"
)

for (module in modules) {
  message("Running ", module)
  source(file.path(script_dir, "scripts", module), local = TRUE)
}

writeLines(capture.output(sessionInfo()), file.path(project_root, "code", "results", "sessionInfo.txt"))
