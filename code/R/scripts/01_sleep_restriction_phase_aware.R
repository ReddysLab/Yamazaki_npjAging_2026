# Phase-aware human sleep-restriction source-data audit.
files <- c(
  sr_gene_level = file.path(paths$github_ready, "Figure1", "Figure_1_SR_gene_level_results.csv"),
  reactome_ora = file.path(paths$github_ready, "Figure1", "Figure_1F_Reactome_ORA_full.csv")
)

checks <- data.frame(
  analysis = "human_sleep_restriction_phase_aware",
  file = names(files),
  path = unname(files),
  present = file.exists(files),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(paths$tables, "sleep_restriction_source_check.csv"), row.names = FALSE)
if (!all(checks$present)) stop("Missing sleep-restriction source table.", call. = FALSE)
