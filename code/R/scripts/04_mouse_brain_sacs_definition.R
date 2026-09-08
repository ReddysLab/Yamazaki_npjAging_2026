# SACS definition source-data audit.
files <- c(
  all_gene_scores = file.path(paths$sacs_outputs, "SACS_acuteSD_chronicSF_attenuation_all_gene_scores.csv"),
  convergent_genes = file.path(paths$sacs_outputs, "SACS_2267_convergent_genes_attenuation_scores.csv")
)
fallback <- list.files(paths$sacs_outputs, pattern = "SACS.*(definition|final|convergent|attenuation).*\\.csv$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
files <- unique(c(files[file.exists(files)], fallback))

checks <- data.frame(
  analysis = "mouse_brain_sacs_definition",
  path = files,
  present = file.exists(files),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(paths$tables, "sacs_definition_source_check.csv"), row.names = FALSE)
if (length(files) == 0) stop("No SACS definition source tables found.", call. = FALSE)
