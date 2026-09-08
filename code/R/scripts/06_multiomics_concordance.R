# RNA/protein concordance and tissue Reactome network source-data audit.
files <- list.files(file.path(paths$github_ready, "Figure4"), pattern = "\\.(csv|xlsx)$", recursive = TRUE, full.names = TRUE)

checks <- data.frame(
  analysis = "multiomics_concordance",
  path = files,
  present = file.exists(files),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(paths$tables, "multiomics_source_check.csv"), row.names = FALSE)
if (nrow(checks) == 0) stop("No Figure 4 multi-omics source tables found.", call. = FALSE)
