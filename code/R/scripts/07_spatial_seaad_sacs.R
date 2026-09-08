# Spatial transcriptomics and SEA-AD SACS source-data audit.
files <- list.files(file.path(paths$github_ready, "Figure5"), pattern = "\\.(csv|xlsx)$", recursive = TRUE, full.names = TRUE)
files <- c(files, list.files(paths$sacs_outputs, pattern = "SEA|SACS_Fig5|spatial", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
files <- unique(files[grepl("\\.(csv|xlsx)$", files)])

checks <- data.frame(
  analysis = "spatial_seaad_sacs",
  path = files,
  present = file.exists(files),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(paths$tables, "spatial_seaad_source_check.csv"), row.names = FALSE)
if (nrow(checks) == 0) stop("No spatial/SEA-AD source tables found.", call. = FALSE)
