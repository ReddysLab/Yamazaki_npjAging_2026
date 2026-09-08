# Independent young-vs-old acute SD inducibility analysis audit.
files <- list.files(paths$source_materials, pattern = "GSE128770|aging_attenuated_SD|age_attenuated|brain_SR_aging_overlap", recursive = TRUE, full.names = TRUE)
files <- files[grepl("\\.(csv|xlsx)$", files)]

checks <- data.frame(
  analysis = "young_old_sd_inducibility",
  path = files,
  present = file.exists(files),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(paths$tables, "age_attenuated_sd_response_source_check.csv"), row.names = FALSE)
if (nrow(checks) == 0) stop("No young-vs-old SD response source tables found.", call. = FALSE)
