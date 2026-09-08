# Human and mouse aging source-data audit.
files <- list.files(paths$github_ready, pattern = "aging|Aging|trajectory|Trajectory", recursive = TRUE, full.names = TRUE)
files <- files[grepl("\\.(csv|xlsx)$", files)]

checks <- data.frame(
  analysis = "aging_trajectories",
  path = files,
  present = file.exists(files),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(paths$tables, "aging_source_check.csv"), row.names = FALSE)
if (nrow(checks) == 0) stop("No aging source tables found.", call. = FALSE)
