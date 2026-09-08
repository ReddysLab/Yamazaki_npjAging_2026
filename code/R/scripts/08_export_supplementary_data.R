# Record the Supplementary Data workbook generated for submission.
workbook <- file.path(project_root, "Supplementary_Data_1_full_analysis_outputs.xlsx")
summary <- data.frame(
  workbook = workbook,
  present = file.exists(workbook),
  bytes = if (file.exists(workbook)) file.info(workbook)$size else NA_real_,
  note = "Workbook assembled from copied source/result CSV tables; see workbook README and Table_Index sheets.",
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(paths$tables, "supplementary_workbook_check.csv"), row.names = FALSE)
if (!summary$present) stop("Supplementary Data workbook not found.", call. = FALSE)
