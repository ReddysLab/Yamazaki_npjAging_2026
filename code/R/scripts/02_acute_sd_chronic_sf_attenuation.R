# Acute SD induction and chronic SF attenuation source-data audit.
files <- c(
  attenuation_stats = file.path(paths$sacs_outputs, "SACS_acuteSD_up_x_positive_attenuation_venn_stats.csv"),
  attenuation_scores = file.path(paths$sacs_outputs, "SACS_acuteSD_chronicSF_attenuation_all_gene_scores.csv")
)

checks <- data.frame(
  analysis = "acute_sd_chronic_sf_attenuation",
  file = names(files),
  path = unname(files),
  present = file.exists(files),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(paths$tables, "acute_chronic_attenuation_source_check.csv"), row.names = FALSE)
if (!all(checks$present)) stop("Missing acute/chronic attenuation source table.", call. = FALSE)
