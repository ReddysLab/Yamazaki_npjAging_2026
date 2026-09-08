inventory <- inventory_source_tables()
write.csv(inventory, file.path(paths$tables, "source_table_inventory.csv"), row.names = FALSE)

required_files <- c(
  GSE39445_matrix = file.path(paths$source_data, "GSE39445", "GSE39445_series_matrix.txt.gz"),
  GSE39445_phase_map = file.path(paths$source_data, "GSE39445", "curated_metadata", "GSE39445_sample_phase_mapping.csv"),
  JenAge_aging_results = file.path(paths$source_data, "GSE134080_or_JenAge", "JenAge_01_healthy_blood_age_log2RPKM_results.csv"),
  Takasugi_multiomics = file.path(paths$source_data, "Takasugi_2024", "Figure_4_multiomics_long.csv.gz"),
  Takasugi_features = file.path(paths$source_data, "Takasugi_2024", "Figure_4_trajectory_features.csv"),
  Schaum_h5ad = Sys.getenv("NPJAGING_SCHAUM_H5AD", unset = file.path(paths$source_data, "GSE132040", "tabula_muris_senis.h5ad"))
)

required <- data.frame(
  item = names(required_files),
  path = unname(required_files),
  present = file.exists(required_files),
  stringsAsFactors = FALSE
)
write.csv(required, file.path(paths$tables, "required_source_file_check.csv"), row.names = FALSE)

if (!all(required$present)) {
  stop(
    "Missing required input(s): ",
    paste(required$item[!required$present], collapse = ", "),
    ". See source_data/SOURCE_MANIFEST.csv and set NPJAGING_SCHAUM_H5AD for the Schaum object.",
    call. = FALSE
  )
}
