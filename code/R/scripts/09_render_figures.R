# Final figure source-data and plot-file audit.
plot_files <- list.files(paths$github_ready, pattern = "\\.(png|pdf|tiff|tif|svg)$", recursive = TRUE, full.names = TRUE)
table_files <- list.files(paths$github_ready, pattern = "\\.(csv|xlsx|tsv)$", recursive = TRUE, full.names = TRUE)

summary <- data.frame(
  type = c(rep("plot", length(plot_files)), rep("source_table", length(table_files))),
  path = c(plot_files, table_files),
  bytes = file.info(c(plot_files, table_files))$size,
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(paths$tables, "figure_asset_inventory.csv"), row.names = FALSE)
if (length(plot_files) == 0 || length(table_files) == 0) {
  stop("Missing figure plots or source tables under GitHub_ready.", call. = FALSE)
}
