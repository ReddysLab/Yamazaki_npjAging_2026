args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
this_file <- if (!is.na(file_arg)) normalizePath(file_arg, mustWork = FALSE) else NA_character_
code_root <- if (!is.na(this_file)) normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE) else normalizePath(getwd(), mustWork = FALSE)
project_root <- Sys.getenv("NPJAGING_PROJECT_ROOT", unset = normalizePath(file.path(code_root, ".."), mustWork = FALSE))

paths <- list(
  code = file.path(project_root, "code"),
  source_data = file.path(project_root, "source_data"),
  curated_metadata = file.path(project_root, "curated_metadata"),
  processed = file.path(project_root, "code", "results", "processed"),
  figures = file.path(project_root, "code", "results", "figures"),
  tables = file.path(project_root, "code", "results", "tables"),
  final_dual_fdr = file.path(project_root, "code", "results", "final_dual_fdr")
)

# Backward-compatible aliases used by the non-SR analysis modules. New code
# should use source_data/curated_metadata directly.
paths$source_materials <- paths$source_data
paths$github_ready <- file.path(paths$tables, "final_source_data")
paths$sacs_outputs <- file.path(paths$source_data, "SACS")
paths$aging_revision_upload <- file.path(paths$source_data, "revision_inputs")
paths$human_healthy_blood <- file.path(paths$source_data, "GSE134080_or_JenAge")
paths$raw_aging <- paths$source_data

thresholds <- list(
  acute_sd_logfc_min = 0,
  acute_sd_fdr_max = 0.05,
  chronic_sf_attenuation_min = 0,
  aging_slope_max = 0,
  human_sr_phase_fdr_max = 0.05,
  human_sr_mesor_fdr_max = 0.05,
  reactome_fdr_max = 0.05,
  reactome_nominal_p_max = 0.05
)

required_packages <- c(
  "GEOquery", "Biobase", "limma", "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble",
  "data.table", "rhdf5", "ggplot2", "ggrepel", "ggforce", "igraph", "ggraph", "openxlsx", "janitor",
  "msigdbr", "clusterProfiler", "ReactomePA", "org.Mm.eg.db", "jsonlite"
)

check_packages <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

inventory_source_tables <- function() {
  roots <- c(paths$source_data, paths$curated_metadata)
  files <- unlist(lapply(roots[file.exists(roots)], function(root) {
    list.files(root, pattern = "\\.(csv|tsv|xlsx)$", recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  data.frame(
    file = files,
    relative_path = sub(paste0("^", project_root, "/?"), "", files),
    bytes = file.info(files)$size,
    stringsAsFactors = FALSE
  )
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing required input for ", label, ": ", path, call. = FALSE)
  }
  path
}
