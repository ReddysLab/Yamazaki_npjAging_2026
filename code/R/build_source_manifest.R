args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
source(file.path(script_dir, "config.R"))

files <- list.files(paths$source_data, recursive = TRUE, full.names = TRUE, all.files = FALSE)
files <- files[file.info(files)$isdir %in% FALSE]
files <- files[basename(files) != "SOURCE_MANIFEST.csv"]

rel <- sub(paste0("^", normalizePath(paths$source_data), "/"), "", normalizePath(files))
dataset <- sub("/.*$", "", rel)
accession <- unname(c(
  GSE39445 = "GSE39445",
  GSE113754 = "GSE113754",
  PRJNA757396 = "PRJNA757396",
  GSE128770 = "GSE128770",
  GSE132040 = "GSE132040",
  GSE134080_or_JenAge = "JenAge; GSE103232; GSE75337",
  Mongrain = "Mongrain et al. published supplementary data",
  Takasugi_2024 = "Takasugi et al. 2024",
  SEA_AD_reference_or_links = "SEA-AD middle temporal gyrus"
)[dataset])
accession[is.na(accession)] <- dataset[is.na(accession)]

sha256 <- vapply(files, function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1]])
}, character(1))

manifest <- data.frame(
  dataset = dataset,
  accession_or_reference = accession,
  bundled_path = rel,
  bytes = file.info(files)$size,
  sha256 = sha256,
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$dataset, manifest$bundled_path), ]
write.csv(manifest, file.path(paths$source_data, "SOURCE_MANIFEST.csv"), row.names = FALSE, na = "")
message("Wrote source_data/SOURCE_MANIFEST.csv with ", nrow(manifest), " bundled files.")
