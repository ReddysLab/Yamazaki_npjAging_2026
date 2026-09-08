#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

if (!exists("paths", inherits = TRUE)) {
  args <- commandArgs(FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)][1])
  script_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg, mustWork = FALSE)) else getwd()
  source(file.path(dirname(script_dir), "config.R"))
}

multiomics_path <- file.path(paths$source_data, "Takasugi_2024", "Figure_4_multiomics_long.csv.gz")
original_sr_path <- file.path(paths$source_data, "GSE39445", "Figure_1E_1F_gene_level_phase_adjusted_limma.csv")
participant_sr_path <- file.path(paths$tables, "moller_levet_participant_aware", "MollerLevet_participant_aware_phase_cosinor_gene_results.csv")
out_dir <- file.path(paths$tables, "multiomics_aging_definition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(multiomics_path), file.exists(original_sr_path), file.exists(participant_sr_path))

message("Calculating replicate-level tissue/assay aging slopes")
dt <- fread(multiomics_path, select = c("value", "age", "tissue", "assay", "feature_upper"))
dt <- dt[is.finite(value) & is.finite(age) & !is.na(feature_upper) & feature_upper != ""]
dt[, `:=`(expr = log2(value + 1), gene = toupper(feature_upper), assay = as.character(assay), tissue = as.character(tissue))]

slope_stats <- dt[, {
  n_obs <- .N
  n_age <- uniqueN(age)
  if (n_age < 3L || n_obs < 6L || var(age) == 0 || var(expr) == 0) {
    .(n = n_obs, n_age = n_age, slope = NA_real_, p = NA_real_)
  } else {
    fit <- tryCatch(summary(lm(expr ~ age)), error = function(e) NULL)
    if (is.null(fit)) {
      .(n = n_obs, n_age = n_age, slope = NA_real_, p = NA_real_)
    } else {
      co <- fit$coefficients
      .(n = n_obs, n_age = n_age, slope = unname(co["age", "Estimate"]), p = unname(co["age", "Pr(>|t|)"]))
    }
  }
}, by = .(tissue, assay, gene)]
slope_stats[, fdr := p.adjust(p, method = "BH"), by = .(tissue, assay)]

rna <- slope_stats[assay == "RNA", .(tissue, gene, rna_slope = slope, rna_p = p, rna_fdr = fdr)]
protein <- slope_stats[assay %in% c("WTL", "LSF")]
protein_gene <- protein[, .(
  has_protein = any(is.finite(slope)),
  protein_any_down = any(slope < 0, na.rm = TRUE),
  protein_layers_down = paste(sort(unique(assay[slope < 0])), collapse = ";")
), by = .(tissue, gene)]

aging <- merge(rna, protein_gene, by = c("tissue", "gene"), all = FALSE)
aging[, aging_decline := is.finite(rna_slope) & rna_slope < 0 & has_protein & protein_any_down]

# The Fig. 4/Supplementary Fig. 5 universe is fixed to genes tested in the
# original GSE39445 phase-adjusted analysis and measured in RNA plus protein.
original_sr_genes <- read_csv(original_sr_path, show_col_types = FALSE) |>
  transmute(gene = toupper(gene_symbol)) |>
  filter(!is.na(gene), gene != "") |>
  distinct() |>
  pull(gene)

participant_sr <- read_csv(participant_sr_path, show_col_types = FALSE) |>
  transmute(
    gene = toupper(gene),
    phase_logFC = participant_aware_phase_logFC,
    phase_p = participant_aware_phase_p,
    phase_FDR = participant_aware_phase_FDR,
    delta_MESOR = delta_mesor,
    MESOR_p = mesor_p,
    MESOR_FDR = mesor_FDR,
    participant_aware_SR_down = phase_logFC < 0 & phase_FDR < 0.05 & delta_MESOR < 0 & MESOR_FDR < 0.05
  ) |>
  arrange(phase_p, MESOR_p) |>
  distinct(gene, .keep_all = TRUE)

membership <- as_tibble(aging) |>
  filter(has_protein, gene %in% original_sr_genes) |>
  left_join(participant_sr, by = "gene") |>
  mutate(
    participant_aware_SR_down = coalesce(participant_aware_SR_down, FALSE),
    convergence = aging_decline & participant_aware_SR_down
  )

fisher_one <- function(df) {
  a <- df$participant_aware_SR_down
  b <- df$aging_decline
  mat <- matrix(c(sum(a & b), sum(a & !b), sum(!a & b), sum(!a & !b)), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  tibble(
    universe_n = nrow(df),
    sr_n = sum(a),
    aging_n = sum(b),
    overlap_n = sum(a & b),
    sr_only_n = sum(a & !b),
    aging_only_n = sum(!a & b),
    fisher_or = unname(ft$estimate),
    fisher_p = ft$p.value
  )
}

summary <- membership |>
  group_by(tissue) |>
  group_modify(~fisher_one(.x)) |>
  ungroup() |>
  mutate(
    sr_threshold = "participant-aware phase logFC < 0 and FDR < 0.05; delta MESOR < 0 and MESOR FDR < 0.05",
    aging_threshold = "replicate-level log2(value+1) RNA age slope < 0 and at least one WTL/LSF protein age slope < 0",
    background_universe = "genes tested in the original GSE39445 phase-adjusted analysis and measured in tissue RNA plus protein"
  )

fwrite(slope_stats, file.path(out_dir, "multiomics_replicate_level_age_slope_stats.csv"))
write_csv(as_tibble(aging), file.path(out_dir, "multiomics_RNAprotein_aging_decline_flags.csv"))
write_csv(membership, file.path(out_dir, "F4_SF5_participant_aware_tissue_membership.csv"))
write_csv(summary, file.path(out_dir, "F4c_SF5_participant_aware_overlap_summary.csv"))

brain_n <- summary$aging_n[summary$tissue == "Brain"]
if (length(brain_n) != 1L || brain_n != 1634L) stop("Brain RNA+protein aging-decline count did not reproduce n=1634; observed ", brain_n)
message("Verified brain RNA+protein aging-decline set: n=", brain_n)
