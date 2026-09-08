# Export Supplementary Fig. 6e from the finalized submission workbook.

suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
})

workbook_path <- Sys.getenv(
  "NPJAGING_SUPPLEMENTARY_OUTPUT",
  unset = file.path(project_root, "Supplementary_Data_1_round2.xlsx")
)
workbook_path <- require_file(workbook_path, "final Supplementary Data 1 workbook")

summary_row <- readWorkbook(
  workbook_path, sheet = "SF6e", rows = 5L, cols = 1L:12L,
  colNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE
)
values <- as.numeric(summary_row[1L, 1L:9L])
names(values) <- c(
  "universe", "sacs", "aging_attenuated", "overlap", "sacs_only",
  "aging_attenuated_only", "background_only", "odds_ratio", "p_value"
)

stopifnot(
  values[["sacs"]] == 553,
  values[["aging_attenuated"]] == 1404,
  values[["overlap"]] == 98,
  values[["sacs_only"]] == 455,
  values[["aging_attenuated_only"]] == 1306,
  abs(values[["odds_ratio"]] - 2.269942) < 1e-6,
  abs(values[["p_value"]] - 3.37301e-11) < 1e-15
)

circle <- function(cx, cy = 0, radius = 1.12, group) {
  theta <- seq(0, 2 * pi, length.out = 500)
  data.frame(
    x = cx + radius * cos(theta),
    y = cy + radius * sin(theta),
    group = group
  )
}

circles <- rbind(
  circle(0, group = "Aging-attenuated SD response"),
  circle(1.28, group = "SACS")
)

panel <- ggplot() +
  geom_polygon(
    data = circles, aes(x, y, group = group, fill = group, colour = group),
    linewidth = 0.8, alpha = 0.14
  ) +
  scale_fill_manual(values = c(
    "Aging-attenuated SD response" = "#4B9CD3",
    "SACS" = "#5DBB63"
  )) +
  scale_colour_manual(values = c(
    "Aging-attenuated SD response" = "#1F78B4",
    "SACS" = "#2A9D45"
  )) +
  annotate("text", x = -0.36, y = 0, label = "1306", size = 4.2) +
  annotate("text", x = 0.64, y = 0, label = "98", size = 4.2) +
  annotate("text", x = 1.64, y = 0, label = "455", size = 4.2) +
  annotate(
    "text", x = -0.48, y = 1.43,
    label = "Aging-attenuated\nSD response\nn=1404",
    lineheight = 0.9, fontface = "bold", size = 3.2
  ) +
  annotate(
    "text", x = 1.76, y = 1.43,
    label = "SACS\nn=553",
    lineheight = 0.9, fontface = "bold", size = 3.2
  ) +
  annotate(
    "label", x = 0.64, y = -1.42,
    label = "Fisher~OR==2.27*';'~~italic(P)==3.37%*%10^{-11}",
    parse = TRUE, size = 3.0, linewidth = 0.25, label.padding = grid::unit(0.12, "lines"),
    fill = "white"
  ) +
  coord_fixed(xlim = c(-1.45, 2.73), ylim = c(-1.72, 1.78), clip = "off") +
  guides(fill = "none", colour = "none") +
  theme_void(base_size = 10) +
  theme(plot.margin = margin(5, 7, 5, 7))

output_dir <- Sys.getenv(
  "NPJAGING_SF6E_OUTPUT_DIR",
  unset = file.path(paths$figures, "supplementary_figure_6")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(output_dir, "SF6e_SACS_aging_attenuated_SD_overlap_exact_Fisher.pdf")
png_path <- file.path(output_dir, "SF6e_SACS_aging_attenuated_SD_overlap_exact_Fisher.png")
ggsave(pdf_path, panel, width = 4.5, height = 3.5, units = "in", device = "pdf", useDingbats = FALSE)
ggsave(png_path, panel, width = 4.5, height = 3.5, units = "in", dpi = 600, bg = "white")

message("Saved ", pdf_path)
message("Saved ", png_path)
