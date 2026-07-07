# =============================================================================
# atac_theme.R — shared publication theme, colour-blind-safe palettes, and
# figure-saving helpers for every R figure script in the pipeline.
# =============================================================================
# Source this from other scripts:  source("workflow/scripts/atac_theme.R")
# Design goals (Nature/Cell-ready): vector PDF + high-DPI PNG, minimal chrome,
# colour-blind-safe categorical hues assigned in a fixed order, a perceptually
# uniform sequential ramp, and a balanced blue-white-red diverging ramp for
# signed quantities (log2FC, z-scores). Palettes are validated colour-blind-safe.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# cairo_pdf gives the nicest font/anti-alias output but needs an R built with
# cairo; fall back to the base pdf device otherwise so figures never hard-fail.
.atac_has_cairo <- isTRUE(tryCatch(capabilities("cairo"), error = function(e) FALSE))
.atac_pdf_device <- function() if (.atac_has_cairo) grDevices::cairo_pdf else "pdf"
.atac_pdf <- function(file, width, height, onefile = FALSE) {
  if (.atac_has_cairo) grDevices::cairo_pdf(file, width = width, height = height,
                                            onefile = onefile)
  else grDevices::pdf(file, width = width, height = height, onefile = onefile)
}

# ---- Colour-blind-safe categorical palette (fixed slot order) ---------------
# Assigned in order; never cycled. Extends to a 12-hue fallback for many groups.
ATAC_CATEGORICAL <- c(
  "#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948",
  "#e87ba4", "#eb6834", "#6d4b9f", "#00a2b3", "#8c8c00", "#a6611a"
)

# Sequential ramp (magnitude; light -> dark blue). Falls back gracefully.
ATAC_SEQ <- colorRampPalette(c("#f5f9ff", "#9ec5f4", "#3987e5", "#184f95", "#0d366b"))

# Diverging ramp (signed; blue - neutral gray - red). Symmetric around zero.
ATAC_DIVERGING <- colorRampPalette(
  c("#2166ac", "#4393c3", "#92c5de", "#f0efec", "#f4a582", "#d6604d", "#b2182b")
)

# Status/highlight colours for up/down/ns categories.
ATAC_STATUS <- c(up = "#d6604d", down = "#2166ac", ns = "#b9b7b0")

# ---- Publication ggplot theme -----------------------------------------------
theme_atac <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6e5df", linewidth = 0.3),
      panel.border     = element_rect(colour = "#9a988f", fill = NA, linewidth = 0.5),
      axis.ticks       = element_line(colour = "#9a988f", linewidth = 0.4),
      axis.text        = element_text(colour = "#2b2b2b"),
      axis.title       = element_text(colour = "#0b0b0b"),
      plot.title       = element_text(face = "bold", hjust = 0, size = rel(1.1)),
      plot.subtitle    = element_text(colour = "#52514e", size = rel(0.9)),
      plot.caption     = element_text(colour = "#898781", size = rel(0.75)),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "#f2f1ec", colour = "#9a988f", linewidth = 0.4),
      strip.text       = element_text(face = "bold")
    )
}

# Map condition levels -> stable colours (named vector), recycling the palette.
# Sorted assignment so a condition keeps its colour across every figure (and
# matches plot_fragment_sizes.py, which sorts likewise).
atac_condition_colours <- function(levels) {
  levels <- sort(as.character(unique(levels)))
  n <- length(levels)
  pal <- if (n <= length(ATAC_CATEGORICAL)) ATAC_CATEGORICAL[seq_len(n)]
         else grDevices::colorRampPalette(ATAC_CATEGORICAL)(n)
  stats::setNames(pal, levels)
}

# ---- Figure saver: writes a vector PDF and, optionally, a MultiQC-ready PNG --
# `mqc = TRUE` also writes "<stem>_mqc.png" so MultiQC auto-embeds the figure as
# its own report section (no extra config needed).
save_fig <- function(plot, path_pdf, width = 7, height = 5.5, dpi = 320,
                     png = TRUE, mqc = FALSE) {
  ggsave(path_pdf, plot, width = width, height = height,
         device = .atac_pdf_device(), limitsize = FALSE)
  if (png) {
    png_path <- sub("\\.pdf$", ".png", path_pdf)
    ggsave(png_path, plot, width = width, height = height, dpi = dpi,
           limitsize = FALSE)
  }
  if (mqc) {
    mqc_path <- sub("\\.pdf$", "_mqc.png", path_pdf)
    ggsave(mqc_path, plot, width = width, height = height, dpi = 200,
           limitsize = FALSE)
  }
  invisible(path_pdf)
}

# Save a base-graphics / pheatmap object (a function or grob) to PDF (+PNG).
save_base_fig <- function(draw_fn, path_pdf, width = 7, height = 6,
                          dpi = 320, png = TRUE, mqc = FALSE) {
  .atac_pdf(path_pdf, width = width, height = height)
  draw_fn()
  grDevices::dev.off()
  render_png <- function(p) {
    grDevices::png(p, width = width, height = height, units = "in", res = dpi)
    draw_fn()
    grDevices::dev.off()
  }
  if (png) render_png(sub("\\.pdf$", ".png", path_pdf))
  if (mqc) render_png(sub("\\.pdf$", "_mqc.png", path_pdf))
  invisible(path_pdf)
}

# Write a tiny "not enough data" placeholder so declared rule outputs exist even
# when a figure cannot be drawn (e.g. no significant peaks).
placeholder_fig <- function(path_pdf, msg, mqc = FALSE) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 5,
                           colour = "#52514e") +
    theme_void()
  save_fig(p, path_pdf, width = 6, height = 4, mqc = mqc)
}
