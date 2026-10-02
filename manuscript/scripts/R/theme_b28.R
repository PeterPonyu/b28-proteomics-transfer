# B28-specific ggplot2 theme. Current scientific figures use Arial.
# cairo_pdf and ragg resolve and embed it by family name. Fail if unavailable,
# rather than silently mixing the figure family with a serif fallback.
#
# FILE-NAME NOTE: leftover-cycle 2026-09-09 renamed this file theme_b1.R ->
# theme_b28.R.  It is a SEPARATE, diverged B28 copy, not B1's theme.  Tests
# and make_b28_figures.R use this path. Legacy exported symbols
# (theme_b1, b1_save, B1_FONT, .b1_pick_font) remain local to this capsule
# for compatibility with its compositor. The historical JPR upload pack
# (manuscript.pdf, cover_letter_JPR.md, MANIFEST.txt, SHA256SUMS.txt,
# submission_checklist.md) does not include this file.
#
# Grid composition (y-label freedom + margin tags) lives in compose_b28.R;
# this file is type, colour, and the PNG/PDF writer.
suppressMessages({
  library(ggplot2)
  library(systemfonts)
  library(ragg)
  library(scales)
})

.b1_pick_font <- function() {
  fams <- unique(systemfonts::system_fonts()$family)
  if (!("Arial" %in% fams)) {
    stop("B28 figures require Arial; no silent fallback", call. = FALSE)
  }
  resolved <- systemfonts::match_fonts(rep("Arial", 2), weight = c("normal", "bold"))
  info <- systemfonts::font_info(path = resolved$path, index = resolved$index)
  if (any(info$family != "Arial") || info$bold[1] || !info$bold[2]) {
    stop("B28 requires exact Arial regular and bold faces; no substitution", call. = FALSE)
  }
  "Arial"
}
B1_FONT <- .b1_pick_font()

# F5/F6 are embedded at 0.92 linewidth. An 8 pt source would print at 7.36 pt.
# Keep every source text at least 8.3 pt, including geom_text (whose size is mm).
B28_MIN_TEXT_PT <- 8.3
b28_text_mm <- function(pt = B28_MIN_TEXT_PT) {
  stopifnot(is.numeric(pt), length(pt) == 1L, pt >= B28_MIN_TEXT_PT)
  pt / ggplot2::.pt
}

# Okabe-Ito colorblind-safe palette (identical hex to the matplotlib scripts).
# oi() returns UNNAMED hex so it can be used as an element inside a named values
# vector (c(label = oi("blue"))) without R contaminating the element name.
.OI <- c(blue = "#0072B2", orange = "#E69F00", verm = "#D55E00", green = "#009E73",
         purple = "#CC79A7", sky = "#56B4E9", yellow = "#F0E442", black = "#000000",
         grey = "#999999")
oi <- function(k) unname(.OI[[k]])
MODCOL <- c("Boltz-2" = oi("blue"), "Boltz-1" = oi("sky"), "Boltz-1x" = oi("green"),
            "Chai" = oi("orange"), "Protenix" = oi("verm"), "AF3" = oi("purple"),
            "AF3-noTmpl" = oi("black"), "RFAA" = oi("grey"))

theme_b1 <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = B1_FONT) +
    theme(
      text = element_text(family = B1_FONT, face = "plain"),
      axis.title = element_text(size = max(base_size, B28_MIN_TEXT_PT)),
      axis.text = element_text(size = max(base_size - 1, B28_MIN_TEXT_PT), colour = "black"),
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = max(base_size - 2, B28_MIN_TEXT_PT)),
      legend.key.size = unit(0.9, "lines"),
      # No legend frame and no legend background: transparent fill + NA colour so an
      # inside-plot legend never draws a box or backing panel over the plotting area.
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      plot.title = element_text(size = max(base_size - 1, B28_MIN_TEXT_PT), face = "plain", hjust = 0.5,
                                family = B1_FONT, margin = margin(b = 5)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      # Single panels retain margin-band tags. compose_b28.R strips these
      # and draws independent letter grobs on a common row baseline.
      plot.tag = element_text(size = 11, face = "bold", family = B1_FONT,
                              hjust = 0, vjust = 1, margin = margin(b = 3, r = 4)),
      plot.tag.position = "topleft",
      plot.tag.location = "margin",
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin = margin(4, 5, 4, 5),
      # Each independent panel keeps this gap from its own ylab to its ticks;
      # it never inherits a wider neighbouring panel's gutter.
      axis.title.y = element_text(margin = margin(r = 3))
    )
}

# Save a ggplot/patchwork object as both a 300-dpi PNG (ragg) and a vector PDF
# (cairo_pdf, fonts embedded) into the figures/ directory.
b1_save <- function(plot, stem, w, h, figdir) {
  # Render at the printed text width (6.5in = \linewidth on letterpaper/1in margins) so the
  # figure is shown 1:1 by \includegraphics[width=\linewidth] and its text prints at base_size
  # (~10pt), matched to the 11pt body. Scale the authored height by 6.5/w so the figure keeps its
  # authored aspect ratio at column width -- otherwise wide multi-panel figures forced to 6.5in
  # wide while keeping their authored height render as tall, skinny panels.
  h <- h * 6.5 / w
  w <- 6.5
  png <- file.path(figdir, paste0(stem, ".png"))
  pdf <- file.path(figdir, paste0(stem, ".pdf"))
  ragg::agg_png(png, width = w, height = h, units = "in", res = 300, background = "white")
  print(plot); invisible(grDevices::dev.off())
  grDevices::cairo_pdf(pdf, width = w, height = h, family = B1_FONT)
  print(plot); invisible(grDevices::dev.off())
  cat("wrote", basename(png), "and", basename(pdf), "\n")
}
