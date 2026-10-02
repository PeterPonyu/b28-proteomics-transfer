# B28 panel compositor — independent grobs, no shared y-gutter.
#
# patchwork's default grid sizes the left gutter of a column to the widest
# tick-label / ylab block in that column.  A panel with short ticks then
# inherits a blank band between its y-axis title and its own axis.  That is
# the over-alignment defect.
#
# free(type = "label") and wrap_elements(full = ggplot) both failed here:
# the former clipped F3B's ylab at the cell edge, the latter still aligned
# sibling ylabs (F5A "AUROC" sat on F5B's "expected calibration error"
# gutter).  Passing a live ggplot into patchwork lets it re-open the gtable.
#
# This compositor tags each atomic ggplot, converts it to a grob, then
# arranges the grobs with cowplot::plot_grid(align = "none").  There is no
# sibling gutter.  Each ylab sits next to its own axis.
#
# Tags occupy an independent, fixed-height grob row. A top/right legend may
# change the panel's own layout, but never the row's A/B baseline or gutter.

B28_TAG_SIZE_PT <- 11
B28_TAG_CLEARANCE_PT <- 4
B28_TAG_BAND_PT <- B28_TAG_SIZE_PT + B28_TAG_CLEARANCE_PT

b28_tagged_cell <- function(panel, tag) {
  stopifnot(inherits(panel, "ggplot"), !inherits(panel, "patchwork"))
  panel <- panel + labs(tag = NULL) + theme(plot.tag = element_blank())
  grob <- ggplot2::ggplotGrob(panel)
  grob <- gtable::gtable_add_rows(grob, grid::unit(B28_TAG_BAND_PT, "pt"), pos = 0)
  gtable::gtable_add_grob(
    grob,
    grid::textGrob(tag, x = grid::unit(5, "pt"),
                   y = grid::unit(1, "npc") - grid::unit(4, "pt"),
                   just = c("left", "top"),
                   gp = grid::gpar(fontfamily = B1_FONT, fontface = "bold",
                                   fontsize = B28_TAG_SIZE_PT)),
    t = 1, l = 1, r = ncol(grob), clip = "off", name = "b28-panel-letter"
  )
}

b28_compose <- function(panels, widths = NULL, heights = NULL, ncol = NULL) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("cowplot required for b28_compose()", call. = FALSE)
  }
  stopifnot(is.list(panels), length(panels) >= 2L)
  n <- length(panels)
  for (i in seq_len(n)) {
    if (!inherits(panels[[i]], "ggplot") || inherits(panels[[i]], "patchwork")) {
      stop("b28_compose() requires atomic ggplot panels, never nested patchwork",
           call. = FALSE)
    }
  }
  if (is.null(ncol)) ncol <- n
  if (is.null(widths)) widths <- rep(1, ncol)
  # grob-first: a live ggplot handed to a grid compositor can be re-aligned.
  # Open cairo so ggplotGrob measures Arial (the default PostScript device
  # does not know the family and emits 50+ warnings plus wrong widths).
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp, width = 6.5, height = 4, family = B1_FONT)
  grobs <- Map(b28_tagged_cell, panels, LETTERS[seq_len(n)])
  grDevices::dev.off()
  unlink(tmp)
  cowplot::plot_grid(plotlist = grobs, nrow = ceiling(n / ncol), ncol = ncol,
                     align = "none", axis = "none",
                     rel_widths = widths,
                     rel_heights = if (is.null(heights)) 1 else heights)
}
