# Ready for QC
#' Save pheatmap fxn

save_pheatmap_pdf <- function(x, filename, width=12, height=25) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  cairo_pdf(filename, width=width, height=height, family = "NotoSans-Condensed", fallback_resolution = 300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf2 <- function(x, filename, width=8, height=11) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  cairo_pdf(filename, width=width, height=height, family = "NotoSans-Condensed", fallback_resolution = 300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
