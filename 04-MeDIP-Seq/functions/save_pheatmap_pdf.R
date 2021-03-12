#' Save pheatmap fxn

save_pheatmap_pdf <- function(x, filename, width = 8, height = 8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  cairo_pdf(filename, width = width, height = height, family = "NotoSans_Condensed")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
