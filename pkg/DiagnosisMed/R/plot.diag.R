#' @rdname diagnosis
#' @import grDevices
#' @export
plot.diagnosis <- function(x, type = c("nomogram","roc"), ...,
                           xlab = "Pre-test probability", ylab = "Post-test proabbility",
                           lines.arg = list(col = "red", lwd = 3),
                           grid = FALSE, auto.shade = TRUE,
                           shade.arg = list(border = par("bg"), col = gray(.8)),
                           auto.legend = TRUE) {
  if (type != "nomogram" && type != "roc") {
    stop("'type' must be either 'nomogram' or 'roc.'")
  }
  if (!is.logical(grid)) {
    stop("'grid' must logical.")
  }
  if (!is.logical(auto.shade)) {
    stop("'auto.shade' must logical.")
  }
  if (!is.logical(auto.legend)) {
    stop("'auto.legend' must logical.")
  }
  if (type[1] == "nomogram") {
    pre.test <- seq(0, 1, by = .01)
    post.test <- ((pre.test / (1 - pre.test)) * x$PLR) / (1 + ((pre.test / (1 - pre.test)) * x$PLR))
    post.test.inf <- ((pre.test / (1 - pre.test)) * x$PLR.inf.cl) / (1 + ((pre.test / (1 - pre.test)) * x$PLR.inf.cl))
    post.test.sup <- ((pre.test / (1 - pre.test)) * x$PLR.sup.cl) / (1 + ((pre.test / (1 - pre.test)) * x$PLR.sup.cl))
    output <- data.frame(pre.test = pre.test, post.test = post.test, post.test.inf = post.test.inf, post.test.sup = post.test.sup)
    plot(NA, NA, type = "n", xlim = c(0,1), ylim = c(0, 1), xlab = xlab, ylab = ylab)
    axis(1, at = seq(0, 1, .05), labels = FALSE, tcl = -.25)
    axis(2, at = seq(0, 1, .05), labels = FALSE, tcl = -.25)
    if (grid) { grid() }
    if (auto.shade) {
      shade.arg$x <- c(pre.test[-101], sort(pre.test[-101], decreasing = TRUE))
      shade.arg$y <- c(post.test.inf[-101], sort(post.test.sup[-101], decreasing = TRUE))
      do.call(polygon, shade.arg)
    }
    lines.arg$y <- output$post.test
    lines.arg$x <- output$pre.test
    do.call(lines, lines.arg)
    if (auto.legend) {
      legend("bottomright", legend = paste("Positive Likelihood Ratio:", sprintf("%.2f",x$PLR)), bty= 'n')
    }
  }
  if (type[1] == "roc") {
    plot(NA, NA, type = "n", xlim = c(0,1), ylim = c(0, 1), xlab = xlab, ylab = ylab)
    axis(1, at = seq(0, 1, .05), labels = FALSE, tcl = -.25)
    axis(2, at = seq(0, 1, .05), labels = FALSE, tcl = -.25)
    if (grid) { grid() }
    lines.arg$x0 <- 0
    lines.arg$y0 <- 0
    lines.arg$x1 <- 1 - x$Sp
    lines.arg$y1 <- x$Se
    do.call(segments, lines.arg)
    lines.arg$x0 <- 1 - x$Sp
    lines.arg$y0 <- x$Se
    lines.arg$x1 <- 1
    lines.arg$y1 <- 1
    do.call(segments, lines.arg)
    # plot(1 - x$Sp, x$Se)
    # segments(0, 0, 1 - x$Sp, x$Se, col = "red")
    # segments(1 - x$Sp, x$Se, 1, 1, col = "red")
    if (auto.legend) {
      legend("bottomright", legend = paste("Area under the ROC curve:", sprintf("%.3f",x$AUC)), bty= 'n')
    }
  }
}
