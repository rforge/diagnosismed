#' Comparing diagnositic tests: a simple graphic using likelihood ratios.
#' @name LRgraph
#'
#' @description \code{LRgraph} graphically compares two or more (all of them with the first test) diagnostic tests with binary results through their likelihood ratios, based on the rationale that the predictive ability of a test is a more interesting characteristic than sensitivity or specificity alone. It is possible to see through the graph that if the tests with smaller sensitivity or specificity may have superior predictive ability, that is, increases the prediction ability with small sensitivity/specificity trade-off.
#'
#' @param tests A list of \code{\link{diagnosis}} objects composed by two or more tests. This object should be created listing two \code{\link{diagnosis}} objects as \code{list(mytest1, mytest2)}. One may insert as many tests as one wishes. See example.
#'
#' @param xlab,ylab Characters indicating the labels of the horizantal and vertical axis. These will be passed to \code{\link[graphics]{plot.default}}.
#'
#' @param xlim,ylim Vectors indicating the limits of the horizantal and vertical axis. These will be passed to \code{\link[graphics]{plot.default}}
#'
#' @param point.arg A \code{\link[base]{list}} of arguments to be passed to \code{\link[graphics]{points}}. Internally the arguments \code{x} and \code{y} are replaced by the values in the \code{tests} argument. These argumentes will be used recurssively to plot points, so the arguments in list should not have a length bigger than \code{tests}.
#'
#' @param lines.arg A \code{\link[base]{list}} of arguments to be passed to \code{\link[graphics]{segments}} to draw the likelihood ratio lines. Internally the argument \code{coef} is replaced by the values in the \code{tests} arguments.
#'
#' @param auto.legend Logical. If \code{TRUE}, it makes a legend of the graph.
#'
#' @param leg.arg A \code{\link[base]{list}} of arguments to be passed to \code{\link[graphics]{legend}} when \code{auto.legend = TRUE}. Internally, the \code{LRgraph} removes the \code{x} and \code{y}, and the duplicated arguments from \code{point.arg} and append them to draw the legend. Therefore, some graphical parameters should be be defined in \code{point.arg}. If one whised to override them, define them in leg.arg as well. It is very important that the \code{index.name} arguments froam all tests are defined when calling the \code{diagnosis} function, as it will be passed to \code{legend} function within \code{LRgraph} and may return an error if otherwise.
#'
#' @param grid Logical. If \code{TRUE}, it calls the \code{\link[graphics]{grid}} function with the deafult arguments.
#'
#' @param SeSp.lines Logical. If \code{TRUE}, a vertical and a horizontal lines will be drawn representing the Sensitivity and Specificity of the test.
#'
#' @param SeSp.lines.arg A \code{\link[base]{list}} of arguments to be passed to \code{\link[graphics]{abline}} when \code{SeSP.lines = TRUE}. Internally the arguments \code{v} and \code{h} are replaced by the values in the \code{tests} arguments.
#'
#' @param ... Other graphical parameters passed to \code{\link[graphics]{plot.default}}
#'
#' @details When a diagnostic test has both sensitivity and specificity higher than a competing test is easy to see that the former is superior than the later. However, sometimes a test may have superior sensitivity and inferior specificity (or the other way around). In this case, a good decision may be toward the test that have a better prediction ability. The graph visually helps the user to see and compare these abilities. The graph is very similar to the ROC graph. The vertical and horizontal axis have the same length as the ROC graph. However, the diagnostic tests are represented as dots instead of curves. The solid line passing through (0,0) is the likelihood ratio positive-line and the solid line passing through (1,1) is the likelihood ratio negative-line. Both negative and positive likelihood are numerically equivalent to the slopes of the solid lines. The solid lines split the graph into four areas (run the example). Also, there are dashed lines representing the sensitivity and specificity of the first test plotted. One may see that there are areas that a test may have superior sensitivity (or specificity) and yet the dot may be below the likelihood solid line. That is because the sensitivity / specificity trade-off is not reasonable, making the test with less predictive ability.
#'
#' @return Returns only a graph which is divided in four areas (by the solid lines representing the likelihood ratios of the firts test). The interpretation of the comparisons will depend on which area the competing tests will fall in. See and run the example to have the idea on how interpretation must be done.
#'
#' @references Biggerstaff, B.J. Comparing diagnostic tests: a simple graphic using likelihood ratios. Statistics in Medicine. 2000; 19(5):649-663
#'
#' @seealso \code{\link{diagnosis}}
#'
#' @examples
#' # Making tests with diagnosis function with different performances for comparison.
#' # mytest5 is the one which all others will be compared with.
#' # It is important to set the index.name argument to later ditinguish the tests at the plot.
#' mytest5 <- diagnosis(TP = 80, FN = 20, FP = 20, TN = 80, index.name = "Original")
#' mytest5
#'
#' # mytest1 has higher sensitivity and specificity.
#' # mytest1 is overall superior compared to mytest5.
#' mytest1 <- diagnosis(TP = 90, FN = 10, FP = 10, TN = 90, index.name = "Test 1")
#' mytest1
#'
#' # mytest2 has lower sensitivity but higher specificity.
#' # mytest2 is better to identify the presence of the target condition compared to mytest5.
#' mytest2 <- diagnosis(TP = 72, FN = 28, FP = 3, TN = 97, index.name = "Test 2")
#' mytest2
#'
#' # mytest3 has higher sensitivity but lower specificity.
#' # mytest3 is better to identify the absence of the target condition compared to mytest5.
#' mytest3 <- diagnosis(TP = 92, FN = 8, FP = 27, TN = 63, index.name = "Test 3")
#' mytest3
#'
#' # mytest4 has higher sensitivity and lower specificity.
#' # Nevertheless, mytest4 is overall inferior compared to mytest5.
#' mytest4 <- diagnosis(TP = 82, FN = 18, FP = 40, TN = 60, index.name = "Test 4")
#' mytest4
#'
#' # But that becomes clear only after ploting the tests.
#' allmytests <- list(mytest5, mytest1, mytest2, mytest3, mytest4)
#' LRgraph(allmytests, grid = FALSE)
#'
#' # The texts below are not part of the function but helps to understand the areas
#' text(x=.5, y =.5, labels ="Area 4: Overall inferior", col=gray(.6),cex=.8)
#' text(x=.5, y =1, labels ="Area 2: Absence", col=gray(.6),cex=.8)
#' text(x=.07, y =.68, labels ="Area 3: Presence", col=gray(.6),cex=.8, xpd = NA)
#' text(x=.1, y =1, labels ="Area 1: Overall superior", col=gray(.6),cex=.8, xpd = NA)
#'
#' rm(mytest5, mytest1, mytest2, mytest3, mytest4, allmytests)
#'
#' @import graphics
#' @import grDevices
#' @export
LRgraph <- function(tests, ..., xlab = "1 - Specificity", ylab = "Sensitivity",
                    ylim = c(0, 1), xlim = c(0, 1),
                    point.arg = list(cex = 2, pch = 1:length(tests), lwd = 2),
                    lines.arg = list(lwd = 2, col = "red"),
                    auto.legend = TRUE,
                    leg.arg = list(x = "top", inset = -.2, xpd = NA, ncol = 2, cex = .9, bty = "n"),
                    SeSp.lines = TRUE,
                    SeSp.lines.arg = list(lty = 2, col = gray(.4)),
                    grid = FALSE){
  if ( length(tests) < 2 ) {
    stop("'tests' should have at least two 2 diagnostic tests results.")
  }
  if ( !is.list(tests) ) {
    stop("'tests' is not a list.")
  }
  cla <- sapply(1:length(tests), function(i) class(tests[[i]]))
  if ( any(cla != "diagnosis") ) {
    stop("'tests' elements should all be from diagnosis function output.")
  }
  if ( !is.logical(auto.legend) ) {
    stop("'auto.lenged' is not logical." )
  }
  if ( !is.logical(SeSp.lines) ) {
    stop("'SeSp.lines' is not logical." )
  }
  if ( !is.logical(grid) ) {
    stop("'grid' is not logical." )
  }
  # Opening the plot
  plot(NA, NA, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  axis(1, at = seq(0, 1, .05), labels = FALSE, tcl = -.25)
  axis(2, at = seq(0, 1, .05), labels = FALSE, tcl = -.25)
  if (grid) { grid() }

  #Drawing the segments
  lines.arg$coef <- c(0, (tests[[1]]$Se / (1 - tests[[1]]$Sp)))
  do.call(abline, lines.arg)
  lines.arg$coef <- c(1 - ((1 - tests[[1]]$Se) / (1 - (1 - tests[[1]]$Sp))), (1 - tests[[1]]$Se) / (1 - (1 - tests[[1]]$Sp)))
  do.call(abline, lines.arg)

  # Drawing the poitns
  point.arg$x <- 1 - sapply(1:length(tests), function(i) tests[[i]]$Sp)
  point.arg$y <- sapply(1:length(tests), function(i) tests[[i]]$Se)
  do.call(points, point.arg)

  # Inserting the legend
  if (auto.legend) {
    point.arg$x <- NULL
    point.arg$y <- NULL
    point.arg <- point.arg[-which(point.arg %in% leg.arg)]
    leg.arg <- append(leg.arg, point.arg)
    leg.arg$legend <- sapply(1:length(tests), function(i) tests[[i]]$index.name)
    do.call(legend, leg.arg)
  }
  if (SeSp.lines) {
    SeSp.lines.arg$v <- 1 - tests[[1]]$Sp
    SeSp.lines.arg$h <- tests[[1]]$Se
    do.call(abline, SeSp.lines.arg)
  }
}
