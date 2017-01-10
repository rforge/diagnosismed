#' Function to plot the mixed barplots of the frequencies of individuals with (1) and without (0) the targeted condition.
#' @name Mixed.Barplot
#' @description This plot function shows the frequencies of the two groups and their overlap in a single graph. This function is intended to be used for tests with 20 or less different test values
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param position.legend the position of the legend. Default"'topright'. See also \code{\link{legend}}.
#' @details The graph shows the overlapping barplots of the frequencies of the observed test results of two groups and their overlap. Many tests of intermediate quality have a considerable overlap. This graph allows for the visual inspection of the frequencies and their overlap. The zero scores (0) represent the norm group, 1 represent the targeted group (abnorm).
#' For studying the overlap, the mixed density plot is recommended: \code{plotMD}.
#' @return No Value returned.
#' @importFrom grDevices rgb
#' @importFrom graphics legend plot
#' @export
#' @seealso \code{\link{UI.ordinal}} for calculating possible Uncertain Intervals for short, ordened tests, \code{\link{plotMD}} for plotting mixed densities.
#' @examples
#' # A short test with 5 ordinal values
#' test0     = rep(1:5, times=c(33,6,6,11,2))
#' test1   = rep(1:5, times=c(3,2,2,11,33))
#' ref = c(rep(0, length(test0)), rep(1, length(test1)))
#' test = c(test0, test1)
#' mixed.barplot(ref, test, position.legend='top')
#' # please note that test score 4 has complete overlapping frequencies.

mixed.barplot <- function(ref, test, position.legend='topright') {

  check.data(ref, test, ordinal = T)

  minmax = c(min(test), max(test))
  x = factor(test[ref==0], levels = c(minmax[1]:minmax[2]))
  y = factor(test[ref==1], levels = c(minmax[1]:minmax[2]))
  tx = table(x)
  ty = table(y)
  m = max(c(tx, ty))
  r1 = 0
  g1 = 0
  b1 = 1
  a1 = 1 / 4
  r2 = 1
  g2 = 0
  b2 = 0
  a2 = 1 / 4
  a3 = 1 - 1 * (1 - a1) * (1 - a2)
  r3 = r1 * a1 / a3 + r2 * a2 * (1 - a1) / a3
  g3 = g1 * a1 / a3 + g2 * a2 * (1 - a1) / a3
  b3 = b1 * a1 / a3 + b2 * a2 * (1 - a1) / a3

  plot(
    x,
    main = 'Mixed Barplot',
    col = rgb(r1, g1, b1, a1),
    border = FALSE,
    xlab = 'Test | Predictor',
    ylim = c(0, m)
  )
  plot(y,
       add = T,
       col = rgb(r2, g2, b2, a2),
       border = FALSE)
  legend(
    position.legend,
    horiz = F,
    border = 0,
    box.col = 'transparent' ,
    legend = c('0', 'overlap', '1'),
    fill = c(rgb(r1, g1, b1, a1), rgb(r3, g3, b3, a3), rgb(r2, g2, b2, a2)),
    cex = .8
  )
}
