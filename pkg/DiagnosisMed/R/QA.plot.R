#' Quality assesment or Risk of Bias plot for systematic reviews
#'
#' @description This function produces a plot with the classification of the risk of bias (or quality assessment) of individual studies included in a systematic review. Originally, it was developed inspired in QUADAS-2, but it could work equally well with any risk of bias assessment tool such as MINORT, or Cochrane`s, as long the tool classifies the studies into one of three classes (or less).
#'
#' @param tab A \code{data.frame}, with the studies IDs and all risk of bias assessment. The study ID must always be in the first column on the left and the remaining columns should have three levels representing the risk of bias assessment (see example).
#'
#' @param class A chracter vector of length three, with the classes of the quality assessment to match the classes in the \code{tab} argument.
#'
#' @param mar The margin parameter that will be passed to \code{\link[graphics]{par}}.
#'
#' @param study.arg A list of arguments to format the studies IDs. These arguments will be passed to \code{\link[graphics]{axis}}. Internally, the arguments \code{side}, \code{at} and \code{labels} will be replaced in the list.
#'
#' @param pch1,pch2 These are the symbols to use in the plot for the classes 1, 2 and 3 in the \code{class} argument. \code{pch2} will overplot the \code{pch1}. These arguments will be passed to the \code{\link[graphics]{points}} and \code{\link[graphics]{legend}} functions. \code{pch1} and \code{pch2} must be of length 3.
#'
#' @param col1,col2 These are the colors of the symbols ploted (the underlying - col1 - and the overploted - col2 - symbol). Each one should be of length 3. These arguments will be passed to the \code{col} argument in the \code{\link[graphics]{points}} and \code{\link[graphics]{legend}} functions.
#'
#'@param pt.cex1,pt.cex2 These are numerical values giving the amount by which plotting symbols should be magnified relative to the default. They represent the undelying  - pt.cex1 - and the overploted - pt.cex2 - symbols. Each one should be of length 2. These arguments will be passed to the \code{cex} argument in the \code{\link[graphics]{points}} and \code{pt.cex} in \code{\link[graphics]{legend}} functions.
#'
#' @param top.lab.arg A list of arguments to plot the columns labels at the top of the plot. This list of arguments will be passed to \code{\link[graphics]{text}}. By default, variables labels stored as attribute of the \code{data.frame} are expected. If there is no labels there, and no labels are provided in this list, \code{NULL} will be passed to the \code{\link[graphics]{text}} will treat it accordingly.
#'
#' @param auto.legend Logical. Default is \code{TRUE}. If \code{FALSE}, no legend is ploted.
#'
#' @param legend.arg A list of arguments to plot the legend. This list of arguments will be passed to \code{\link[graphics]{legend}}. Internally, the arguments \code{legend}, \code{pch}, \code{col} and \code{pt.cex} will be replaced.
#'
#' @references
#'
#' Whiting PF, Rutjes AWS, Westwood ME, Mallett S, Deeks JJ, Reitsma JB, et al. QUADAS-2: A Revised Tool for the Quality Assessment of Diagnostic Accuracy Studies. Annals of Internal Medicine. 18 de outubro de 2011;155(8):529-36.
#'
#' Higgins JPT, Green S. Cochrane Handbook for Systematic Reviews of Interventions. 1 edition. Chichester, England; Hoboken, NJ: Wiley; 2008. 672 p.
#'
#' @examples
#' # Simulating a dataset
#' set.seed(12345)
#' Study = paste0(sample(c("Goodman", "White", "Fring", "Ehrmantraut"), 10, TRUE),", ", sample(LETTERS, 10, TRUE),"-", sample(1998:2015, 10, TRUE))
#' mydata <- data.frame(StudyID = Study,
#'                v1 = sample(c("Low risk", "Unclear risk", "High risk"), 10, TRUE),
#'                v2 = sample(c("Low risk", "Unclear risk", "High risk"), 10, TRUE),
#'                v3 = sample(c("Low risk", "Unclear risk", "High risk"), 10, TRUE),
#'                v4 = sample(c("Low risk", "Unclear risk", "High risk"), 10, TRUE),
#'                v5 = sample(c("Low risk", "Unclear risk", "High risk"), 10, TRUE),
#'                v6 = sample(c("Low risk", "Unclear risk", "High risk"), 10, TRUE),
#'                v7 = sample(c("Low risk", "Unclear risk", "High risk"), 10, TRUE))
#'
#' # Setting variables labels
#' attr(mydata, "var.labels") <- c("Study unique ID", "Patient selection",
#'                              "Index test", "Reference standard",
#'                              "Flow and timing", "Patient selection-aplicability",
#'                              "Index test-aplicability", "Reference standard-aplicability")
#'
#' # Making a barplot with the fractions of studies with each risk of bias
#' # Setting the bars collors
#' ramp <- colorRamp(c("darkgreen", "white"))
#'
#' # Making a table with the classification proportions
#' tabtmp <- sapply(2:5, function(i) prop.table(table(mydata[, i])))
#' colnames(tabtmp) <- c("Patient selection", "Index test", "Reference standard", "Flow and timing")
#'
#' # Warning! Resizing the window or setting different margins may be required
#' # Setting the plot with two windows
#' par(mar = c(5.1, 13.1, 2.1, 5.1), mfrow = c(2, 1))
#'
#' # Calling the superior plot
#' barplot(tabtmp, horiz = TRUE, cex.names = .9,
#'   xlab = 'Proportion of studies with low, high \n or unclear risk of bias',
#'   cex.lab = .9, las = 1, legend.text = FALSE, col = rgb(ramp(seq(0, 1, length = 3)),
#'   max = 255))
#'
#' # Making a table with the classification proportions for the inferior window
#' tabtmp <- sapply(6:8, function(i) prop.table(table(mydata[, i])))
#' colnames(tabtmp) <- c("Patient selection - applicability",
#'         "Index test  - applicability", "Reference standard - applicability")
#'
#' # Calling the inferior window plot
#' barplot(tabtmp, legend.text = rownames(tabtmp), horiz = TRUE, cex.names = .9,
#'   args.legend = list(x = 'top', inset = -.4, xpd = NA, horiz = TRUE, cex = .8),
#'   xlab = 'Proportion of studies with low, high, or unclear \n concerns regarding applicability',
#'   cex.lab = .9, las = 1, space = 1,
#'   col = rgb(ramp(seq(0, 1, length = 3)), max = 255))
#'
#' # Calling the QA.plot with default arguments
#' par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1, 1))
#' QA.plot(mydata)
#'
#' # The same plot but with squares, different colors, and fonts.
#' # Changing the 'class' argument order to match with the bars colors.
#' QA.plot(mydata, class = c("High risk", "Low risk", "Unclear risk"),
#'      pch2 = c("x","+","?"), pch1 = rep(15, 3),
#'      col1 = rgb(ramp(seq(0, 1, length = 3)), max = 255),
#'      study.arg = list(tick = FALSE, cex = .9, las = 1, font = 3, col = "lightgray"))
#'
#' rm(mydata, Study, tabtmp, ramp)
#'
#' @export
QA.plot <- function(tab, class = c("Low risk", "High risk", "Unclear risk"),
                    mar = c(3,10.1,7.1,8.1),
                    study.arg = list(tick = FALSE, cex = .9, las = 1),
                    pch1 = rep(16, 3), pch2 = c("+","x","?"),
                    col1 = c("green", "red", "orange"), col2 = rep("black", 3),
                    pt.cex1 = rep(2.5, 3), pt.cex2 = rep(.85, 3),
                    top.lab.arg = list(labels = attr(tab, "var.labels")[2:ncol(tab)], srt = 45, adj = 0, cex = .75, xpd = NA),
                    auto.legend = TRUE,
                    legend.arg = list(x = "bottom", bty = "n", inset = -0.10, horiz = TRUE, xpd = NA)) {
  # Conditions section
  if (!is.data.frame(tab) && !is.matrix(tab) && !is.table(tab)) {
    stop("'tab' is not a data.frame, table or matrix.")
  }
  if (length(class) != 3) {
    stop("'class' must be of length 3.")
  }
  if (!is.character(class)) {
    stop("'class' must be a character vector.")
  }
  if (length(class) != 3) {
    stop("'class' must be of length 3.")
  }
  if (length(mar) != 4) {
    stop("'mar' must be of length 4.")
  }
  if (!is.numeric(mar)) {
    stop("'mar' is not numeric.")
  }
  if (length(pch1) != 3) {
    stop("'pch1' must be of length 3.")
  }
  if (length(pch2) != 3) {
    stop("'pch2' must be of length 3.")
  }
  if (length(col1) != 3) {
    stop("'col1' must be of length 3.")
  }
  if (length(col2) != 3) {
    stop("'col2' must be of length 3.")
  }
  if (length(pt.cex1) != 3) {
    stop("'pt.cex1' must be of length 3.")
  }
  if (length(pt.cex1) != 3) {
    stop("'pt.cex1' must be of length 3.")
  }
  if (!is.numeric(pt.cex1)) {
    stop("'pt.cex1' must be numeric.")
  }
  if (!is.numeric(pt.cex2)) {
    stop("'pt.cex2' must be numeric.")
  }
  if (!is.logical(auto.legend)) {
    stop("'auto.legend' must be either TRUE or FALSE.")
  }
  if (any(is.null(top.lab.arg$labels))) {
    stop("At least one of 'labels' in the 'top.lab.arg' list is NULL.")
  }
  # Holdingn default graphical parameter
  opar <- par(mar = par()$mar)
  on.exit(par(opar))

  # Settingn an empty plot
  par(mar = mar)
  plot(x = c(1, (ncol(tab) - 1)), y = c(1, nrow(tab)), type = "n", axes = FALSE, xpd = NA, xlab = "", ylab = "")

  # Making an legend
  if (auto.legend) {
    legend.arg$legend <- class
    legend.arg$pch <- pch1
    legend.arg$col <- col1
    legend.arg$pt.cex <- pt.cex1
    do.call("legend", legend.arg)

    legend.arg$pch <- pch2
    legend.arg$col <- col2
    legend.arg$pt.cex <- pt.cex2
    legend.arg$pt.cex <- pt.cex2
    legend.arg$text.col <- "transparent"
    do.call("legend", legend.arg)
  }

  # Settingn a text for the columns (risk of bias dimensionsn names)
  top.lab.arg$x <- 1:(ncol(tab)-1)
  top.lab.arg$y <- nrow(tab) + 1
  do.call(text, top.lab.arg)

  # Ploting the studies IDS at the left axis
  study.arg$side <- 2
  study.arg$at <- 1:nrow(tab)
  study.arg$labels <- tab[, 1]
  do.call(axis, study.arg)
  # axis(4, 1:nrow(tab), tab[,9], F, hadj = adj.study, pos = pos.study, cex.axis = cex.study, col = col.study, font = font.study, line=study.line, las = 1)

  # PLoting the symbols with the risk assessment
  for (i in 1:nrow(tab)) {
    for (j in 2:ncol(tab)) {
      if (tab[i, j] == class[1]) {
        points(j - 1, i, type = "p", pch = pch1[1], col = col1[1], cex = pt.cex1[1], xpd = NA)
        points(j - 1, i, type = "p", pch = pch2[1], col = col2[1], cex = pt.cex2[1], xpd = NA)
      }
      if (tab[i, j] == class[2]) {
        points(j - 1, i, type = "p", pch = pch1[2], col = col1[2], cex = pt.cex1[2], xpd = NA)
        points(j - 1, i, type = "p", pch = pch2[2], col = col2[2], cex = pt.cex2[2], xpd = NA)
      }
      if (tab[i, j] == class[3]){
        points(j - 1, i, type = "p", pch = pch1[3], col = col1[3], cex = pt.cex1[3], xpd = NA)
        points(j - 1, i, type = "p", pch = pch2[3], col = col2[3], cex = pt.cex2[3], xpd = NA)
      }
    }
  }
}
