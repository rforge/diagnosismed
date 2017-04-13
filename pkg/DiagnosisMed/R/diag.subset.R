#' @name display
#'
#' @title Diagnostic tests display with forest plots.
#'
#' @description This group of functions have the purpose to graphically display and compare either the same diagnostic test across several subsets (e.g. groups formed by clinical characteristics), several diagnostic tests or a combination of both. The \code{\link{diagnosis}} function is invocked and tests perfomances are stacked in a \code{data.frame}, and forest plots may be called from these.
#'
#' \code{diag.subset} and \code{forest.diag.subset} display test performance across one or several groups represented by an additional columns in the same \code{data.frame}.
#'
#' \code{diag.stack} and \code{forest.diag.stack} display two or more tests performances against the very same reference standard in the same sample. The index tests are represented by an additional columns in the same \code{data.frame}.
#'
#' \code{diag.by} and \code{forest.diag.by} display two or more tests performances against the very same reference standard, across subsets defined by a single variable. The index tests and the subset variable are represented by an additional columns in the same \code{data.frame}.
#'
#' @param ref A character vector of length 1 representing the reference standard variable in data. If it is formated as character, the function will convert it to a factor, then use the reference factor level as the absence of the condition. See \code{\link{diagnosis}}
#'
#' @param test A character vector representing the test under evaluation variable in data. If the test variables are formated as character, the function will convert it to a factor, then use the reference factor level as the absence of the condition. See \code{\link{diagnosis}}. In \code{diag.subset} this must be o length 1, in \code{diag.by} and \code{diag.stack} this can be of length 2 or more.
#'
#' @param group.var A character vector representing the names of categorical variables in data to make subsets of data. In \code{diag.by}, this argument must be a single variable, while in \code{diag.subset} it could two or more variables.
#'
#' @param data A \code{data.frame} contining the required variables
#'
#' @param CL Confidece interval limit. Default is 0.95. Must be a number between 0 and 1. See \code{\link{binom.CI}}
#'
#' @param CL.type Method to estimate the confidece interval. Default is "wilson". See \code{\link{binom.CI}}
#'
#' @param var.labels For \code{diag.subset} and \code{diag.stack}, this is the labels (in the same order) of the variables in the \code{group.var} and \code{test} arguments, respectively. If not informed, the function will try to retrieve the labels from the "var.labels" attribute of the original data.frame. If there is no labels in this attribute, the function will return the variable names.
#'
#' @param group.var.labels,test.var.labels For \code{diag.by}, these are the labels of the variables in the \code{group.var} and \code{test} arguments (respectively), respectively. If not informed, the function will try to retrieve the labels from the "var.labels" attribute of the original data.frame. If there is no labels in this attribute, the function will return the variable names.
#'
#' @param x An object resulting from the function \code{diag.subset}, \code{diag.stack} or \code{diag.by}.
#'
#' @param type Type of forest plot. \code{"SeSp"} is a three window plot, showing Sensitivity and Specificty in the right windows. \code{"PLR"} (not yet implemented) is a two window plot showing the Positive likelihood ratio in right window.
#'
#' @param mar1,mar.se,mar.sp,mar.plr These are the margins that will be passed to \code{par()$mar} in the left window, the sensitivity window, the specificity window and the positive likelyhood ratio window respectively. See \code{\link[graphics]{par}}.
#'
#' @param seg.list A list of arguments that will be passed to \code{\link[graphics]{segments}} to plot the sensitivity, specificity and positive likelihood ratio confidence intervals. Internally the x and y coordenates are automatically replaced.
#'
#' @param points.list A list of arguments that will be passed to \code{\link[graphics]{points}} to plot the sensitivity, specificity and positive likelihood ratio point estimate. Internally the x and y coordenates are automatically replaced.
#'
#' @param num.list A list of arguments that will be passed to \code{\link[graphics]{text}} to plot in the left window the samples N and Prevalences. Internally the y coordenates are automatically replaced.
#'
#' @param var.list A list of arguments that will be passed to \code{\link[graphics]{text}} to plot variables or major group names (or labels) in the left window. Internally the y coordenates are automatically replaced.
#'
#' @param cat.list A list of arguments that will be passed to \code{\link[graphics]{text}} to plot category names (or labels) in the left window . Internally the y coordenates are automatically replaced.
#'
#' @param x.n,x.prev x axis coordenates of the samples N and Prevalences in the left window (from 0 to 1). Internally, these will replace the x argument in the \code{cat.list}.
#'
#' @param adj.se,adj.sp,adj.plr The adjustment of the labels argument passed to \code{\link[graphics]{text}}, to plot the sensitivity, specificity and positive likelihood ratio (respectively). 0 means left adjustment and 1 means right adjustment.
#'
#'@param se.xlim,sp.xlim,plr.xlim Limits of the horizontal axis of each window. This will be pased to \code{\link[graphics]{plot.default}}. This may receive either \code{"auto"} or a numeric vector of length 2, representing the lower and upper limits of the horizontal axis. If auto, the function picks the lower values of the lower confidence limit and the upper value of the upper confidence limits.

#' @param se.xlab,sp.xlab,plr.xlab Labels of the horizontal axis of each window. This will be pased to \code{\link[graphics]{plot.default}}
#'
#' @param se.pos,sp.pos,plr.pos This is added to the horizontal positioning to the left of the lower limit of the axis (in case of sensitivity and positive lieklihood ratio) to plot the text of the point and confidence interval estimates. In the case of specificity, this is added to the right of th upper limit of the horizontal axis. This is intended to be an argument to make the text scoot over from the plot area and avoid overplot with the confidence interval lines in the edge of the graph.
#'
#' @param se.col.label,sp.col.label,plr.col.label Character strings representing the labels of the sensitivity, specificity and positive likelihood ratios columns.
#'
#' @param digits Number of digits to plot all numbers. Should an integer number.
#'
#' @param grid Logical. If \code{TRUE}, the default value, it invokes the \code{\link[graphics]{grid}} function with its default arguments.
#'
#' @return \code{diag.by}, \code{diag.stack}, \code{diag.subset} return data.frames with each strata sample size, sample prevalence, the sensitivities, the specificities and positive likelihood ratios (and their confidence intervals). The \code{forest.diag.by}, \code{forest.diag.stack}, \code{forest.diag.subset} return plots.
#'
#' @seealso \code{\link{diagnosis}}, \code{\link{LRgraph}}
#'
#' @export
diag.subset <- function(ref, test, group.var, data, CL = .95, CL.type = c("wilson", "exact", "approximate"), var.labels = attr(data, "var.labels")[match(group.var, names(data))]){
  if (!is.character(ref)) {
    stop("'ref' should be a character vector.")
  }
  if (length(ref) != 1) {
    stop("'ref' should have length 1.")
  }
  if (!is.character(test)) {
    stop("'test' should be a character vector.")
  }
  if (length(test) != 1) {
    stop("'test' should have length 1.")
  }
  if (!is.character(group.var)) {
    stop("'group.var' should be a character vector.")
  }
  if (!is.data.frame(data)) {
    stop("'data' should be a data.frame.")
  }
  if (!all(group.var %in% names(data))) {
    stop(paste0("The following variables of 'group.var' are not in 'data': ",group.var[!(group.var %in% names(data))], collapse = ","))
  }
  if (!is.null(var.labels)) {
    if (length(var.labels) != length(group.var)) {
      warning("'group.var' and 'var.labels' do not have the same lenth. Check carefully if labels match correctly.")
    }
  }
  Variables <- c("Overall", unlist(sapply(1:length(group.var),
                function(i) rep(group.var[i], nlevels(data[, group.var[i]])))))
  Levels <- c(NA, unname(unlist(sapply(data[, group.var], levels))))
  output <- data.frame(Variables, Levels, n = NA, prevalence = NA, Sensitivity = NA,
            Se.inf.cl = NA, Se.sup.cl = NA, Specificity = NA, Sp.inf.cl = NA, Sp.sup.cl = NA,
            PLR = NA, PLR.inf.cl = NA, PLR.sup.cl = NA)
  output[which(output$Variables == "Overall"), 3:13] <- diagnosis(tab = table(data[, test], data[, ref]), CL = CL, CL.type = CL.type[1])[c(2:3,6:14)]
  # i = 7
  for (i in 2:nrow(output)) {
    cond <- which(data[, as.character(output[i, "Variables"])] == as.character(output[i, 2]))
    tab <- table(data[cond, test], data[cond, ref])
    if (all(dim(tab) == 2)) {
      output[i, 3:13] <- diagnosis(tab = tab, CL = CL, CL.type = CL.type[1])[c(2:3,6:14)]
    } else {
      warning(paste0("Not able to make a 2x2 table in the subset with the variable '", as.character(output[i, "Variables"]),"' in level ", as.character(output[i, "Levels"])))
    }
  }
  output$Variables <- as.character(output$Variables)
  for (i in length(output$Variables):2) {
    if (output$Variables[i-1] == output$Variables[i]) {output$Variables[i] <- NA}
  }
  for (i in 2:nrow(output)) {
    if (all(is.na(output[i,3:13]))) {output <- output[-i,]}
  }
  if (!is.null(var.labels))
    for (i in 2:length(output$Variables)) {
      output$Variables[which(!is.na(output$Variables))][-1] <- var.labels
    }
  class(output) <- c("data.frame", "diag.subset")
  output
}

