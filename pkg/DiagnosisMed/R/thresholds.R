#' Collection of functions to estimate decision thresholds for diagnostic tests.
#' @name thresholds
#'
#' @description This collection is intended to be a intermediate function, not be used by the end user. One may wish to use \code{\link{TGROC}} instead, as it calls \code{thresholds} and other functions at once for a more complete analysis including a flexible plot function.
#'
#' By a variety of methods, this collection may either dichotomize or threechotomize a diagnostic test scale, finding the thresholds to classify subjects with and without a condition. Along with the threshold, these functions return the sensitivity, specificity and the positive likelihood ratios (and their respective confidence interval) if the test is dichotomized at that threshold. The \code{thresholds} function is a wrapper for all methods at once. The methods currently available are (see details for the formulas):
#' \itemize{
#'   \item \code{Se.equals.Sp} The threshold which Sensitivity is equal to Specificity.
#'   \item \code{Max.Accuracy} The threshold which maximize the accuracy.
#'   \item \code{Max.DOR} The threshold which maximize the diagnostic odds ratio.
#'   \item \code{Min.Error} The threshold which minimizes the error rate.
#'   \item \code{Max.Accuracy.area} The threshold which maximize the accuracy area.
#'   \item \code{Max.Youden} The threshold which maximize the Youden J index.
#'   \item \code{Min.ROC.dist} The threshold which minimize the distance between the curve and the upper left corner of the graph.
#'   \item \code{Max.Efficiency} The threshold which maximize the efficiency.
#'   \item \code{Min.MCT} The threshold which minimize the misclassification cost term.
#'   \item \code{thresholds} Wrapper for all the above methods.
#'   \item \code{inc.limits} Threechotomizes the test according to a minimum required sensitivity and specificity.
#' }
#'
#' @param x is the output from any of the \code{\link{SS}}, \code{\link{BN.SS}} or \code{\link{NN.SS}} functions.
#'
#' @param Inconclusive This is a value that ranges from 0 to 1 that will identify the test range where the performance of the test is not acceptable and thus considered inconclusive. It represents the researcher tolerance of how good the test should be. If it is set to 0.95 (which is the default value), test results that have less than 0.95 sensitivity and specificity will be in the inconclusive range. Also known as the minimum required sensitivity and specificity.
#'
#' @param Cost Cost = cost(FN)/cost(FP). Use in the \code{Min.MCT} (minimizing the misclassification cost term) function. It is a value in a range from 0 to infinite. Could be financial cost or a health outcome with the perception that FN are more undesirable than FP (or the other way around). Cost = 1 means FN and FP have even cost. Cost = 0.9 means FP are 10 percent more costly. Cost = 0.769 means that FP are 30 percent more costly. Cost = 0.555 means that FP are 80 percent more costly. Cost = 0.3 means that FP are 3 times more costly. Cost = 0.2 means that FP are 5 times more costly. Also, it can be more easily inserted as any ratio such as 1/2.5 or 1/4.
#'
#' @param pop.prevalence The disease prevalence informed by the user. If not informed, it will be the same as the sample prevalence. This will be passed to \code{Max.Efficiency} and \code{Min.MCT}. Particularly interesting if the test will be applied to a population with a different condition prevalence.
#'
#' @details Occasionally the dichotomizing methods may find ties, i.e. more than one threshold matching the criteria. This may occur particularly in small sample sizes or with the non-parametric analysis. If this is the case, the functions will return a warning and automatically pick the median value, or the value closest to the median.
#'
#' Similar phenomena may occur with the \code{inc.limits} function. If this is the case, for specificity, the function will pick the lowest threshold, and for sensitivity, the function will pick the highest threshold. Additionally, one may notice that frequently the sensitivity and specificity output does not exactly match the Inconclusive argument, e.g. 0.90. That depends on the data, that may have small or big jumps of sensitivity and specificity from one threshold to another. It seems this phenomena is more frequent with small sample sizes, and with the non-parametric analysis. If this is the case, the \code{inc.limits} will always pick the threshold with the sensitivity and specificity nearest to the Inconclusive argument.
#'
#' The formulas for the methods used are:
#' \itemize{
#'   \item \code{Max.Accuracy} \eqn{(TN+TP)/sample size}
#'   \item \code{Max.DOR} \eqn{(TP*TN)/(FN*FP)}
#'   \item \code{Min.Error} \eqn{(FN+FP)/sample size}
#'   \item \code{Max.Accuracy.area} \eqn{(TP*TN)/((TP+FN)*(FP+TN))}
#'   \item \code{Max.Youden} \eqn{Se+Sp-1}
#'   \item \code{Min.ROC.dist} \eqn{(Sp-1)^2+(1-Se)^2}
#'   \item \code{Max.Efficiency} \eqn{Se*prevalence+(1-prevalence)*Sp}
#'   \item \code{Min.MCT} \eqn{(1-prevalence)*(1-Sp)+Cost*prevalence(1-Se)}
#' }
#'
#' @return
#'A data.frame with the threshold (test.value), sensitivity, specificity and positive likelihood ratio (with their respective confidence interval). The row names will match the methods.
#'
#' @seealso \code{\link{SS}}, \code{\link{TGROC}}
#'
#' @examples
#' data(rocdata)
#'
#' # Thresholds and inconclusive limits from non-parametric ROC analysis
#' x <- SS(ref = rocdata$Gold, test = rocdata$test2)
#' thresholds(x)
#' inc.limits(x)
#'
#' # Thresholds and inconclusive limits from smoothed analysis with NN
#' y <- NN.SS(x)
#' thresholds(y)
#' inc.limits(y)
#'
#' data("tutorial")
#'
#' # Thresholds and inconclusive limits at a very accurate test.
#' x <- SS(ref = ifelse(tutorial$Gold == "pos", 1, 0), test = tutorial$Test_A)
#' # Notice that ties occurs.
#' thresholds(x)
#' # Notice that the "Lower  inconclusive" test.values is higher than the "Upper inconclusive"
#' # This may occur if the test is very accurate, or argument "Inconclusive" is too low.
#' # In this case a inconclusive range may not be applicable.
#' inc.limits(x)
#' # But, for this data, increasing the "Inconclusive" solves the issue.
#' inc.limits(x, Inconclusive = .99)
#'
#' # When smoothing the analysis, the ties do not occur.
#' y <- NN.SS(x)
#' thresholds(y)
#' # Same issue with the parametric analysis.
#' # But it returns very different threshold values.
#' inc.limits(y, Inconclusive = .99)
#'
#' rm(rocdata, tutorial, x, y)
#'
#' @export
# Se=Sp threshold where x is a SS object----------------------------------------
Se.equals.Sp <- function(x){
  x$table$Se.equals.Sp <- abs(x$table$Specificity - x$table$Sensitivity)
  condition <- which(x$table$Se.equals.Sp == min(x$table$Se.equals.Sp))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning("Se = Sp is reached at the following test values: ", toString(warn),". \n  The one nearest to median was chosen!")
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Se = Sp"
  output
}
# Se.equals.Sp(x) ; Se.equals.Sp(NN.SeSp)

#' @rdname thresholds
#' @export
# Maximizing the accuracy threshold where x is a SS object----------------------
Max.Accuracy <- function(x){
  x$table$Accuracy <- (x$table$TN + x$table$TP) / x$sample.size
  condition <- which(x$table$Accuracy == max(x$table$Accuracy))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Accuracy reaches its maximum at the following test values: ",toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max Accuracy"
  output
}
# Max.Accuracy(x) ; Max.Accuracy(NN.SeSp)

#' @rdname thresholds
#' @export
# Maximizing Diagnostic Odds Ratio threshold where x is a SS object------------
Max.DOR <- function(x){
  x$table$DOR <- ((x$table$TN)*(x$table$TP))/((x$table$FP)*(x$table$FN))
  x$table$DOR <- ifelse(x$table$DOR == Inf, NA, x$table$DOR)
  condition <- which(x$table$DOR == max(x$table$DOR, na.rm = T))
  if(length(condition) > 1){
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("DOR reaches its maximum at the following test values: ",toString(warn),".  \n The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max DOR"
  output
}
# Max.DOR(x)

#' @rdname thresholds
#' @export
# Minimizing error rate threshold where x is a SS object------------------------
Min.Error <- function(x){
  x$table$Error.rate <- ((x$table$FP) + (x$table$FN)) / x$sample.size
  condition <- which(x$table$Error.rate == min(x$table$Error.rate))
  if  (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Error rate reaches its minimum at the following test values: ",toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Min Error rate"
  output
}
# Min.Error(x)

#' @rdname thresholds
#' @export
# Maximizing the accuracy area threshold where x is a SS object------------------
Max.Accuracy.area <- function(x){
  # D and ND are constants will not make any difference in the final result
  # removing them will make it easier for smoothed SS objects
  # D <- sum(x$table$D)
  # ND <- sum(x$table$ND)
  x$table$Accuracy.area <- ((x$table$TP)*(x$table$TN)) # / (D * ND)
  condition <- which(x$table$Accuracy.area == max(x$table$Accuracy.area))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Accuracy area reaches its maximum at the following test values: ",
            toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl",
                           "Specificity","Sp.inf.cl","Sp.sup.cl","PLR",
                           "PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max Accuracy area"
  output
}
# Max.Accuracy.area(x); Max.Accuracy.area(NN.SeSp)

#' @rdname thresholds
#' @export
# Maximizing the Youden J index threshold where x is a SS object----------------
Max.Youden <- function(x){
  x$table$Youden <- x$table$Sensitivity + x$table$Specificity - 1
  condition <- which(x$table$Youden == max(x$table$Youden))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Youden J index reaches its maximum at the following test values: ",
            toString(warn),".\n  The one closest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl",
                           "Specificity","Sp.inf.cl","Sp.sup.cl","PLR",
                           "PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max Youden"
  output
}
# Max.Youden(x)

#' @rdname thresholds
#' @export
# Minimizing The ROC 0 1 distance threshold where x is a SS object-------------
 Min.ROCdist <- function(x){
  x$table$MinRocDist <- (x$table$Specificity - 1) ^ 2 +
                        (1 - x$table$Sensitivity) ^ 2
  condition <- which(x$table$MinRocDist == min(x$table$MinRocDist))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Minimum ROC distance reaches its minimum at the following test values: ",
                   toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl",
                           "Specificity","Sp.inf.cl","Sp.sup.cl","PLR",
                           "PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Min ROC distance"
  output
}
# Min.ROCdist(x)

 #' @rdname thresholds
 #' @export
# maximizing Efficiency threshold where x is a SS object------------------------
Max.Efficiency <- function(x, pop.prevalence = NULL){
  if (is.null(pop.prevalence)) {
    pop.prevalence <- x$sample.prevalence
  }
  x$table$Efficiency <- (x$table$Sensitivity * (pop.prevalence)) + ((1 -
                        (pop.prevalence)) * x$table$Specificity)
  condition <- which(x$table$Efficiency == max(x$table$Efficiency))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Maximum Efficiency reaches its maximum at the following test values: ",
                   toString(warn), ".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl",
                           "Specificity","Sp.inf.cl","Sp.sup.cl","PLR",
                           "PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max Efficiency"
  output
}
# Max.Efficiency(x)

#' @rdname thresholds
#' @export
# Minimizing MissClassificatio Cost Term threshold where x is a SS object-------
Min.MCT <- function(x, pop.prevalence = NULL, Cost = 1){
  if (is.null(pop.prevalence)) {
    pop.prevalence <- x$sample.prevalence
  }
  x$table$MCT <- (1 - (pop.prevalence)) * (1 - x$table$Specificity) + (Cost *
                  (pop.prevalence)) * (1 - x$table$Sensitivity)
  condition <- which(x$table$MCT == min(x$table$MCT))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Minimum Misclassification cost term reaches its minimum at
                   the following test values: ",toString(warn),
                  ".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl",
                           "Specificity","Sp.inf.cl","Sp.sup.cl","PLR",
                           "PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Min MCT"
  output
}
# Min.MCT(x)

#' @rdname thresholds
#' @export
# Wraping all threshold functions where x is a SS object------------------------
thresholds <- function(x, pop.prevalence = NULL, Cost = 1){
  output <- rbind(Max.Youden(x),
                  Max.Accuracy(x),
                  Max.Accuracy.area(x),
                  Max.DOR(x),
                  Min.Error(x),
                  Min.ROCdist(x),
                  Se.equals.Sp(x),
                  Min.MCT(x, pop.prevalence = pop.prevalence, Cost = Cost),
                  Max.Efficiency(x, pop.prevalence = pop.prevalence)
                  )
  output
}
# thresholds(x, Cost = 2) ; thresholds(NN.SeSp)

#' @rdname thresholds
#' @export
# Inconclusive thresholds (threechotomization) ---------------------------------
inc.limits <- function(x, Inconclusive = .95){
  # Checking the Se values
  condition <- which.min(abs(Inconclusive - x$table$Sensitivity))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    warning(paste0("Sensitivity matches the minimum required at the following test values: ",toString(warn),". Highest one was chosen!"))
  }
  Se.pos <- condition[length(condition)]

  # Checking the Sp values
  condition <- which.min(abs(Inconclusive - x$table$Specificity))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    warning(paste0("Specificity matches the minimum required at the following
                   test values: ",toString(warn),". First one was chosen!"))
  }
  Sp.pos <- condition[1]
  # Extracting the test.values and corresponding validity
  output <- rbind(
            x$table[Se.pos,c("test.values","Sensitivity","Se.inf.cl",
                             "Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl",
                             "PLR","PLR.inf.cl","PLR.sup.cl")],
            x$table[Sp.pos,c("test.values","Sensitivity","Se.inf.cl",
                             "Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl",
                             "PLR","PLR.inf.cl","PLR.sup.cl")]
            )
  rownames(output) <- c("Lower inconclusive","Upper inconclusive")
  output
}
# inc.limits(x)
