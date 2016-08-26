#' Area under the receiver operating charachteristic (ROC) curve
#'
#' \code{np.auROCc} estimates the non-parametric (DeLong method) area under the ROC curve,  its standard error and confidence interval. It confidence interval has the purpose to test if the AUC is different from 0.5. It is not appropriate to use this confidence interval to compare with other auc due to lack of adjustment of possible correlation.
#'
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param CL Confidence limit. The limits of the confidence interval. Must be coded as number in range from 0 to 1. Default value is 0.95
#' @param reverse "auto" (deafult), TRUE or FALSE are acceptable values. ROC analysis assumes that higher values of the test are from subjects with the condition and lower values are from subjects without the condition. If it occurs the other way around, the ROC analysis and its interpretation must be reversed. If "auto", \code{np.auROCc} internally checks if the auc is below 0.5. If this is the case, it returns a warning and makes the reversion as 1 - auc, then conducts the standard error and confidence interval estimation.
#'
#' @seealso TGROC
#'
#' @export
np.auROCc <- function(ref, test, CL = 0.95, reverse = "auto") {
  # Warning section ...
  if (any(is.na(test) | is.na(ref))) {
    stop('It seems there are NAs either in the index test or in the reference test. Consider imputing or removing NAs!')
  }
  if (any(levels(as.factor(ref)) != c(0,1))) {
    stop("Your reference standard must be coded as 0 (absence) and 1 (presence). Check reference categories!")
  }
  if (reverse != "auto" && !is.logical(reverse)) {
    stop("reverse must be either 'auto', TRUE or FALSE.")
  }
  if (!is.numeric(CL) || CL < 0 || CL > 1) {
    stop("Confidence limit (CL) must be nemeric between 0 and 1.")
  }

  # Estimating the AUC
  X <- sort(test[which(ref == 0)])
  Y <- sort(test[which(ref == 1)])
  m <- length(X)
  n <- length(Y)
  AUC <- (m * n + (m * (m + 1)) / 2 - sum(rank(test, ties.method = "average")[which(ref == 0)])) / (m * n)

  # Reversing the AUC
  if (reverse == "auto") {
    if (AUC < .5) {
      reverse <- TRUE
    } else {
      reverse <- FALSE
    }
  }
  if (reverse) {
    AUC <- 1 - AUC
    warning("The area under the ROC curve was reversed!")
  }

  # Estimatingn the AUC SE and IC
  D10X <- function(Xi) { (1 / n) * sum(Y >= Xi[1]) }
  D01Y <- function(Yi) { (1 /m) * sum(Yi[1] >= X) }
  VAR.AUC <- sum((tapply(X, X, "D10X") - AUC) ^ 2) / (m * (m - 1)) + sum((tapply(Y, Y ,"D01Y") - AUC) ^ 2) / (n * (n - 1))
  SE.AUC <- sqrt(VAR.AUC)
  alpha <- 1 - CL
  output <- c(AUC, SE.AUC, AUC - qnorm(1 - alpha / 2) * SE.AUC, AUC + qnorm(1 - alpha / 2) * SE.AUC)
  names(output) <- c("AUC", "AUC.SE", "AUC.lower.CL", "AUC.upper.CL")
  output
}