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