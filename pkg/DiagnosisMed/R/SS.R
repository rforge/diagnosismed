#' ROC analysis - Sensitivity and Specificity trade-off
#' @name SS
#'
#' @description \code{SS} collection is intended to be intermediate functions, not be use by the end user. One may wish to use \code{\link{TGROC}} instead, as it calls \code{SS}, \code{BN.SS} and \code{NN.SS} and other functions at once for a more complete analysis including a flexible plot function.
#'
#'  \code{SS}, \code{BN.SS} and \code{NN.SS} compute validity measures for each decision threshold of a continuous scale diagnostic test with their respective confidence intervals. It shows the trade-off of the Sensitivity and Specificity values with progressive changes of the threshold. \code{SS} does it non-parametrically and uses \code{\link{binom.CI}} to estimate the confidence intervals. \code{NN.SS} does it by fitting a feed forward neural network with the \pkg{AMORE} package. The underling idea is that this analysis is a robust way to smooth the Sensitivity and Specificity trade-off and represent the population, as neural networks may approximate any population function distribution. One may notice that running the neural network more than once may retun slightly different values. This is expected as it depends on the fit for each run. \code{BN.SS} does the same thing but assuming that the test values from subjects with and without the condition have Gaussian distribution (bi-normally distributed test values). Both \code{BN.SS} and \code{NN.SS} use a Gaussian confidence interval estimation.
#'
#' @param  ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition).
#'
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#'
#' @param x For the \code{NN.SS} function, x is the output of the \code{SS} function.
#'
#' @param reverse "auto" (default), TRUE or FALSE are the acceptable values. ROC analysis assumes that higher values of the test are from subjects with the condition, and lower values are from subjects without the condition. If it occurs the other way around, the ROC analysis and its interpretation must be reversed. If "auto", \code{SS} internally checks if the mean (or median) test values are higher among subject without the condition. If this is the case, it returns a warning and sets the reverse = TRUE. The reversion is simply done by multiplying all test values by -1, make all the computations and returning the absolute values.
#'
#' @param CL Confidence limit. The limits of the confidence interval. Must be coded as number in range from 0 to 1. Default value is 0.95.
#'
#' @param binom.conf Method of binomial confidence interval. "wilson" (default), "exact" and "approximate" are acceptable. See \code{\link{binom.CI}}.
#'
#' @param pop.prevalence Population condition prevalence. If this values is not NULL, the sample prevalence is internally replaced by the population prevalence to estimate the confidence intervals and, for \code{NN.SS} and \code{BN.SS} to estimate the TP, TN, FP, and FN fractions needed in other functions to estimate decision \code{\link{thresholds}}. So, use it wisely. Particularly interesting with data from case-control design.
#'
#' @param t.min,t.max,precision Test minimum, maximum and intervals to simulate the parametric estimation as seq(t.min, t.max, precision). The \code{NN.SS} and \code{BN.SS} functions need test values to which simulate the Sensitivity and Specificity values. If left NULL, internally the function will pick the test maximum and minimum and the default precision. However, it does not need to be the same values observed in data. To give a smooth appearance and allow the \code{NN.SS} to fit nicely, the sequence should have hundreds of tests values to simulate, say at least 200 values. A good rationale to create this sequence is to simulate all possible values of the test results. If this rationale does not have 200 values, than, perhaps, every half value should do. At the another extreme, creating a sequence with too many values (e.g. 2000) may create unrealistic test values, increase computational time (or even explode memory) and may not increase the smoothness to the parametric analysis.
#'
#' @param n.neurons Numeric vector containing the number of neurons of each layer. See \code{\link[AMORE]{newff}}.
#'
#' @param learning.rate.global Learning rate at which every neuron is trained. See \code{\link[AMORE]{newff}}.
#' @param momentum.global Momentum for every neuron. See \code{\link[AMORE]{newff}}.
#'
#' @param error.criterium Criteria used to measure to proximity of the neural network prediction to its target. See \code{\link[AMORE]{newff}}.
#'
#' @param Stao Stao parameter for the TAO error criteria. See \code{\link[AMORE]{newff}}.
#'
#' @param hidden.layer Activation function of the hidden layer neurons. See \code{\link[AMORE]{newff}}.
#'
#' @param output.layer Activation function of the hidden layer neurons. See \code{\link[AMORE]{newff}}.
#'
#' @param method Preferred training method. See \code{\link[AMORE]{newff}}.
#'
#' @param report Logical value indicating whether the training function should keep quiet. See \code{\link[AMORE]{train}}.
#'
#' @param show.step Number of epochs to train non-stop until the training function is allow to report. See \code{\link[AMORE]{train}}.
#'
#' @param n.shows Number of times to report (if report is TRUE). See \code{\link[AMORE]{train}}.
#'
#' @details Tests results matching the cut-off values will be considered a positive test. \code{SS}, as all ROC analysis, assumes that subjects with higher values of the test are with the target condition, and those with lower values are without the target condition. Tests that behave like glucose (middle values are supposed to be normal and extreme values are supposed to be abnormal) will not be correctly analyzed. This may be the degenerated data for ROC analysis. If for a particular tests higher values are from subjects without the condition and lower values are from subject with the condition, the analysis must be reversed. In this case, \code{SS} does it by multiplying the test results by -1 before analysis, and returning its absolute value before output.
#'
#' @return
#' \code{table} A dataset with the number of subjects with and without the condition, the TP, TN, FP, and FP, Sensitivity, Specificity, predictive values and likelihood ratios for each threshold.
#'
#' \code{sample.size} The sample size.
#'
#' \code{sample.prevalence} The sample prevalence. However, if the argument pop.prevalence is not NULL, it returns the pop.prevalence despite not replacing the name.
#'
#' @seealso \code{\link{binom.CI}}, \code{\link{np.auROCc}}, \code{\link{thresholds}}, \code{\link{TGROC}}
#'
#' @examples
#' data("rocdata")
#' # Artificially forcing the reversion (this is NOT correct for this data)
#' x <- SS(ref = rocdata$Gold, test = rocdata$test1, reverse = TRUE)
#' # Printing just a subset of the table
#' # In this case sensitivity has higher values at higher test values
#' tail(x$table[,c("test.values","Sensitivity","Specificity")])
#' # And specificity has higher values at lower test values
#' head(x$table[,c("test.values","Sensitivity","Specificity")])
#'
#' # The same analysis without forcing the reversion
#' x <- SS(ref = rocdata$Gold, test = rocdata$test1)
#' # Printing just a subset of the table
#' # In this case sensitivity has higher values at lower test values
#' head(x$table[,c("test.values","Sensitivity","Specificity")])
#' # And specificity has higher values at higher test values
#' tail(x$table[,c("test.values","Sensitivity","Specificity")])
#'
#' # Smoothingn with bi-normal function
#' # It will be easier to check the fit graphically with TGROC
#' # Rejecting the assumption of normality for those without the condition.
#' shapiro.test(rocdata$test1[which(rocdata$Gold == 1)])
#' shapiro.test(rocdata$test1[which(rocdata$Gold == 0)])
#' z <- BN.SS(ref = rocdata$Gold, test = rocdata$test1, t.min = 0.005, t.max = 2, precision = 0.005)
#' head(z$table[,c("test.values","Sensitivity","Specificity")])
#' tail(z$table[,c("test.values","Sensitivity","Specificity")])
#'
#' # Smoothingn with neural network
#' y <- NN.SS(x, t.min = 0.005, t.max = 2, precision = 0.005)
#' head(y$table[,c("test.values","Sensitivity","Specificity")])
#' tail(y$table[,c("test.values","Sensitivity","Specificity")])
#'
#' rm(rocdata, x, y, z)
#' @export
# Non-parametric trade-off (ROC-analysis)--------------------------------------
SS <- function(ref, test, reverse = "auto", CL = 0.95, binom.conf = "wilson",
               pop.prevalence = NULL) {
  # Warning section ...
  if (any(is.na(test) | is.na(ref))) {
    stop("It seems there are NAs either in the index test or in the reference
         test. Consider imputing or removing NAs!")
  }
  if (any(levels(as.factor(ref)) != c(0, 1))) {
    stop("Your reference standard must be coded as 0 (absence) and 1 (presence).
         Check reference categories!")
  }
  if (reverse != "auto" && !is.logical(reverse)) {
    stop("reverse must be either 'auto', TRUE or FALSE.")
  }
  if (!is.numeric(CL) || CL < 0 || CL > 1) {
    stop("Confidence limit (CL) must be nemeric between 0 and 1.")
  }
  if (reverse == "auto") {
    if (mean(test[which(ref == 0)]) > mean(test[which(ref == 1)]) ||
        median(test[which(ref == 0)]) > median(test[which(ref == 1)])) {
      reverse <- TRUE
    } else {
      reverse <- FALSE
    }
  }
  if (reverse) {
    test <- test * -1
    warning("The ROC analysis was reversed!")
  }

  test.table <- table(test, ref)

  D <- length(test[which(ref == 1)])
  ND <- length(test[which(ref == 0)])
  sample.size <- ND + D
  if (!is.null(pop.prevalence)) {
    sample.prevalence <- pop.prevalence
  } else {
    sample.prevalence <- D / sample.size
  }

  # Taking the rownames of the test.table to be results first column
  test.values <- as.numeric(rownames(test.table))
  tab <- as.data.frame(test.values)
  tab$D <- test.table[, 2]
  tab$ND <- test.table[, 1]
  tab$TP <- sapply(1:nrow(test.table), function(i)sum(
    test.table[i:nrow(test.table), 2]))
  if (test.table[1, 2] == 0) {fFN <- 0} else {fFN <- 1}
  tab$FN <- c(as.integer(fFN),
          sapply(2:nrow(test.table),  function(i)sum(test.table[1:(i - 1), 2])))
  tab$FP <- sapply(1:nrow(test.table),
                              function(i)sum(test.table[i:nrow(test.table), 1]))
  if (test.table[1, 1] == 1) {fTN <- 0} else {fTN <- 1}
  tab$TN <- c(as.integer(fTN), sapply(2:nrow(test.table),
                                      function(i)sum(test.table[1:(i - 1), 1])))

  tmp.CI <- binom.CI(tab$TP, D, conf.level = CL, type = binom.conf)
  tab$Sensitivity <- tmp.CI$proportion
  tab$Se.inf.cl <- tmp.CI$lower
  tab$Se.sup.cl <- tmp.CI$upper

  tmp.CI <- binom.CI(tab$TN, ND, conf.level = CL, type = binom.conf)
  tab$Specificity <- tmp.CI$proportion
  tab$Sp.inf.cl <- tmp.CI$lower
  tab$Sp.sup.cl <- tmp.CI$upper

  tmp.CI <- binom.CI(tab$TP, (tab$TP + tab$FP),
                     conf.level = CL, type = binom.conf)
  tab$PPV <- tmp.CI$proportion
  tab$PPV.inf.cl <- tmp.CI$lower
  tab$PPV.sup.cl <- tmp.CI$upper

  tmp.CI <- binom.CI(tab$TN, (tab$TN + tab$FN),
                     conf.level = CL, type = binom.conf)
  tab$NPV <- tmp.CI$proportion
  tab$NPV.inf.cl <- tmp.CI$lower
  tab$NPV.sup.cl <- tmp.CI$upper

  tab$PLR <- tab$Sensitivity / (1 - tab$Specificity)
  tab$PLR.inf.cl <- exp(log(tab$PLR) - (qnorm(1 - ( (1 - CL) / 2),
                                mean = 0, sd = 1)) * sqrt( (1 - tab$Sensitivity)
                                / (D * tab$Specificity) + (tab$Specificity)
                                / (ND * (1 - tab$Specificity))))
  tab$PLR.sup.cl <- exp(log(tab$PLR) + (qnorm(1 - ( (1 - CL) / 2),
                                mean = 0, sd = 1)) * sqrt( (1 - tab$Sensitivity)
                                / (D * tab$Specificity) + (tab$Specificity)
                                / (ND * (1 - tab$Specificity))))

  tab$NLR <- (1 - tab$Sensitivity) / tab$Specificity
  tab$NLR.inf.cl <- exp(log(tab$NLR) - (qnorm(1 - ( (1 - CL) / 2), mean = 0,
                    sd = 1)) * sqrt(tab$Sensitivity / (D * (1 -
                    tab$Sensitivity)) + (1 - tab$Specificity) / (ND *
                    (tab$Specificity))))
  tab$NLR.sup.cl <- exp(log(tab$NLR) + (qnorm(1 - ( (1 - CL) / 2), mean = 0,
                    sd = 1)) * sqrt(tab$Sensitivity / (D * (1 -
                    tab$Sensitivity)) + (1 - tab$Specificity) / (ND *
                    (tab$Specificity))))

  if (reverse) {
    tab$test.values <- abs(tab$test.values)
    tab <- tab[order(tab$test.values), ]
  }

  output <- list(table = tab, sample.size = sample.size,
                 sample.prevalence = sample.prevalence)
  invisible(output)
}
