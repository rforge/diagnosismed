#' @rdname SS
#' @export
# Binormal SS object ----------------------------------------------------------
BN.SS <- function(ref, test, CL = 0.95, t.max = NULL, t.min = NULL,
                  precision = 0.01, pop.prevalence = NULL){
  # Warning section ...
  if (any(is.na(test) | is.na(ref))) {
    stop('It seems there are NAs either in the index test or in the reference
         test. Consider imputing or removing NAs!')
  }
  if (any(levels(as.factor(ref)) != c(0,1)))  {
    stop("Your reference standard must be coded as 0 (absence) and 1 (presence).
         Check reference categories!")
  }
  if (!is.numeric(CL) || CL < 0 || CL > 1) {
    stop("Confidence limit (CL) must be nemeric between 0 and 1.")
  }

  # first get an estimate of the normal curves
  m0  <- mean(test[which(ref == 0)])
  sd0 <- sd(test[which(ref == 0)])
  m1 <- mean(test[which(ref == 1)])
  sd1 <- sd(test[which(ref == 1)])

  if (m0 > m1 || median(test[which(ref == 0)]) > median(test[which(ref == 1)])) {
    warning("The mean test values of subjects with the condition is lower than from those without the condition.")
  }

  D <- length(test[which(ref == 1)])
  ND <- length(test[which(ref == 0)])
  sample.size <- D + ND
  if (!is.null(pop.prevalence)) {
    sample.prevalence <- pop.prevalence
  } else {
    sample.prevalence <- D / sample.size
  }

  # Defining the values where the NN will be simulated
  if (is.null(t.max)) {
    t.max <- max(test)
  }
  if (is.null(t.min)) {
    t.min <- min(test)
  }
  test.values = seq(t.min, t.max, precision)
  if (length(test.values) < 200) {
    warning("The number of tests values to simulate the binormal ROC analysis is too low, below 200. \n Perhaps changing t.max, t.min or precision will solve.")
  }
  # if(length(test.values) > 2000){
  #   warning("The number of tests values to simulate the binormal ROC analysis is unecessarily high, above 2000. \n Perhaps changing t.max, t.min or precision will solve.")
  # }
  # obtain b and a to obtain the scores as bi-normal curves
  b <- sd0 / sd1
  a <- (m1 - m0) / sd1

  # Obtaining the Se and Sp
  Sp <- 1 - pnorm((m0 - test.values) / sd0)
  Se <- pnorm((m1 - test.values) / sd1)

  # Making a diagnostic table Make homogeneuous SS output
  BN <- data.frame(test.values = test.values, Sensitivity = Se)
  BN$Se.inf.cl <- BN$Sensitivity - qnorm(1 - (1 - CL) / 2) * sqrt(((BN$Sensitivity * (1 - BN$Sensitivity))/(sample.size * sample.prevalence)))
  BN$Se.sup.cl <- BN$Sensitivity + qnorm(1 - (1 - CL) / 2) * sqrt(((BN$Sensitivity * (1 - BN$Sensitivity))/(sample.size * sample.prevalence)))
  BN$Specificity <- Sp
  BN$Sp.inf.cl <- BN$Specificity - qnorm(1 - (1 - CL) / 2) * sqrt(((BN$Specificity * (1 - BN$Specificity))/(sample.size * (1 - sample.prevalence))))
  BN$Sp.sup.cl <- BN$Specificity + qnorm(1 - (1 - CL) / 2) * sqrt(((BN$Specificity * (1 - BN$Specificity))/(sample.size * (1 - sample.prevalence))))
  BN$PLR <- BN$Sensitivity/(1 - BN$Specificity)
  BN$PLR.inf.cl <- exp(log(BN$PLR) - (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - BN$Sensitivity) / ((sample.size * sample.prevalence) * BN$Specificity) + (BN$Specificity) / ((sample.size * (1 - sample.prevalence)) * (1 - BN$Specificity))))
  BN$PLR.sup.cl <- exp(log(BN$PLR) + (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - BN$Sensitivity) / ((sample.size * sample.prevalence) * BN$Specificity) + (BN$Specificity) / ((sample.size * (1 - sample.prevalence)) * (1 - BN$Specificity))))
  BN$TP <- BN$Sensitivity * sample.prevalence
  BN$FN <- (1 - BN$Sensitivity) * sample.prevalence
  BN$FP <- (1 - BN$Specificity) * (1 - sample.prevalence)
  BN$TN <- BN$Specificity * (1 - sample.prevalence)

  output <- list(table = BN, sample.size = sample.size, sample.prevalence = sample.prevalence)
  output
  }

