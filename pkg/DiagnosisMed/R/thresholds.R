# Colection of functions to estimate trhesholds

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

# Maximizing the accuracy threshold where x is a SS object----------------------
max.Accuracy <- function(x){
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
# max.Accuracy(x) ; max.Accuracy(NN.SeSp)

# Maximizing Diagnostic Odds Ratio threshold where x is a SS object------------
max.DOR <- function(x){
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
# max.DOR(x)

# Minimizing error rate threshold where x is a SS object------------------------
min.Error <- function(x){
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
# min.Error(x)

# Maximizing the accuracy area threshold where x is a SS object------------------
max.Accuracy.area <- function(x){
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
    warning(paste0("Accuracy area reaches its maximum at the following test values: ",toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max Accuracy area"
  output
}
# max.Accuracy.area(x); max.Accuracy.area(NN.SeSp)

# Maximizing the Youden J index threshold where x is a SS object----------------
max.Youden <- function(x){
  x$table$Youden <- x$table$Sensitivity + x$table$Specificity - 1
  condition <- which(x$table$Youden == max(x$table$Youden))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Youden J index reaches its maximum at the following test values: ",toString(warn),".\n  The one closest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max Youden"
  output
}
# max.Youden(x)

# Minimizing The ROC 0 1 distance threshold where x is a SS object-------------
 min.ROCdist <- function(x){
  x$table$MinRocDist <- (x$table$Specificity - 1) ^ 2 + (1 - x$table$Sensitivity) ^ 2
  condition <- which(x$table$MinRocDist == min(x$table$MinRocDist))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Minimum ROC distance reaches its minimum at the following test values: ",toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Min ROC distance"
  output
}
# min.ROCdist(x) 

# maximizing Efficiency threshold where x is a SS object------------------------
max.Efficiency <- function(x, pop.prevalence = NULL){
  if (is.null(pop.prevalence)) {
    pop.prevalence <- x$sample.prevalence
  }
  x$table$Efficiency <- (x$table$Sensitivity * (pop.prevalence)) + ((1 - (pop.prevalence)) * x$table$Specificity)
  condition <- which(x$table$Efficiency == max(x$table$Efficiency))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Maximum Efficiency reaches its maximum at the following test values: ",toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Max Efficiency"
  output
}
# max.Efficiency(x) 

# Minimizing MissClassificatio Cost Term threshold where x is a SS object-------
min.MCT <- function(x, pop.prevalence = NULL, Cost = 1){
  if (is.null(pop.prevalence)) {
    pop.prevalence <- x$sample.prevalence
  }
  x$table$MCT <- (1 - (pop.prevalence)) * (1 - x$table$Specificity) + (Cost * (pop.prevalence)) * (1 - x$table$Sensitivity)
  condition <- which(x$table$MCT == min(x$table$MCT))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    m <- median(warn)
    pos <- which.min(abs(x$table$test.values - m))
    warning(paste0("Minimum Misclassification cost term reaches its minimum at the following test values: ",toString(warn),".\n  The one nearest to the median was chosen!"))
  } else {
    pos <- condition
  }
  output <- x$table[pos, c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
  rownames(output) <- "Min MCT"
  output
}
# min.MCT(x)

# Wraping all threshold functions where x is a SS object------------------------
thresholds <- function(x, pop.prevalence = NULL, Cost = 1){
  output <- rbind(max.Youden(x),
                  max.Accuracy(x),
                  max.Accuracy.area(x),
                  max.DOR(x),
                  min.Error(x),
                  min.ROCdist(x),
                  Se.equals.Sp(x),
                  min.MCT(x, pop.prevalence = pop.prevalence, Cost = Cost),
                  max.Efficiency(x, pop.prevalence = pop.prevalence)
                  )
  output
}
# thresholds(x, Cost = 2) ; thresholds(NN.SeSp)

# Inconclusive thresholds (threechotomization) ---------------------------------
inc.limits <- function(x, Inconclusive = .95){
  # Checking the Se values
  condition <- which.min(abs(Inconclusive - x$table$Sensitivity))
  if(length(condition) > 1){
    warn <- x$table$test.values[condition]
    warning(paste0("Sensitivity matches the minimum required at the following test values: ",toString(warn),". Highest one was chosen!"))
  }
  Se.pos <- condition[length(condition)]

  # Checking the Sp values
  condition <- which.min(abs(Inconclusive - x$table$Specificity))
  if (length(condition) > 1) {
    warn <- x$table$test.values[condition]
    warning(paste0("Specificity matches the minimum required at the following test values: ",toString(warn),". First one was chosen!"))
  }
  Sp.pos <- condition[1]
  # Extracting the test.values and corresponding validity
  output <- rbind(
            x$table[Se.pos,c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")],
            x$table[Sp.pos,c("test.values","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl")]
            )
  rownames(output) <- c("Lower inconclusive","Upper inconclusive")
  output
}
# inc.limits(x)