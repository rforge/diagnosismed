
SS <- function(ref, test, reverse = "auto", CL = 0.95, binom.conf = "wilson"){
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
  if (reverse == "auto") {
    if (mean(test[which(ref == 0)]) > mean(test[which(ref == 1)]) || median(test[which(ref == 0)]) > median(test[which(ref == 1)])) {
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
  sample.prevalence <- D / sample.size
  
  
  # Taking the rownames of the test.table to be results first column
  test.values <- as.numeric(rownames(test.table))
  test.diag.table <- as.data.frame(test.values)
  test.diag.table$D <- test.table[,2] 
  test.diag.table$ND <- test.table[,1]  
  test.diag.table$TP <- sapply(1:nrow(test.table), function(i)sum(test.table[i:nrow(test.table),2]))
  if (test.table[1, 2] == 0) {fFN <- 0} else {fFN <- 1}
  test.diag.table$FN <- c(as.integer(fFN), sapply(2:nrow(test.table), function(i)sum(test.table[1:(i - 1), 2])))
  test.diag.table$FP <- sapply(1:nrow(test.table), function(i)sum(test.table[i:nrow(test.table),1]))
  if (test.table[1, 1] == 1) {fTN <- 0} else {fTN <- 1}
  test.diag.table$TN <- c(as.integer(fTN), sapply(2:nrow(test.table), function(i)sum(test.table[1:(i - 1), 1])))
  
  tmp.CI <- binom.CI(test.diag.table$TP, D, conf.level = CL, type = binom.conf)
  test.diag.table$Sensitivity <- tmp.CI$proportion
  test.diag.table$Se.inf.cl <- tmp.CI$lower
  test.diag.table$Se.sup.cl <- tmp.CI$upper
  
  tmp.CI <- binom.CI(test.diag.table$TN, ND, conf.level = CL, type = binom.conf)
  test.diag.table$Specificity <- tmp.CI$proportion
  test.diag.table$Sp.inf.cl <- tmp.CI$lower
  test.diag.table$Sp.sup.cl <- tmp.CI$upper
  
  tmp.CI <- binom.CI(test.diag.table$TP, (test.diag.table$TP + test.diag.table$FP), conf.level = CL, type = binom.conf)
  test.diag.table$PPV <- tmp.CI$proportion
  test.diag.table$PPV.inf.cl <- tmp.CI$lower
  test.diag.table$PPV.sup.cl <- tmp.CI$upper
  
  tmp.CI <- binom.CI(test.diag.table$TN, (test.diag.table$TN + test.diag.table$FN), conf.level = CL, type = binom.conf)
  test.diag.table$NPV <- tmp.CI$proportion
  test.diag.table$NPV.inf.cl <- tmp.CI$lower
  test.diag.table$NPV.sup.cl <- tmp.CI$upper
  
  test.diag.table$PLR <- test.diag.table$Sensitivity / (1 - test.diag.table$Specificity)
  test.diag.table$PLR.inf.cl <- exp(log(test.diag.table$PLR) - (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - test.diag.table$Sensitivity) / (D * test.diag.table$Specificity) + (test.diag.table$Specificity) / (ND * (1 - test.diag.table$Specificity))))
  test.diag.table$PLR.sup.cl <- exp(log(test.diag.table$PLR) + (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - test.diag.table$Sensitivity) / (D * test.diag.table$Specificity) + (test.diag.table$Specificity) / (ND * (1 - test.diag.table$Specificity))))
  
  test.diag.table$NLR <- (1 - test.diag.table$Sensitivity) / test.diag.table$Specificity
  test.diag.table$NLR.inf.cl <- exp(log(test.diag.table$NLR) - (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((test.diag.table$Sensitivity) / (D * (1 - test.diag.table$Sensitivity)) + (1 - test.diag.table$Specificity) / (ND * (test.diag.table$Specificity))))
  test.diag.table$NLR.sup.cl <- exp(log(test.diag.table$NLR) + (qnorm(1 - ((1 - CL) / 2), mean = 0,sd = 1)) * sqrt((test.diag.table$Sensitivity) / (D * (1 - test.diag.table$Sensitivity)) + (1 - test.diag.table$Specificity) / (ND * (test.diag.table$Specificity))))
  
  if (reverse) {
    test.diag.table$test.values <- abs(test.diag.table$test.values)
    test.diag.table <- test.diag.table[order(test.diag.table$test.values),]
  }

  output <- list(table = test.diag.table, sample.size = sample.size, sample.prevalence = sample.prevalence)
  invisible(output)
}

# x <- SS(pop_data$reference,pop_data$test) #Ok
# test <- c(rnorm(200,95,15),rnorm(200,75,20))
# ref <- c(rep(1,200),rep(0,200))
# cbind(test,ref)
# x <- SS(ref,test) #Ok
# plot(1-x$Specificity,x$Sensitivity)

# Making a smoothed SS object with a Neural Network modeling--------------------

NN.SS <- function(x,
                  t.max = NULL,
                  t.min = NULL,
                  precision = .01,  
                  n.neurons = c(1,5,1),
                  learning.rate.global = 1e-2,
                  momentum.global = 0.3,
                  error.criterium = "LMS",
                  Stao = NA,
                  hidden.layer = "sigmoid",
                  output.layer = "sigmoid",
                  method = "ADAPTgdwm",
                  report = FALSE,
                  show.step = 5000,
                  n.shows = 1,
                  CL = .95,
                  pop.prevalence = NULL){
  # Fitting the Neural Network
  net <- newff(n.neurons = n.neurons, learning.rate.global = learning.rate.global, momentum.global = momentum.global, error.criterium = error.criterium, Stao = Stao, hidden.layer = hidden.layer, output.layer = output.layer, method = method)
  
  # Training the Sensitivity
  net.Se <- train(net, P = x$table$test.values, T = x$table$Sensitivity, error.criterium = error.criterium, report = report, show.step = show.step, n.shows = n.shows)
  
  # Training the Specificity
  net.Sp <- train(net, P = x$table$test.values, T = x$table$Specificity, error.criterium = error.criterium, report = report, show.step = show.step, n.shows = n.shows)
  
  # Defining the values where the NN will be simulated
  if (!is.null(pop.prevalence)) {
    x$sample.prevalence <- pop.prevalence
  }
  if (is.null(t.max)) {
    t.max <- max(x$table$test.values)
  }
  if (is.null(t.min)) {  
    t.min <- min(x$table$test.values)
  }
  test.values = seq(t.min, t.max, precision)
  if (length(test.values) < 200) {
    warning("The number of tests values to simulate the neural network is too low, below 200. \n Perhaps changing t.max, t.min or precision will solve.")
  }
  # if(length(test.values) > 2000){
  #   warning("The number of tests values to simulate the neural network is unecessarily high, above 2000. \n Perhaps changing t.max, t.min or precision will solve.")
  # }
  NN <- data.frame(test.values = test.values)
  
  # Simulating the NN
  NN$Sensitivity <- as.numeric(sim(net.Se$net, NN$test.values))
  NN$Se.inf.cl <- NN$Sensitivity - qnorm(1 - (1 - CL) / 2) * sqrt(((NN$Sensitivity * (1 - NN$Sensitivity))/(x$sample.size * x$sample.prevalence)))
  NN$Se.sup.cl <- NN$Sensitivity + qnorm(1 - (1 - CL) / 2) * sqrt(((NN$Sensitivity * (1 - NN$Sensitivity))/(x$sample.size * x$sample.prevalence)))
  NN$Specificity <- as.numeric(sim(net.Sp$net, NN$test.values))
  NN$Sp.inf.cl <- NN$Specificity - qnorm(1 - (1 - CL) / 2) * sqrt(((NN$Specificity * (1 - NN$Specificity))/(x$sample.size * (1 - x$sample.prevalence))))
  NN$Sp.sup.cl <- NN$Specificity + qnorm(1 - (1 - CL) / 2) * sqrt(((NN$Specificity * (1 - NN$Specificity))/(x$sample.size * (1 - x$sample.prevalence))))
  NN$PLR <- NN$Sensitivity/(1 - NN$Specificity) 
  NN$PLR.inf.cl <- exp(log(NN$PLR) - (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - NN$Sensitivity) / ((x$sample.size * x$sample.prevalence) * NN$Specificity) + (NN$Specificity) / ((x$sample.size * (1 - x$sample.prevalence)) * (1 - NN$Specificity))))
  NN$PLR.sup.cl <- exp(log(NN$PLR) + (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - NN$Sensitivity) / ((x$sample.size * x$sample.prevalence) * NN$Specificity) + (NN$Specificity) / ((x$sample.size * (1 - x$sample.prevalence)) * (1 - NN$Specificity))))
  #NN$NLR <- (1-NN$Specificity)/NN$Sensitivity
  #NN$NLR.inf.cl <- exp(log(NN$NLR)-(qnorm(1-((1-(1-CL))/2),mean=0,sd=1))*
  #     sqrt((NN$Sensitivity)/((x$sample.size * x$sample.prevalence)*(1-NN$Sensitivity))+(1-NN$Specificity)/
  #     ((x$sample.size * (1-x$sample.prevalence))*(NN$Specificity))))
  #NN$NLR.sup.cl <- exp(log(NN$NLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((NN$Sensitivity)/
  #     ((x$sample.size * x$sample.prevalence)*(1-NN$Sensitivity))+(1-NN$Specificity)/((x$sample.size * 
  #     (1-x$sample.prevalence))*(NN$Specificity))))
  NN$TP <- NN$Sensitivity * x$sample.prevalence
  NN$FN <- (1 - NN$Sensitivity) * x$sample.prevalence
  NN$FP <- (1 - NN$Specificity) * (1 - x$sample.prevalence)
  NN$TN <- NN$Specificity * (1 - x$sample.prevalence)
  output <- list(table = NN, sample.size = x$sample.size, sample.prevalence = x$sample.prevalence)
  output
}
# NN.SeSp <- NN.SS(mytest, precision = .009)
# plot(NN.SeSp$table$test.values,NN.SeSp$table$Sensitivity,type = "l") ; lines(NN.SeSp$table$test.values,NN.SeSp$table$Specificity)
# lines(mytest$table$test.values,mytest$table$Sensitivity) ; lines(mytest$table$test.values,mytest$table$Specificity)


# Binormal SS object -------------------------------------------------------------------
BN.SS <- function(ref, test, CL = 0.95, t.max = NULL, t.min = NULL, precision = 0.01){
  # Warning section ...
  if(any(is.na(test) | is.na(ref))){
    stop('It seems there are NAs either in the index test or in the reference test. Consider imputing or removing NAs!')
  }
  if(any(levels(as.factor(ref)) != c(0,1))){
    stop("Your reference standard must be coded as 0 (absence) and 1 (presence). Check reference categories!")
  }
  if(!is.numeric(CL) || CL < 0 || CL > 1){
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
  sample.prevalence <- D / sample.size

  # Defining the values where the NN will be simulated
  if (is.null(t.max)) {
    t.max <- max(test)
  }
  if(is.null(t.min)){  
    t.min <- min(test)
  }
  test.values = seq(t.min, t.max, precision)
  if(length(test.values) < 200){
    warning("The number of tests values to simulate the binormal ROC analysis is too low, below 200. \n Perhaps changing t.max, t.min or precision will solve.")
  }
  if(length(test.values) > 2000){
    warning("The number of tests values to simulate the binormal ROC analysis is unecessarily high, above 2000. \n Perhaps changing t.max, t.min or precision will solve.")
  }
  # obtain b and a to obtain the scores as bi-normal curves
  b <- sd0 / sd1
  a <- (m1 - m0) / sd1
  # x <- seq(0, 1, length = 500)
  # y <- 1 - pnorm(b * qnorm(1 - x) - a)
  # cp = seq(min(m0 - 3 * sd0, m1 - 3 * sd1), max(m0 + 3 * sd0, m1 + 3 * sd1), length.out = 500) # potential cut-points
  
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
  BN$PLR <- BN$Sensitivity/(1-BN$Specificity) 
  BN$PLR.inf.cl <- exp(log(BN$PLR) - (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - BN$Sensitivity) / ((sample.size * sample.prevalence) * BN$Specificity) + (BN$Specificity) / ((sample.size * (1 - sample.prevalence)) * (1 - BN$Specificity))))
  BN$PLR.sup.cl <- exp(log(BN$PLR) + (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - BN$Sensitivity) / ((sample.size * sample.prevalence) * BN$Specificity) + (BN$Specificity) / ((sample.size * (1 - sample.prevalence)) * (1 - BN$Specificity))))
  BN$TP <- BN$Sensitivity * sample.prevalence
  BN$FN <- (1 - BN$Sensitivity) * sample.prevalence
  BN$FP <- (1 - BN$Specificity) * (1 - sample.prevalence)
  BN$TN <- BN$Specificity * (1 - sample.prevalence)
  
  output <- list(table = BN, sample.size = sample.size, sample.prevalence = sample.prevalence)
  output
}

