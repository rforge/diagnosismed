#' Diagnostic test accuracy evaluation
#'
#' @name diagnosis
#'
#' @description \code{diagnosis} estimate sensitivity, specificity, predictive values, likelihood ratios, area under ROC curve and other validity measures for binary diagnostic test evaluation. It accepts as input either columns from a dataset or vectors, a 2 x 2 table or numbers representing true positives, false negatives, false positives and true negatives. \code{plot} for \code{diagnosis} draw a simple nomogram and \code{summary} for \code{diagnosis} creates an output in a table format that allows the output to be easily exported to a spreadsheet.
#'
#' @param TP,FN,FP,TN A number representing True Positives from a 2x2 table. Also, TP could either a 2x2 table (see below) or a column in a dataset representing the reference standard.
#'
#' @param tab A number representing False Negatives from a 2x2 table. If TP is a column in a dataset FN should also be a columns in a dataset, however representing the index test.
#'
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#'
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#'
#' @param CL Confidence limits for confidence intervals. Must be a numeric value between 0 and 1. Default is 0.95.
#'
#' @param CL.type Type of confidence limit. Accepted values are "wilson", "exact", approximate". See \code{\link{binom.CI}}
#'
#' @param print If TRUE, diagnosis will print in the output window the statistics resulted from the 2x2 table. For plot, this will print in the output window a table with all pre-test and its corresponding post-test probabilities.
#'
#' @param plot There are two options of plot. If plot is called in \code{diagnosis}, then a ROC curve of the test under evaluation plotted. If plot is called from an object storing diagnosis output (see example) than a nomogram is plotted. These plots may later be edited, as any other plot, with title, legends etc. Default is FALSE.
#'
#' @param x For plot, print and summary functions, x is an object assigned with diagnois output.
#'
#' @param ... Other options passed to \code{\link[base]{print}}, \code{\link[graphics]{plot.default}} or \code{\link[base]{summary}}.
#'
#' @details In \command{diagnosis}, the values entered must be eather two variables in a data frame, a 2 x 2 table or numbers corresponding to 2 x 2 table cells. If vectors or columns from a dataset, the first one should be the gold standard and the second should be the index test. These two variables must be coded either as numeric - 0 for negative and 1 for a positive test - or with the words "positive" and "negative" or "presence" and "absence". In a older version, there was diagnosisI function that was replace by diagnosis function. The values of a 2 x 2 table can be inputted as: TP is true positive; TN is true negative; FP is false positive and FN is false negative. Sensitivity, Specificity, Predictive values and Accuracy confidence limits rely on binomial distribution, which does not give result outside [0:1] such as normal distribution or asymptotic theory. DOR, Likelihood  ratios and Youden J index confidence limits rely on normal approximation (Wald method for likelihoods). The AUC (area under the ROC curve) is estimated by trapezoidal method (see below). Also, this functions have a summary function wich creates a matrix as a result (identical to the default print option) wich allows to easily export the results to a spreadsheet or to a odt file (with OdfWeave) in a table format (see example). If the input is a 2 x 2 table it should be formated as:
#'
#' \tabular{cll}{
#' \tab TN \tab FN \cr
#' \tab FP \tab TP \cr
#' }
#'
#' \code{plot.diag} will draw a very simple nomogram as many examples from wikipedia \url{http://en.wikipedia.org/wiki/Nomogram}. This is not a generic nomogram as shown in many evidenced based medicine texts, because this one shows only pre-test and post-test variations with a fixed positive likelihood ratio. This likelihood is a statistic from an object created by \code{diagnosis} function. Its usage is the same as applying the Bayes theorem where the pre-test odds times positive likelihood ratio equals the pos-test odd (transforming the odds to probabilities). To use it, draw, with a rule, a vertical line from a desired pre-test  probability, and to find the corresponding post-test probability, draw a horizontal line from the intersection of the curve and the vertical line toward the vertical axis.
#'
#' @return A 2x2 table from which the validity measures are calculated.
#' \itemize{
#' \item Sample size. The number of subjects analyzed.
#' \item Prevalence. The proportion classified as with the target condition by the reference standard.
#' \item Sensitivity. The probability of the test to correctly classify subjects with the target condition (TP/(TP+FN)).
#' \item Specificity. The probability of the test to correctly classify subjects without the target condition (TN/(TN+FP)).
#' \item Predictive values. The probabilities of being with (positive predictive value) (TP/(TP+FP)) or without (negative predictive value) the target condition given a test result (TN/(TN+FN)).
#' \item Likelihood ratios. The probability of test a result in people with the target condition, divided by the probability of the same test result in people without the target condition (PLR = Se/(1-Sp); NLR = (1-Sp)/Se).
#' \item Diagnostic odds ratio. Represents the overall discrimination of a dichotomous test, and is equivalent to the ratio of PLR and NLR.
#' \item Error trade off. Is the amount of false positives traded with false negatives for each decision threshold; here expressed as an odd - for binary results there is only one threshold.
#' \item Error rate. Expresses how many errors we make when we diagnose patients with an abnormal test result as diseased, and those with a normal test result as non-diseased ((FP+FN)/sample size).
#' \item Accuracy. Overall measure that express the capacity of the test to correctly classify subjects with and without the target condition ((TP+TN)/(sample size)).
#' \item Area under ROC curve. Overall measure of accuracy - here the method is the trapezoidal. It gives identical results as (Se+SP)/2.
#' }
#'
#' @references
#' Knotterus. The Evidence Based Clinical Diagnosis; BMJBooks, 2002.
#'
#' Xiou-Hua Zhou, Nancy A Obuchowsky, Donna McClish. Statistical Mehods in diagnostic Medicine; Wiley, 2002.
#'
#' Simel D, Samsa G, Matchar D (1991). Likelihood ratios with confidence: Sample size estimation for diagnostic test studies. Journal of Clinical Epidemiology 44: 763 - 770
#'
#' @seealso \code{\link{LRgraph}}, \code{\link{binom.CI}}
#'
#' @examples
#'  # Simulating a dataset
#'  mydata <- as.data.frame(rbind(
#'      cbind(rep(c("positive"),18),rep(c("negative"),18)),
#'      cbind(rep(c("positive"),72),rep(c("positive"),72)),
#'      cbind(rep(c("negative"),25),rep(c("positive"),25)),
#'      cbind(rep(c("negative"),149),rep(c("negative"),149))
#'      ))
#'  colnames(mydata) <- c('culture','serology')
#'
#'  # A little description of the data set to check if it is ok!
#'  str(mydata)
#'  # Attaching the data set and checking
#'  attach(mydata)
#'
#'  # Running the diagnosis analysis
#'  diagnosis(culture,serology)
#'  detach(mydata)
#'
#'  #Simulating a table
#'  mytable <- matrix(c(149,18,25,72), nrow = 2, ncol=2, byrow=TRUE,
#'                      dimnames = list(enzyme=c('absent','present'),
#'                      citology=c('absent','present')))
#'
#'  # Running analysis from a 2 x 2 table
#'  diagnosis(mytable)
#'
#'  #Inserting values as isolated numbers
#'  diagnosis(72,18,25,149)
#'
#'  #---------------------------------
#'  # Export results to a spreadsheet:
#'  #---------------------------------
#'
#'  # Assing diagnosis to an object
#'  mytest <- diagnosis(364,22,17,211,print=FALSE)
#'
#'  # Assign the summary to an object
#'  mt.sum <- summary(mytest)
#'
#'  # Export to a spreadsheet using csv format - could also work to text with OdfWeave export.
#'  # write.csv(mt.sum,'MytestResults.csv',quote = F,na='')
#'
#'  # Draw a nomogram from a test
#'  plot(mytest)
#'
#'  rm(mydata,mytable,mytest,mt.sum)
#' @import stats
#' @export
diagnosis <- function(tab = NULL, ref = NULL, test = NULL, TN = NULL, FN = NULL, FP = NULL, TP = NULL, reference.name = NULL, index.name = NULL, CL = 0.95, CL.type = c("wilson", "exact", "approximate"), plot = FALSE){
  # Verificando se todos são NULL
  if (all(is.null(TP), is.null(TN), is.null(FP), is.null(FN), is.null(ref), is.null(test), is.null(tab))) {
    stop("A combination of TP and FP and FN and TN, or ref and test, or a valid tab must be provided.")
  }

  # Verificando se a entrada de TP + TN + FP + FN é valida e fazendo a tabela---
  if (!is.null(TP) && !is.null(TN) && !is.null(FP) && !is.null(FN)){
    if (!is.numeric(c(TP, TN, FN, FP)) || length(FP) != 1 || length(TP) != 1 || length(FN) != 1 || length(TN) != 1) {
      stop("When TP and FP and FN and TN are provided, they all must be numeric, each of length 1.")
    }
    if (is.null(reference.name)) {
      reference.name <- 'Not informed'
    }
    if (is.null(index.name)) {
      index.name <- 'Not informed'
    }
    tab <- as.table(cbind(rbind(TN, FP), rbind(FN, TP)))
    dimnames(tab) <- list(index.test = c("Negative", "Positive"), reference.standard = c("Negative", "Positive"))
    # TN <- d
    # FN <- b
    # FP <- c
    # TP <- a
  }

  # Verifcando a entrada de duas variáveis e fazendo a tabela
  if (!is.null(ref) && !is.null(test)) {
    if (any(is.na(ref), is.na(test))) {
      stop('There are NAs either in index test or reference standard. Consider removing or inputing!')
    }
    if (is.factor(ref) || is.factor(test)) {
      if(nlevels(ref) != 2 || nlevels(test) != 2){
        stop('It seems there are more than two levels.')
      }
    }
    if (!is.factor(ref) || !is.factor(test)) {
      if(nlevels(as.factor(ref)) != 2 || nlevels(as.factor(test)) != 2){
        stop('It seems there are more than two levels.')
      }
    }
    if (is.null(reference.name)) {
      reference.name <- deparse(substitute(ref))
    }
    if (is.null(index.name)) {
      index.name <- deparse(substitute(test))
    }
    tab <- table(test, ref, dnn = c(index.name, reference.name))
    TN <- tab[1, 1]
    FN <- tab[1, 2]
    FP <- tab[2, 1]
    TP <- tab[2, 2]
  }

  # Verificando a entrada como tabela
  if (!is.null(tab)) {
    if (any(!is.table(tab) & !is.matrix(tab))) {
      stop("'tab' should be a table or a matrix.")
    }
    if (any(dim(tab) != c(2, 2))) {
      stop("'tab' should be a 2x2 table or a matrix.")
    }
    if (is.null(reference.name)) {
      reference.name <- 'Not informed'
    }
    if (is.null(index.name)) {
      index.name <- 'Not informed'
    }
    if (is.na(names(dimnames(tab)[2])) || is.null(names(dimnames(tab)[2]))) {
      names(dimnames(tab)[2]) <- reference.name
    }
    if (is.na(names(dimnames(tab)[1])) || is.null(names(dimnames(tab)[1]))) {
      names(dimnames(tab)[1]) <- index.name
    }
    TN <- tab[1, 1]
    FN <- tab[1, 2]
    FP <- tab[2, 1]
    TP <- tab[2, 2]
    if (!is.numeric(c(TN, FN, FP, TP))) {
      stop("At least one of the table cells are not numeric.")
    }
  }

  tabmarg <- addmargins(tab)
  Conf.limit <- CL
  # sample size
  n <- sum(tab)
  # prevalence
  p <- (TP + FN) / n
  # sensitivity and confidence limits
  tmp <- binom.CI(TP, TP + FN, conf.level = CL, type = CL.type)
  Se <- tmp$proportion
  Se.inf.cl <- tmp$lower
  Se.sup.cl <- tmp$upper
  # especificity and confidence limits
  tmp <- binom.CI(TN, FP + TN, conf.level = CL, type = CL.type)
  Sp <- tmp$proportion
  Sp.inf.cl <- tmp$lower
  Sp.sup.cl <- tmp$upper
  # positive and negative likelyhood ratios and confidence limits
  PLR <- Se / (1 - Sp)
  # LR confidence limists inspired in epi.tests{epiR}
  PLR.term <- (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((1 - Se) / ((TP + FN) * Sp) + (Sp) / ((FP + TN) * (1 - Sp)))
  PLR.inf.cl <- exp(log(PLR) - PLR.term)
  PLR.sup.cl <- exp(log(PLR) + PLR.term)
  NLR <- (1 - Se) / Sp
  NLR.term <- (qnorm(1 - ((1 - CL) / 2), mean = 0, sd = 1)) * sqrt((Se) / ((TP + FN) * (1 - Se)) + (1 - Sp) / ((FP + TN) * (Sp)))
  NLR.inf.cl <- exp(log(NLR) - NLR.term)
  NLR.sup.cl <- exp(log(NLR) + NLR.term)
  #accuracy and confidence limits
  tmp <- binom.CI(TP + TN, n, conf.level = CL, type = CL.type)
  accu <- tmp$proportion
  accu.inf.cl <- tmp$lower
  accu.sup.cl <- tmp$upper
  # positive and negative predictive values and confidence limits
  tmp <- binom.CI(TP, TP + FP, conf.level = CL, type = CL.type)
  PPV <- tmp$proportion
  PPV.inf.cl <- tmp$lower
  PPV.sup.cl <- tmp$upper
  tmp <- binom.CI(TN, TN + FN, conf.level = CL, type = CL.type)
  NPV <- tmp$proportion
  NPV.inf.cl <- tmp$lower
  NPV.sup.cl <- tmp$upper
  # diagnostic odds ratio and confidence limits
  OR <- fisher.test(tab, conf.level = CL)
  DOR <- unname(OR$estimate)
  #DOR<-(TP*TN)/(FP*FN)
  DOR.inf.cl <- OR$conf.int[1]
  DOR.sup.cl <- OR$conf.int[2]
  # error rate and error trade
  #ER<-((FN/(FN+TN))*p)+(((FP/(FP+TP))*(TN+FP))
  tmp <- binom.CI(FN + FP, n, conf.level = CL, type = CL.type)
  ER <- tmp$proportion
  ER.inf.cl <- tmp$lower
  ER.sup.cl <- tmp$upper
  # ET <- (FN/FP)
  # pre-test and pos-test odds (to do)
  # area under ROC curve
  AUC <- (Se + Sp) / 2
  # Younden J index
  Youden <- Se + Sp - 1
  Youden.term <- qnorm(1 - ((1 - CL) / 2)) * sqrt(((Se * (1 - Se)) / (TP + FN) + ((Sp * (1 - Sp)) / (TN + FP))))
  Youden.inf.cl <-Youden - Youden.term
  Youden.sup.cl <-Youden + Youden.term
  if(plot){
    plot(1 - Sp, Se, xlim = c(0, 1), ylim = c(0, 1))
    segments(0, 0, 1 - Sp, Se, col = "red")
    segments(1 - Sp, Se, 1, 1, col = "red")
    grid()
  }
  results <- data.frame(Estimate = c(n, p, Se, Sp, PPV, NPV, PLR, NLR, DOR, ER, accu, Youden, AUC),
             lower.cl = c(NA, NA, Se.inf.cl, Sp.inf.cl, PPV.inf.cl, NPV.inf.cl, PLR.inf.cl, NLR.inf.cl, DOR.inf.cl, ER.inf.cl, accu.inf.cl, Youden.inf.cl, NA),
             upper.cl = c(NA, NA, Se.sup.cl, Sp.sup.cl, PPV.sup.cl, NPV.sup.cl, PLR.sup.cl, NLR.sup.cl, DOR.sup.cl, ER.sup.cl, accu.sup.cl, Youden.sup.cl, NA))
  rownames(results) <- c('Sample size:','Prevalence:', 'Sensitivity:', 'Specificity:',
                         'Postive predictive value:', 'Negative predictive value:',
                         'Positive likelihood ratio:', 'Negative likelihood ratio:',
                         'Diagnostic Odds Ratio:', 'Error rate:',
                         'Accuracy:', 'Youden J index:', 'Area under ROC curve:')
  # results evaluations
  output <- list(tabmarg = tabmarg, n = n, p = p, Se = Se, Se.inf.cl = Se.inf.cl, Se.sup.cl = Se.sup.cl,
                 Sp = Sp, Sp.inf.cl = Sp.inf.cl, Sp.sup.cl = Sp.sup.cl,
                 PLR = PLR, PLR.inf.cl = PLR.inf.cl, PLR.sup.cl = PLR.sup.cl,
                 NLR = NLR, NLR.inf.cl = NLR.inf.cl, NLR.sup.cl = NLR.sup.cl,
                 accu = accu, accu.inf.cl = accu.inf.cl, accu.sup.cl = accu.sup.cl,
                 PPV = PPV, PPV.inf.cl = PPV.inf.cl, PPV.sup.cl = PPV.sup.cl,
                 NPV = NPV, NPV.inf.cl = NPV.inf.cl, NPV.sup.cl = NPV.sup.cl,
                 DOR = DOR, DOR.inf.cl = DOR.inf.cl, DOR.sup.cl = DOR.sup.cl,
                 ER = ER, ER.inf.cl = ER.inf.cl, ER.sup.cl = ER.sup.cl,
                 Youden = Youden, Youden.inf.cl = Youden.inf.cl, Youden.sup.cl = Youden.sup.cl,
                 AUC = AUC, Conf.limit = Conf.limit, reference.name = reference.name, index.name = index.name, results = results)
  class(output) <- "diag"
  output
}
