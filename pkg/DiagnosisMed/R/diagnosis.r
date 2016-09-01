#' Diagnostic test accuracy evaluation
#'
#' @name diagnosis
#'
#' @description \code{diagnosis} estimate sensitivity, specificity, predictive values, likelihood ratios, area under ROC curve and other validity measures for binary diagnostic test evaluation. It accepts as input either columns from a dataset or vectors, a 2 x 2 table or numbers representing true positives, false negatives, false positives and true negatives. \code{plot} for \code{diagnosis} draw a simple nomogram and \code{summary} for \code{diagnosis} creates an output in a table format that allows the output to be easily exported to a spreadsheet.
#'
#' @param a A number representing True Positives from a 2x2 table. Also, TP could either a 2x2 table (see below) or a column in a dataset representing the reference standard.
#'
#' @param b A number representing False Negatives from a 2x2 table. If TP is a column in a dataset FN should also be a columns in a dataset, however representing the index test.
#'
#' @param c A number representing False Positives from a 2x2 table.
#'
#' @param d A number representing True Negatives from a 2x2 table.
#'
#' @param CL Confidence limits for confidence intervals. Must be a numeric value between 0 and 1. Default is 0.95.
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
#'
#' @export
diagnosis <- function(a,b=NULL,c=NULL,d=NULL,CL=0.95,print=TRUE,plot=FALSE){
  #require(epitools)
  if(is.numeric(a)){
       if(all(length(a)==1 & length(b)==1 & length(c)==1 & length(d)==1 & !is.matrix(a))){
         reference.name <- 'Not informed'
         index.name <- 'Not informed'
         tab<-as.table(cbind(rbind(d,c),rbind(b,a)))
         dimnames(tab)<-list(index.test=c("negative","positive"),reference.standard=c("negative","positive"))
         TN<-d
         FN<-b
         FP<-c
         TP<-a
       }
       if(all(length(a) > 1 & length(b) > 1 & is.null(c) & is.null(d) & !is.matrix(a))){
          if(any(is.na(a),is.na(b))){stop('There are NAs either in index test or reference standard. Consider removing or inputing!')}
          if(nlevels(as.factor(a))!=2 | nlevels(as.factor(b))!=2){
             stop('It seems there are more levels then 0 and 1.')
          }
          if(!all(levels(as.factor(a))==c(0,1) & levels(as.factor(b))==c(0,1))){
             stop('Either the index test or the reference test is not correctly coded. 0 and 1 were expected!')
          }
          else{reference.name <- deparse(substitute(a))
               index.name <- deparse(substitute(b))
               tab<-table(b,a,dnn=c(deparse(substitute(b)),deparse(substitute(a))))
               TN<-tab[1,1]
               FN<-tab[1,2]
               FP<-tab[2,1]
               TP<-tab[2,2]
          }
       }
       if(any(is.table(a) | is.matrix(a))){
          if(!all(dim(a)==c(2,2))){
             stop('It seems the inputed table is not 2x2. Check your table output.')
          }
          else{tab<-a
              if(is.null(dimnames(tab))){
                 reference.name <- 'Not informed'
                 index.name <- 'Not informed'
                 dimnames(tab) <- list(index.test = c('negative','positive'),reference.standard = c('negative','positive'))
              }
              else{
                 reference.name <- names(dimnames(tab)[2])
                 index.name <- names(dimnames(tab)[1])
              }
              TN<-tab[1,1]
              FN<-tab[1,2]
              FP<-tab[2,1]
              TP<-tab[2,2]
          }
       }
  }
  if(all(any(is.factor(a) | is.character(a)) & any(is.factor(b) | is.character(b)) & !is.matrix(a))){
     if(any(is.na(a) | is.na(b))){
        stop('There seem to be NAs either in the reference standard or index test. Consider removing or inputing!')
     }
     if(nlevels(as.factor(a))!=2 | nlevels(as.factor(b))!=2){
        stop('It seems there are more levels then negative/absence and positive/presence.')
     }
     if(!all(levels(as.factor(a))==c("negative","positive") & levels(as.factor(b))==c("negative","positive")) &
        !all(levels(as.factor(a))==c("absence","presence") & levels(as.factor(b))==c("absence","presence"))){
          stop('It seems categories are not correctly coded in either the reference or index test.')
     }
     else{reference.name <- deparse(substitute(a))
          index.name <- deparse(substitute(b))
          tab<-table(b,a,dnn=c(deparse(substitute(b)),deparse(substitute(a))))
          TN<-tab[1,1]
          FN<-tab[1,2]
          FP<-tab[2,1]
          TP<-tab[2,2]
     }
  }

  tabmarg<-addmargins(tab)
  Conf.limit<-CL
  # sample size
  n<-sum(tab)
  # prevalence
  p<-(TP+FN)/n
  # sensitivity and confidence limits
  Se<-TP/(TP+FN)
  Se.cl<-as.numeric(binom.CI(TP, TP+FN, conf.level = CL)[4:5])
  # especificity and confidence limits
  Sp<-TN/(FP+TN)
  Sp.cl<-as.numeric(binom.CI(TN, FP+TN, conf.level = CL)[4:5])
  # positive and negative likelyhood ratios and confidence limits
  PLR<-Se/(1-Sp)
  # LR confidence limists inspired in epi.tests{epiR}
  PLR.inf.cl<-exp(log(PLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-Se)/(
    (TP+FN)*Sp)+(Sp)/((FP+TN)*(1-Sp))))
  PLR.sup.cl<-exp(log(PLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-Se)/(
    (TP+FN)*Sp)+(Sp)/((FP+TN)*(1-Sp))))
  NLR<-(1-Se)/Sp
  NLR.inf.cl<-exp(log(NLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((Se)/((TP+
    FN)*(1-Se))+(1-Sp)/((FP+TN)*(Sp))))
  NLR.sup.cl<-exp(log(NLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((Se)/((TP+
    FN)*(1-Se))+(1-Sp)/((FP+TN)*(Sp))))
  #accuracy and confidence limits
  accu<-(TP+TN)/n
  accu.cl<-as.numeric(binom.CI(TP+TN, n, conf.level = CL)[4:5])
  # positive and negative predictive values and confidence limits
  PPV<-TP/(TP+FP)
  PPV.cl<-as.numeric(binom.CI(TP, TP+FP, conf.level = CL)[4:5])
  NPV<-TN/(TN+FN)
  NPV.cl<-as.numeric(binom.CI(TN, TN+FN, conf.level = CL)[4:5])
  # diagnostic odds ratio and confidence limits
  OR<-fisher.test(tab,conf.level=CL)
  DOR<-unname(OR$estimate)
  #DOR<-(TP*TN)/(FP*FN)
  DOR.inf.cl<-OR$conf.int[1]
  DOR.sup.cl<-OR$conf.int[2]
  # error rate and error trade
  #ER<-((FN/(FN+TN))*p)+(((FP/(FP+TP))*(TN+FP))
  ER<-(FN+FP)/n
  ER.cl<-as.numeric(binom.CI(FN+FP, n, conf.level = CL)[4:5])
  ET<-(FN/FP)
  # pre-test and pos-test odds (to do)
  # area under ROC curve
  AUC<-(Se+Sp)/2
  if(plot==TRUE){
    plot(1-Sp,Se,xlim=c(0,1),ylim=c(0,1))
    segments(0,0,1-Sp,Se,col="red")
    segments(1-Sp,Se,1,1,col="red")
    grid()
  }
  #if(plot==FALSE)
  #  {ROC<-roc.from.table(tab, graph = FALSE)}
  # gives same results as AUC<-(Se+Sp)/2
  Youden<-Se+Sp-1
  Youden.inf.cl<-Youden-qnorm((1-CL)/2)*sqrt(((Se * (1 - Se))/(TP+FN) + ((Sp * (1 - Sp))/(TN+FP))))
  Youden.sup.cl<-Youden+qnorm((1-CL)/2)*sqrt(((Se * (1 - Se))/(TP+FN) + ((Sp * (1 - Sp))/(FP+TN))))
#  rm(ROC)
  rm(tab)
  # results evaluations
  output <- list(tabmarg=tabmarg,n=n,p=p,Se=Se,Se.cl=Se.cl,Sp=Sp,Sp.cl=Sp.cl,PLR=PLR,
    PLR.inf.cl=PLR.inf.cl,PLR.sup.cl=PLR.sup.cl,NLR=NLR,NLR.inf.cl=NLR.inf.cl,
    NLR.sup.cl=NLR.sup.cl,accu=accu,accu.cl=accu.cl,PPV=PPV,PPV.cl=PPV.cl,NPV=NPV,NPV.cl=NPV.cl,
    DOR=DOR,DOR.inf.cl=DOR.inf.cl,DOR.sup.cl=DOR.sup.cl,ET=ET,ER=ER,ER.cl=ER.cl,
    Youden=Youden,Youden.inf.cl=Youden.inf.cl,Youden.sup.cl=Youden.sup.cl,AUC=AUC,
    Conf.limit=Conf.limit,reference.name=reference.name,index.name=index.name)
  class(output) <- "diag"
  if(print==TRUE)  {print(output)}
  output
}
