#' Function for the description of the qualities of one or two decision thresholds or threshold.
#'
#' @name quality.threshold
#'
#' @description This function can be used for both dichotomization (single threshold or cut-point) methods and for trichotomization (two thresholds or cut-points) methods. In the case of the Uncertain Interval trichotomization method, it provides descriptive statistics for the test scores outside the Uncertain Interval. For the TG-ROC trichotomization method it provides the descriptive statistics for TG-ROC's Valid Range.
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param threshold The decision threshold of a dichotomization method, or the lower decision threshold of a trichotomization method.
#' @param threshold.upper (default = NULL). The upper decision threshold of a trichotomization method. When NULL, the test scores are dichotomized.
#' @param model The model to use. Default = 'kernel'.
#'
#' @return{ A list of}
#' \describe{
#'   \item{table}{The confusion table of {diag x ref}, where diag is the diagnosis based on the test, when applying the threshold(s). The reference standard (ref) has categories 0 and 1, while the diagnosis based on the test scores (diag) has categories 0 and 1 in the case of applying a single threshold (dichotomization), and the categories 0, NA and 1 in the case of trichotomization. In the case of the Uncertain Interval trichotomization method, the row NA shows the count of test scores within the Uncertain Interval. When applying the trichotomization method TG-ROC, the row NA shows the count of the test scores within the Intermediate Range. Table cell {0, 0} shows the True Negatives (TN), cell {0, 1} shows the False Negatives (FN), cell {1, 0} shows the False Positives (FP), and cell {1, 1} shows the True Positives (TP).}
#'   \item{cut}{The values of the threshold(s).}
#'   \item{indices}{A named vector, with the following statistics for the test-scores with diagnosis 0 or 1:
#'     \itemize{
#'       \item{prevalence: }{Diagnosable patients with the targeted condition / Total diagnosable sample = (TP+FN)/(TN+FP+FN+TP)}
#'       \item{correct.classification.rate (or Accuracy): }{(TP+TN)/(TN+FP+FN+TP)}
#'       \item{balance.correct.incorrect : }{(TP+TN)/(FP+FN)}
#'       \item{specificity: }{TN/(TN+FN)}
#'       \item{sensitivity: }{TP/(TP+FN)}
#'       \item{negative.predictive.value: }{TN/(TN+FN)}
#'       \item{positive.predictive.value: }{TP/(TN+FN)}
#'       \item{neg.likelihood.ratio: }{(1-sensitivity)/specificity}
#'       \item{pos.likelihood.ratio: }{sensitivity/(1-specificity)}
#'       \item{concordance: }{The probability that a random chosen patient with the condition is correctly ranked higher than a randomly chosen patient without the condition. Equal to AUC, with for the more certain interval a higher outcome than the overall concordance.}
#'
#'   }
#' }
#' }
#'
#' @examples
#' # A simple test
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1))
#' ua = ui.nonpar(ref, test)
#' quality.threshold(ref, test, threshold=ua[1], threshold.upper=ua[2])

#' @export
quality.threshold <- function(ref, test, threshold, threshold.upper=NULL, model = c('kernel', 'binormal', 'ordinal')){
  # .Deprecated('quality.treshold', msg = 'Deprecated. Replaced by quality.mci')

  model <- match.arg(model)
  df=check.data(ref, test, model=model)
  ref=df$ref
  test=df$test
  ia = !is.null(threshold.upper)
  # raw=T
  threshold=unname(unlist(threshold)) # threshold=ua[1]
  threshold.upper=unname(unlist(threshold.upper)) # threshold.upper=NULL
  certain.sel=rep(TRUE, length(ref))
  y.hat=rep(1, length(ref)) # init

  ND0=0
  ND1=0
  if (ia) {
    if (threshold.upper < threshold) {
      temp=threshold; threshold=threshold.upper; threshold.upper=temp
    }
    certain.sel = (test < threshold) | (test > threshold.upper)
    uncertain.obs = sum(!certain.sel)
    # if (!raw){
    #   uncertainty.all = mean((ref-test)^2)
    #   uncertainty.certain.interval = mean((ref[certain.sel]-test[certain.sel])^2)
    #   uncertainty.uncertain.interval = mean((ref[!certain.sel]-test[!certain.sel])^2)
    # }
    y.hat[!certain.sel]= NA
    ND0 = sum(is.na(y.hat) & ref==0)
    ND1 = sum(is.na(y.hat) & ref==1)
  }

  y.hat[test < threshold]=0  # change <=
  # table(y.hat,ref)

  TN = sum(y.hat==0 & ref==0, na.rm=T)
  FN = sum(y.hat==0 & ref==1, na.rm=T)
  FP = sum(y.hat==1 & ref==0, na.rm=T)
  TP = sum(y.hat==1 & ref==1, na.rm=T)

  prevalence=(TP+FN)/(TP+FP+FN+TN)
  sensitivity = TP/(TP+FN)
  specificity = TN/(FP+TN)
  positive.predictive.value=TP/(TP+FP)
  negative.predictive.value=TN/(FN+TN)
  correct.classification.rate=(TP+TN)/(TP+FP+FN+TN)
  balance.correct.incorrect=(TP+TN)/(FP+FN)
  likelihood.ratio.pos = sensitivity /(1-specificity)
  likelihood.ratio.neg = (1-sensitivity)/specificity

  o = outer(test[certain.sel & ref==1], test[certain.sel & ref==0], "-")
  cstat=mean((o>0) + .5*(o==0))

  if (ia) {
    t = matrix(data=c(TN, FN, ND0, ND1, FP, TP), ncol=2, byrow=T,
               dimnames=list(diag=c('0', 'NA', '1'), ref=c('0', '1')))
    cut=c(threshold.lower=threshold, threshold.upper=threshold.upper)
    # if (!raw) {uncertainty = c(uncertain.obs=uncertain.obs, uncertainty.all=uncertainty.all,
    #                            uncertainty.certain.interval=uncertainty.certain.interval,
    #                            uncertainty.uncertain.interval=uncertainty.uncertain.interval)
    # }else{
    #   uncertainty=c(uncertain.obs=uncertain.obs)
    # }
  }else{
    t = matrix(data=c(TN, FN, FP, TP), ncol=2, byrow=T,
               dimnames=list(y.hat=c('0', '1'), ref=c('0', '1')))
    cut=c(threshold=threshold)
    uncertainty=NA
  }
  return(list(table=addmargins(t), cut=cut,
              # uncertainty = uncertainty,
              indices=c(prevalence=prevalence,
                        correct.classification.rate=correct.classification.rate,
                        balance.correct.incorrect=balance.correct.incorrect,
                        specificity =specificity,
                        sensitivity =sensitivity,
                        negative.predictive.value=negative.predictive.value,
                        positive.predictive.value=positive.predictive.value,
                        neg.likelihood.ratio = likelihood.ratio.neg,
                        pos.likelihood.ratio = likelihood.ratio.pos,
                        concordance=cstat)))
}
