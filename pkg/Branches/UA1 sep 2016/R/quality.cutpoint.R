#' Function for the description of te qualities of one or two decision thresholds or cutpoint. 
#'
#' @name quality.cutpoint 
#'
#' @description This function can be used for both dichotomization (single cut-point) methods and for trichomization (two cut-points) methods. In the case of the Uncertain Area trichotomization method, it provides descriptive statitistics for the test scores outside the Uncertain Area. For the TG-ROC trichomization method it provides te descriptive statistics for TG-ROC's Valid Range.
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param cutpoint The decision threshold of a dichotomization method, or the lower decision threshold of a trichotomization method.
#' @param cutpoint.higher (default = NULL). The upper decision threshold of a trichotomization method. When NULL, the test cores are dicotomized.
#'
#' @return{ A list of}
#' \describe{
#'   \item{table}{The confusion table of {diag x ref}, where diag is the diagnosis based on the test, when applying the cutpoint(s). The reference standard (ref) has categoires 0 and 1, while the diagnosis based on the test scores (diag) has categories 0 and 1 in the case of applying a single cut-point (dischotomization), and the categories 0, NA and 1 in the case of trichomization. In the case of the Uncertain Area trichomization method, the row NA shows the count of test scores within the Uncertain Area. When applying the trichimization method TG-ROC, the row NA shows the count of the test scores within the Intermediate Range. Table cell {0, 0} shows the True Negatives (TN), cell {0, 1} shows the False Negatives (FN), cell {1, 0} shows the False Positives (FP), and cell {1, 1} shows the True Positives (TP).}
#'   \item{cut}{The values of the cutpoint(s).}
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
#'       
#'   }
#' }
#' }
#' @export
#'
#' @examples
#' # A test with considerable overlap, resulting in a relatively large Uncertain Area
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,2))
#' ua = uncertain.area(ref, test) # with warning message!
#' quality.cutpoint(ref, test, ua[1], ua[2])
quality.cutpoint <- function(ref, test, cutpoint, cutpoint.higher=NULL){
  
  df=check.data(ref, test)
  ref=df$ref
  test=df$test
  ia = !is.null(cutpoint.higher)
  # raw=T
  cutpoint=unname(unlist(cutpoint)) # cutpoint=ua[1]
  cutpoint.higher=unname(unlist(cutpoint.higher)) # cutpoint.higher=NULL
  certain.sel=rep(TRUE, length(ref))
  y.hat=rep(1, length(ref)) # init
  
  ND0=0
  ND1=0
  if (ia) {
    if (cutpoint.higher < cutpoint) {
      temp=cutpoint; cutpoint=cutpoint.higher; cutpoint.higher=temp
    }
    certain.sel = (test < cutpoint) | (test > cutpoint.higher)
    uncertain.obs = sum(!certain.sel)
    # if (!raw){
    #   uncertainty.all = mean((ref-test)^2)
    #   uncertainty.certain.area = mean((ref[certain.sel]-test[certain.sel])^2)
    #   uncertainty.uncertain.area = mean((ref[!certain.sel]-test[!certain.sel])^2)
    # }
    y.hat[!certain.sel]= NA
    ND0 = sum(is.na(y.hat) & ref==0)
    ND1 = sum(is.na(y.hat) & ref==1)
  }
  
  y.hat[test < cutpoint]=0  # change <=
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
  if (ia) {
    t = matrix(data=c(TN, FN, ND0, ND1, FP, TP), ncol=2, byrow=T,
               dimnames=list(diag=c('0', 'NA', '1'), ref=c('0', '1')))
    cut=c(cutpoint.lower=cutpoint, cutpoint.higher=cutpoint.higher)
    # if (!raw) {uncertainty = c(uncertain.obs=uncertain.obs, uncertainty.all=uncertainty.all,
    #                            uncertainty.certain.area=uncertainty.certain.area,
    #                            uncertainty.uncertain.area=uncertainty.uncertain.area)
    # }else{
    #   uncertainty=c(uncertain.obs=uncertain.obs)
    # }
  }else{
    t = matrix(data=c(TN, FN, FP, TP), ncol=2, byrow=T,
               dimnames=list(y.hat=c('0', '1'), ref=c('0', '1')))
    cut=c(cutpoint=cutpoint)
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
                        pos.likelihood.ratio = likelihood.ratio.pos)))
}