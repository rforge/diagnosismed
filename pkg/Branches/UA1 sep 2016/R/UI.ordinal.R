#' Experimental! Function to calculate all possible uncertain intervals of ordinal test results of individuals with (1) and without (0) the targeted condition.
#' @name UI.ordinal
#' @importFrom grDevices rgb
#' @importFrom graphics legend plot
#' @export
#' @description This function calculates the uncertain interval of the test results of the two groups. This function is intended to be used for tests with 20 or less ordened test values
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param adjust The adjust factor that is transferred to the density function. Default value is 2.
#' @param ... Further parameters that are transferred to the density function.
#' @return List of values:
#' \describe{
#'   \item{$Youden}{A vector of statistics concerning the maximized Youden index:}
#'      \itemize{
#'       \item{max.Youden: }{The value of the Maximized Youden Index (= max(tpr - fpr)).}
#'       \item{threshold: }{The threshold associated with the Maximized Youden Index. Test values >= threshold indicate the targeted condition.}
#'       \item{Sp: }{The Specificity of the test when this threshold is applied}
#'       \item{Se: }{The Sensitivity of the test when this threshold is applied}
#'       \item{Acc: }{The Accuracy of the test when this threshold is applied}
#'       }
#'   \item{data.table}{A data.frame with the following columns:}
#'      \itemize{
#'       \item{test: }{The test scores.}
#'       \item{d0: }{The frequencies of the test scores of the norm group.}
#'       \item{d1: }{The frequencies of the test scores of the group with the targeted condition.}
#'       \item{tot: }{The total frequency of each test scores.}
#'      \item{TP: }{The number of True Positives when this test score is used as threshold.}
#'      \item{FP: }{The number of False Positives when this test score is used as threshold.}
#'      \item{fpr: }{The false positive rate when this test score is used as threshold.}
#'      \item{tpr: }{The true positive rate when this test score is used as threshold.}
#'      \item{Y: }{The Youden Index (= tpr - fpr) when this test score is used as threshold.}
#'      }
#'   \item{intersection}{The (rounded) intersection for the distributions of the two groups. Most often these distributions have no true point of intersection and the rounded intersection is an approximation. Often, this equals the Maximized Youden threshold (see Schisterman 2005). }
#'   \item{Uncertain.Interval}{Data frame with the statistics of all possible bounds of the uncertain interval. The columns are the following: }
#'     \itemize{
#'       \item{lowerbound: }{Lower bound of the possible uncertain interval.}
#'       \item{upperbound: }{Upper bound of the possible uncertain interval.}
#'       \item{UI.Sp: }{Specificity of the test scores between and including the lower and upper boundary. Closer to .5 is 'better', that is, more uncertain.}
#'       \item{UI.Se: }{Sensitivity of the test scores between and including the lower and upper boundary. Closer to .5 is 'better', that is, more uncertain.}
#'       \item{UI.ratio: }{The ratio between the number of patients in the uncertain area with and without the condition. Closer to one is better; (.45/.55) < UI.ratio < (.55/.45) or 0.82 < UI.ratio < 1.22 as a rule of fist.}
#'       \item{UI.n: }{Number of patients with test scores between and including the lower and upper boundary.}
#'       \item{MCI.Sp: }{Specificity of the more certain interval, i.e., the test scores lower than the lower boundary and higher than the upper boundary.}
#'       \item{MCI.Se: }{Sensitivity of the test scores lower than the lower boundary and higher than the upper boundary.}
#'       \item{MCI.Acc: }{Accuracy of the test scores lower than the lower boundary and higher than the upper boundary.}
#'       \item{MCI.n: }{Number of patients with test scores lower than the lower boundary and higher than the upper boundary.}
#'       \item{Loss: }{Loss of the test scores lower than the lower boundary and higher than the upper boundary. The total loss is the sum of the loss of the three areas:
#'       lower MCI: the summed number of false positives for each test score, divided by the number of persons that have received that test score.
#'       upper MCI: the summed number of false negatives, divided by the number of persons that have received that test score.
#'       uncertain interval: the sum of the absolute differences in the number of people out of the norm group d0 and the number of persons in the group with the targeted condition (d1) per test score, divided by the total number of persons.}
#'       }
#'   \item{candidates: }{Candidates with a loss lower than the Youden loss which might be considered for the Uncertain Interval }
#'   }
#'   @references{ Schisterman, E. F., Perkins, N. J., Liu, A., & Bondell, H. (2005). Optimal cut-point and its corresponding Youden Index to discriminate individuals using pooled blood samples. Epidemiology, 73-81.
#'
#'   Landsheer, J. A. (2016). Interval of Uncertainty: An alternative approach for the determination of decision thresholds, with an illustrative application for the prediction of prostate cancer. PLOS One.

#'   }
#' @details Warning: This function has not been sufficiently tested or studied.
#'
#' Due to the limited possibilities of short scales, it is more difficult to determine a suitable uncertain interval. For any threshold determination, one needs a representative sample of sufficient size (200 or larger). If there are no test scores below the intersection in the candidate uncertain area, Sp of the Uncertain Interval (UI.Sp) is not available, while UI.Se equals 1.
#' The essential question is always whether the patients with the test scores inside the uncertain interval can be sufficiently distinguished.
#' Discussion of the first example (please run the code first): Visual inspection of the Mixed Barplot shows that distinguishing patients with and without the targeted condition is almost impossible for test scores 2, 3 and 4.
#' Sensitivity and Specificity of the uncertain interval should be not too far from .5. In the first example, the first interval (3:3) has no lower scores than the intersection (3), and therefore Ui.Sp is not available and UI.Se = 1. The UI.ratio indicates that the number of patients with and without the condition are equal in this interval. For these 110 patients, a diagnosis of uncertainty is probably the best choice.
#' The second interval (3:4) has an UI.Sp of .22, which is a large deviation from .5. In this slightly larger interval, the patients with a test score of 3 have a slightly larger probability to belong to the group without the condition. UI.Se is .8. UI.ratio is close to 1, which makes it a feasible candidate.
#' The third interval (2:4) has an UI.Sp of .35 and an UI.Se of .70 and an UI.ratio still close to one.
#' The other intervals show either Se or Sp that deviate strongly from .5, which makes them unsuitable choices.
#' Probably the easiest way to determine the uncertain interval is the interval with minimum loss. This is interval (2:4).
#' The Loss formula is (created by https://www.codecogs.com/latex/eqneditor.php):
#' \deqn{L =\frac{ \left (\sum_{i=l}^{u} \left |d0_{i}-d1_{i}  \right | + \sum_{i=u+1}^{h} d1_{i}+ \sum_{i=1}^{l-1}d0_{i}\right )}{N}}{ L = 1/N * (sum(abs(d0[u:l] - d1[u:l])) + sum(d1[1:(l-1)]) + sum(d0[(u+1):h]))}
#' where \emph{d0} represents the test scores of the norm group, \emph{d1} represents the test scores of the targeted patient group, \emph{l} is the lower limit of the uncertain interval, \emph{u} the upper limit, the first test score is enumerated 1 and the last test score is enumerated \emph{h}. \emph{N} is the total number of all persons with test scores.
#'  \itemize{
#'  \item{\eqn{\sum_{i=l}^{u} \left |d0_{i}-d1_{i}  \right |}{sum(abs(d0[u:l] - d1[u:l])}}{ is the loss in the uncertain interval, that is, the total deviation from equality.}
#'  \item{\eqn{\sum_{i=u+1}^{h} d1_{i}}{sum(d1[1:(l-1)])}}{ is the loss in the lower More Certain Interval, that is, the total of False Negatives, the number of patients with the targeted condition with a test score lower than \emph{l}, and}
#'  \item{\eqn{\sum_{i=u+1}^{h} d0_{i}}{sum(d0[(u+1):h])}}{ is the loss in the lower More Certain Interval, that is, the total of False Positives, the number of patients without the targeted condition with a test score higher than \emph{u}.}
#'
#' The Loss is higher when the deviation from equality is higher in the uncertain area, higher when the number of False Negatives is higher, and higher when the number of False Positives is higher. The loss of a single threshold method equals 1 - its Accuracy.
#' In this example, the minimum Loss is found with interval (2:4). As this agrees with acceptable values for UI.Se, UI.Sp and UI.ratio, this seems the most suitable choice. In this case, the number of patients with test scores 2 to 4 are almost as likely to come from either population.
#' The remaining cases outside the uncertain interval (2:4) show high Accuracy, Specificity and Sensitivity.
#'
#' }
#' @seealso{  \code{\link{plotMD}} for plotting the mixed densities of the test values.
#'             \code{link{stats::density}} for the parameters of the density function. }
#' @examples
#' # A short test with 5 ordinal values
#' test0     = rep(1:5, times=c(165,14,16,55, 10)) # test results norm group
#' test1     = rep(1:5, times=c( 15,11,13,55,164)) # test results of patients
#' ref = c(rep(0, length(test0)), rep(1, length(test1)))
#' test = c(test0, test1)
#' table(ref, test)
#' plotMD(ref, test, breaks=c(min(test)-.5, seq(min(test), max(test), by=1)+.5),
#'   adjust=2) # visual inspection
#' UI.ordinal(ref, test) # determination of the uncertain interval

UI.ordinal <- function(ref, test, adjust=2, ...){
  # pr= sum(ref==1)/length(ref) # prevalence

  tab=as.matrix(table( test, ref))
  # print(addmargins(tab))
  totpos=sum(tab[,2])                          # The total number of positives (one number)
  totneg=sum(tab[,1])                          # The total number of negatives (one number)

  d=data.frame(test=unique(sort(test)), d0=tab[,'0'], d1=tab[,'1'], row.names=1:nrow(tab))
  d$tot=rowSums(tab)                             # Number of patients w/ each test result
  d$TP=unname(rev(cumsum(rev(tab[,2]))))     # Number of true positives, start with maximum
  d$FP=unname(rev(cumsum(rev(tab[,1]))))     # Number of false positives, start with maximum
  # d$TN=unname(cumsum(tab[,1]))
  d$fpr=d$FP/totneg  # 1-d$FP/totneg         # 1 − specificity (false positives)
  d$tpr=d$TP/totpos          # Sensitivity (fraction true positives)
  d$Y = d$tpr-d$fpr # Youden index
  Yt= which.max(d$Y) # Youden threshold as row number
  Acc.Y = (sum(d$d0[1:(Yt-1)])+sum(d$d1[(Yt):nrow(d)]))/ (totpos+totneg)
#  Loss.Y = sum(d$d1[1:(Yt-1)])/totpos + sum(d$d0[(Yt):nrow(d)])/totneg #/ (totpos+totneg)
  # Acc.Y + Loss.Y # d$d0[1:0]

  lr = nrow(tab)
  t=prop.table(addmargins(table(test, ref),2),2)
  p0=p1=.5
  R = cbind(r0 = p0*t[,'0']/(p0*t[,'0']+p1*t[,'1']),
            r1 = p1*t[,'1']/(p1*t[,'1']+p0*t[,'0']))
  temp = data.frame(cbind(tab, t, R))
  temp$ll = temp$Sum*temp$r1
  temp$lu = temp$Sum*abs(temp$r1-temp$r0)
  temp$lh = temp$Sum*temp$r0
  temp$lmin=pmin(temp$ll, temp$lu, temp$lh)
  rank(temp$r0)
  seq_along(temp$r0)

  Loss1 = function(A) {
    lr = nrow(t)
    sum(t[-(A:lr), 'Sum'] * (1 - R[-(A:lr), 'r0'])) +  # gewogen kans dat testwaarden NIET uit D0 komen
      sum(t[A:lr, 'Sum'] * (1 - R[A:lr, 'r1']))   # kans dat testwaarden NIET uit D1 komen
  }
  Loss.Y = Loss1(Yt)

  # is = Yt # is = get.intersection(ref = ref, test = test, adjust=adjust)
  is = get.intersection(ref = ref, test = test, adjust=adjust, ...) # closest test value
  isr = round(is[length(is)]) # tail has highest density
  wr = which(d$test==isr) # row number

  dl=nrow(d) # upperbound=lowerbound=5; tab[lowerbound:(wr),]
  # ifelse(lowerbound < (wr-1), tab[lowerbound:(wr-1),], NA)/tab[lowerbound:upperbound,]
  upperbound = c(wr:dl) # row numbers!
  lowerbound = c(wr:1)
  uncertain.interval= as.data.frame(expand.grid(lowerbound,upperbound))
  colnames(uncertain.interval) <- c('lowerbound', 'upperbound')

  fUI.Sp = function (x){ifelse(x[1] < wr, sum(tab[x[1]:(wr-1),1]), NA)/sum(tab[x[1]:x[2],1])}
  fUI.Se = function (x){sum(tab[wr:x[2],2])/sum(tab[x[1]:x[2],2])}
  fUI.ratio = function (x){(sum(tab[x[1]:x[2],2])/totpos)/(sum(tab[x[1]:x[2],1])/totneg)}
  uncertain.interval$UI.Sp = apply(uncertain.interval, 1, fUI.Sp) # x = unlist(uncertain.interval[1,])
  uncertain.interval$UI.Se = apply(uncertain.interval, 1, fUI.Se)
  uncertain.interval$UI.ratio = apply(uncertain.interval, 1, fUI.ratio)
  # which.min(abs(uncertain.interval$UI.ratio-1))


  f.UI.n = function(x){ sum(d$tot[x[1]:x[2]]) }
  uncertain.interval$UI.n = apply(uncertain.interval, 1, f.UI.n)

  f.Sp = function(x) {
    sum(d$d0[-(x[1]:dl)])/sum(d$d0[-(x[1]:x[2])]) }
  f.Se = function(x) {
    sum(d$d1[-(1:x[2])])/sum(d$d1[-(x[1]:x[2])]) }
  f.Acc = function(x) {
    (sum(d$d0[-(x[1]:dl)])+sum(d$d1[-(1:x[2])]))/
      sum(d$d0[-(x[1]:x[2])]+d$d1[-(x[1]:x[2])]) }
  # f.Loss = function(x) {
  #   sum(abs(d$d0[x[1]:x[2]]/totneg-d$d1[x[1]:x[2]]/totpos)) + # sum(abs(t[3:5,1]-t[3:5,2]))
  #    sum(d$d0[-(1:x[2])])/sum(d$d0) +  # FN # sum(t[6:10,1])
  #       +sum(d$d1[-(x[1]:dl)])/sum(d$d1) }  # FP # sum(t[1:2, 2])
  f.Loss = function(x) {
    sum(t[-(x[1]:lr), 'Sum'] * (1 - R[-(x[1]:lr), 'r0'])) +  # gewogen kans dat testwaarden NIET uit D0 komen
      sum(t[x[1]:x[2], 'Sum'] * abs(R[x[1]:x[2], 'r1'] - R[x[1]:x[2], 'r0'])) + # gewogen kans dat testwaarden niet verschillen
      sum(t[-(1:x[2]), 'Sum'] * (1 - R[-(1:x[2]), 'r1']))   # kans dat testwaarden NIET uit D1 komen
  }
  f.MCI.n = function(x){ sum(d$tot[-(x[1]:x[2])]) }

  uncertain.interval$MCI.Sp = apply(uncertain.interval, 1, f.Sp)
  uncertain.interval$MCI.Se = apply(uncertain.interval, 1, f.Se)
  uncertain.interval$MCI.Acc = apply(uncertain.interval, 1, f.Acc) # x = unlist(uncertain.interval[1,])
  uncertain.interval$MCI.n = apply(uncertain.interval, 1, f.MCI.n) # sum(d$tot)
  uncertain.interval$Loss = apply(uncertain.interval, 1, f.Loss) # x = uncertain.interval[1,]

  uncertain.interval$Random.gain=((round(uncertain.interval$MCI.Acc*uncertain.interval$MCI.n)+
                                 .5*uncertain.interval$UI.n)/ (totpos+totneg)) - Acc.Y

  # replace row numbers by test values
  uncertain.interval$lowerbound =  d$test[uncertain.interval$lowerbound]
  uncertain.interval$upperbound =  d$test[uncertain.interval$upperbound]

  sel = uncertain.interval$Loss <= Loss.Y
  uncertain.interval2 = uncertain.interval[sel,]
  uncertain.interval2 = uncertain.interval2[order(uncertain.interval2$Loss),]
  # sel = min(uncertain.interval$Total.Loss)==uncertain.interval$Total.Loss
  #sel = uncertain.interval2$UI.ratio > .45/.55 & uncertain.interval2$UI.ratio < .55/.45
  #uncertain.interval2 = uncertain.interval2[sel,]

  return(
    list(Youden=c(max.Youden=d$Y[Yt], # max Youden index
    threshold=d$test[Yt], # cut-point test values
    Sp=1-d$fpr[Yt], Se=d$tpr[Yt],Acc=Acc.Y, Loss=Loss.Y),
    data.table = d,
    intersection=is,
    Uncertain.Interval=uncertain.interval,
    candidates=uncertain.interval2[sel,])
    )
}


