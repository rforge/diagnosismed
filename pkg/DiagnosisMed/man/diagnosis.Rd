\name{diagnosis}
\alias{diagnosis}
\alias{diagnosisI}
\alias{print.diag}
\title{Diagnostic test accuracy evaluation}
\description{diagnosis estimate sensitivity, specificity, predictive values, likelihood ratios, area under ROC curve and other validity measures for binary diagnostic test evaluation.
}
\usage{
diagnosis(gold, test, CL = 0.95, print = TRUE, plot = FALSE)
diagnosisI(TP,FN,FP,TN, CL = 0.95, print = TRUE, plot = FALSE)
}
\arguments{
  \item{gold}{The reference test (gold standard)}
  \item{test}{The index test (test under evaluation)}
  \item{TP}{A number representing True Positives from a 2x2 table}
  \item{FN}{A number representing False Negatives from a 2x2 table}
  \item{FP}{A number representing False Positives from a 2x2 table}
  \item{TN}{A number representing True Negatives from a 2x2 table}
  \item{print}{If TRUE, diagnosis will print in the output window the statistics resulted from the 2x2 table. Default is TRUE.}
  \item{plot}{If TRUE, diagnosis will plot a ROC curve of the test under evaluation. This plot may later be edited, as any other plot, with title, legends etc. Default is FALSE.}
  \item{CL}{Confidence limits for confidence intervals. Must be a numeric value between 0 and 1. Default is 0.95.}
}
\details{
In diagnosis, the values entered must be two variables in a data frame. The first one should be the gold standard and the second should be the index test. These two variables must be coded either as numeric - 0 for negative and 1 for a positive test - or with the words "positive" and "negative" or "presence" and "absence". Never code the variables as TRUE or FALSE; diagnosis will level this categories incorrectly. diagnosisI executes the very same statistics, but here, the values of a 2x2 table are inputted. TP is true positive; TN is true negative; FP is false positive and FN is false negative. Sensitivity, Specificity, Predictive values and Accuracy confidence limits rely on binomial distribution, which does not give result outside [0:1] such as normal distribution or asymptotic theory. DOR, Likelihood  ratios and Youden J index confidence limits rely on normal approximation (Wald method for likelihoods). The AUC (area under the ROC curve) is estimated by trapezoidal (or Mann-Witney method). See formula below.
}
\value{
  A 2x2 table from which the validity measures are calculated.
  \item{Sample size}{The number of subjects analyzed.}
  \item{Prevalence}{The proportion classified as with the target condition by the
        reference standard}
  \item{Sensitivity}{the probability of the test to correctly classify subjects
        with the target condition (TP/(TP+FN))}
  \item{Specificity}{the probability of the test to correctly classify subjects
        without the target condition (TN/(TN+FP))}
  \item{Predictive values}{the probabilities of being with (positive predictive
        value) (TP/(TP+FP)) or without (negative predictive value) the target
         condition given a test result (TN/(TN+FN)).}
  \item{Likelihood ratios}{the probability of test a result in people with the
       target condition, divided by the probability of the same test result in
       people without the target condition (PLR = Se/(1-Sp); NLR = (1-Sp)/Se)}
  \item{Diagnostic odds ratio}{represents the overall discrimination of a
        dichotomous test, and is equivalent to the ratio of PLR and NLR.}
  \item{Error trade off}{Is the amount of false positives traded with false negatives
        for each decision threshold; here expressed as an odd - for binary results
        there is only one threshold}
  \item{Error rate}{Expresses how many errors we make when we diagnose patients with an abnormal test result as diseased, and those with a normal test result as non-diseased ((FP+FN)/sample size).}
  \item{Accuracy}{overall measure that express the capacity of the test to correctly classify subjects with and without the target condition ((TP+TN)/(sample size))}
  \item{Area under ROC curve}{overall measure of accuracy - here the method is the trapezoidal (same as Mann-Whitney rank sum test or DeLong method). It gives identical results as (Se+SP)/2.}
  }
\references{

Knotterus. The Evidence Based Clinical Diagnosis; BMJBooks, 2002.

Xiou-Hua Zhou, Nancy A Obuchowsky, Donna McClish. Statistical Mehods in diagnostic Medicine; Wiley, 2002.

Simel D, Samsa G, Matchar D (1991). Likelihood ratios with confidence: Sample size estimation for diagnostic test studies. Journal of Clinical Epidemiology 44: 763 - 770
}
\author{Pedro Brasil - \email{diagnosismed-list@lists.r-forge.r-project.org}}
\note{Bug reports, malfunctioning, or suggestions for further improvements or contributions can be sent, preferentially, through the DiagnosisMed email list, or R-Forge website \url{https://r-forge.r-project.org/projects/diagnosismed/}.
}
\seealso{\code{\link{plot.diag}}, \code{\link{ROC}},\code{\link{LRgraph}},\link[epitools]{binom.conf.int},
 \link[epibasix]{sensSpec},\link[epiR]{epi.tests}}
\examples{
# Simulating a dataset
mydata<-as.data.frame(rbind(c("positive","positive"),c("positive","negative"),
        c("negative","positive"),c("negative","negative")))
# Auxiliary vector to indicate the number of times to expand the lines
Freq<-as.numeric(c(173,13,27,638))
# Binding the dataset to the auxiliary vector
mydata<-as.data.frame(cbind(mydata,Freq))
# Naming the columns of the data set
colnames(mydata)<-c("Gold","mytest","Freq")
# A little description of the data set to check if it is ok!
str(mydata)
# Removing unwanted objects
rm(Freq)
# Expanding the four lines dataset to an adequate format
mydata<-expand(mydata, index.var = "Freq")
# Attaching the data set
attach(mydata)
# Running the diagnosis analysis
test<-diagnosis(Gold,mytest,print=FALSE)
print(test)
# The same as
diagnosis(Gold,mytest)
# Draw a nomogram from this test
test<-diagnosis(Gold,mytest,print=FALSE)
plot(test)
#Different from - draw a ROC curve from this test
diagnosis(Gold,mytest,plot=TRUE,print=FALSE)
#Inserting values from a 2x2 table
diagnosisI(364,22,17,211)
rm(test)
rm(mydata)
}
\keyword{univar}
\keyword{htest}

