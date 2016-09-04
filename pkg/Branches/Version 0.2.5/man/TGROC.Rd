\name{TGROC}
\alias{TGROC}
\alias{print.TGROC}
\alias{plot.TGROC}
\title{TG-ROC - Two Graphic Receiver Operating Characteristic}
\description{
TGROC draws a graph of sensitivity and specificity with the variations of a diagnostic test scale. Also, it demonstrates which cut-offs (or decision thresholds) may trichotomize the test results into a range where the test is good to identify those with the target condition, a inconclusive range and a range where the test is good to identify those without the target condition according with the researcher tolerance. Also, it estimates and graphically demonstrates good cut-offs by different methods. TGROC estimates non-parametric statistics and uses the AMORE package to simulate the parametric curve and values with a neural network.
}
\usage{
TGROC(gold,
      test,
      Cost=1,
      CL=0.95,
      Inconclusive=0.95,
      Prevalence=0,
      t.max=NULL,
      t.min=NULL,
      precision=.0001,  
      n.neurons=c(1,5,1),
      learning.rate.global=1e-2,
      momentum.global=0.3,
      error.criterium="LMS",
      Stao=NA,
      hidden.layer="sigmoid",
      output.layer="sigmoid",
      method="ADAPTgdwm",
      report=FALSE,
      show.step=5000,
      n.shows=1,
      Plot="Both",
      Plot.inc.range=TRUE,
      Plot.Cl=FALSE,
      Plot.cutoff="None",
      cex=0.5,
      cex.sub=0.85,
      Print=TRUE)
\method{plot}{TGROC}(x,...,
      Plot="Both",
      Plot.inc.range=TRUE,
      Plot.Cl=FALSE,
      Plot.cutoff="None",
      cex=0.5,
      cex.sub=0.85)            
}
\arguments{
  \item{gold}{The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must have two levels: must be coded either as 0 - without target condition - or 1 - with the target condition;  or could be coded \link[base]{as.factor} with the words "negative" - without target condition - and "positive"  - with the target condition.}
  \item{test}{ The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale. It may also work with discrete ordinal scale.}
  \item{Cost}{Cost = cost(FN)/cost(FP). MCT (misclassification cost term) will be used to estimate a good cut-off. It is a value in a range from 0 to infinite. Could be financial cost or a health outcome with the perception that FN are more undesirable than FP (or the other way around). This item will run into MCT - (1-prevalence)*(1-Sp)+Cost*prevalence(1-Se). Cost=1 means FN and FP have even cost. Cost = 0.9 means FP are 10 percent more costly. Cost = 0.769 means that FP are 30 percent more costly. Cost = 0.555 means that FP are 80 percent more costly. Cost = 0.3 means that FP are 3 times more costly. Cost = 0.2 means that FP are 5 times more costly. Also, it can be inserted as any ratio such as 1/2.5 or 1/4.}
  \item{CL}{Confidence limit. The limits of the confidence interval. Must be coded as number in range from 0 to 1. Default value is 0.95}
  \item{Inconclusive}{Inconclusive is a value that ranges from 0 to 1 that will identify the test range where the performance of the test is not acceptable and thus considered inconclusive. It represents the researcher tolerance of how good the test should be. If it is set to 0.95 (which is the default value), test results that have less than 0.95 sensitivity and specificity will be in the inconclusive range.}
  \item{Prevalence}{Prevalence of the disease in the population who the test will be performed. It must be a value from 0 to 1. If left 0 (the default value), this will be replaced by the disease prevalence in the sample. This value will be used in the MCT and Efficiency formulas to estimate good cut-offs.}
  \item{t.max}{Test upper range limit to be set as numeric. It will be used to simulate the parametric curve. If left NULL TGROC will assume that the sample maximum value is the upper limit of the test.}
  \item{t.min}{Test lower range limit to be set as numeric. It will be used to simulate the parametric curve. If left NULL TGROC will assume that the sample minimum value is the lower limit of the test.}
  \item{precision}{The test precision is the unit of variation of the test and should be set as numeric. It will be used to simulate the parametric curve. It will express how many estimations the network will do between each test unit. It is interesting the precision to be something between 1/2 to 1/10 of the test unit. The higher the precision, smoother the parametric curve will look. However, if too much precision is set the function may give an error as a result. This also may depends on the amount of observations in the dataset}
  \item{n.neurons}{Numeric vector containing the number of neurons of each layer. See \link[AMORE]{newff}.}
  \item{learning.rate.global}{Learning rate at which every neuron is trained. See \link[AMORE]{newff}.}
  \item{momentum.global}{Momentum for every neuron. See \link[AMORE]{newff}.}
  \item{error.criterium}{Criteria used to measure to proximity of the neural network prediction to its target. See \link[AMORE]{newff}.}
  \item{Stao}{Stao parameter for the TAO error criteria. See \link[AMORE]{newff}.}
  \item{hidden.layer}{Activation function of the hidden layer neurons. See \link[AMORE]{newff}.}
  \item{output.layer}{Activation function of the hidden layer neurons. See \link[AMORE]{newff}.}
  \item{method}{Preferred training method. See \link[AMORE]{newff}.}
  \item{report}{Logical value indicating whether the training function should keep quiet. See \link[AMORE]{train}.}
  \item{show.step}{Number of epochs to train non-stop until the training function is allow to report. See \link[AMORE]{train}.}
  \item{n.shows}{Number of times to report (if report is TRUE).See \link[AMORE]{train}.}
  \item{Plot}{Possible values are: "None", "Both", "Parametric" and "Non-parametric". TGROC may plot parametric, non-parametric, both or no plot at all depending of this option. Default is to plot both curves.}
  \item{Plot.inc.range}{Plot inconclusive range. If TRUE, the lines representing the limits of the inconclusive range will be displayed. Default is TRUE. Parametric inconclusive range will be displayed if Plot = "Parametric" or Plot = "Both" and non-parametric inconclusive range otherwise. If Plot is FALSE then Plot.inc.range is not considered.}
  \item{Plot.Cl}{Plot confidence limits. If TRUE, confidence bands for sensitivity and specificity curves will be displayed. If Plot = "Parametric" or Plot = "Both" than parametric bands are displayed and non-paramentric otherwise. Default is FLASE. If Plot is FALSE than Plot.Cl is not considered.}
  \item{Plot.cutoff}{A line representing the estimated best cut-off (threshold) will be displayed. If Plot is FALSE then Plot.cutoff is not considered. If Plot = "Parametric" or Plot = "Both" than the parametric values are represented and non-parametric otherwise. Default is "None". Possible values are:
        "Se=Sp" - the cut-off which Sensitivity is equal to Specificity;

        "Max.Efficiency" - the cut-off which maximize the efficiency;

        "Min.MCT" - the cut-off which minimizes the misclassification cost term.}
  \item{cex}{See \link[graphics]{par}. A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.}
  \item{cex.sub}{See \link[graphics]{par}. Controls the font size in the subtitle. If Plot is FALSE than cex.sub is not considered.}
  \item{Print}{If FALSE, statistics estimated by TGROC will not be displayed in the output window. Default is TRUE.}
  \item{x}{For the plot function, x is an object created TGROC function.}  
  \item{...}{Other plot parameters form \link[graphics]{plot.default}}
}
\details{
There are two main advantages of TG-ROC over ROC analysis: (1) for the uninitiated is much easier to understand how sensitivity and specificity changes with different cut-offs; (2) and because of the graphical display is much easier to understand and estimate reasonable inconclusive test ranges. Occasionally the MCT or Efficiency cut-offs may be set outside the inconclusive range. This may happens with extreme values of Cost and population prevalence. If this is the case, perhaps the inconclusive range may not be of interest or not applicable. Also, if the test is too good or inconclusive range tolerance is set too low, then there may be no inconclusive range at all, because sensitivity and specificity may not be below this tolerance at the same time. If this is the case, setting a higher inconclusive tolerance may work. Tests results matching the cut-off values will be considered a positive test.  TGROC assumes that subjects with higher values of the test are with the target condition and those with lower values are without the target condition. Tests that behave like glucose (middle values are supposed to be normal and extreme values are supposed to be abnormal) and immunefluorescence (lower values - higher dilutions - are suppose to be abnormal) will not be correctly analyzed. In the latter, multiplying the test results by -1 or other transformation before analysis could make it work. The validity measures such as Sensitivity, Specificity and Likelihood ratios and its confidence limits are estimated as  in \code{\link{diagnosis}} function. MCT and Efficiency are estimated as in \code{\link{ROC}} function. Non-parametric confidence bands are estimated by binomial confidence interval and parametric with normal confidence interval. The parametric curve and validity measures are estimated with a neural network strategy using the AMORE package. Usually neural networks, uses a subset of the data to estimate weights, a subset to calibrate/validate the weights and a third subset to simulate the function. TGROC uses only the estimate and simulate steps, therefore there is no stopping rule for the neural network parametric estimation. The only way to check the fit of the neural network is to visually compare with the non-parametric curve. If the curve looks weird or not good enough, than progressive slight changes in the momentum, learning rate, number of layers and / or other parameters should work fine.}
\value{
  \item{Sample size}{Amount of subjects analyzed.}
  \item{Sample prevalence}{Prevalence of target condition in the sample.}
  \item{Population prevalence.}{Informed prevalence in the population.}
  \item{Test summary}{A summary of central and dispersion tendencies of test results.}
  \item{Non-parametric inconclusive limits.}{Estimate of the inconclusive limits of the tests and its corresponding validity measures.}
  \item{Non-parametric best cut-offs.}{The cut-offs estimated by different methods and its corresponding validity measures.}
  \item{Parametric inconclusive limits.}{Estimate of the inconclusive limits of the tests and its corresponding validity measures with the parametric simulation.}
  \item{Parametric best cut-off}{The cut-offs estimated by different methods and its corresponding validity measures with the parametric simulation.}
}
\references{Greiner, M. (1996) Two-graph receiver operating characteristic (TG-ROC): update version supports optimization of cut-off values that minimize overall misclassification costs. J.Immunol.Methods 191:93-94.

M. Greiner (1995) Two-graph receiver operating characteristic (TG-ROC): a Microsoft-EXCEL template for the selection of cut-off values in diagnostic tests. Journal of Immunological Methods. 185(1):145-146.

M. Greiner, D. Sohr, P. Gobel (1995) A modified ROC analysis for the selection of cut-off values and the definition of intermediate results of serodiagnostic tests. Journal of immunological methods. 185(1):123-132.
}
\author{Pedro Brasil; - \email{diagnosismed-list@lists.r-forge.r-project.org}}
\note{Bug reports, malfunctioning, or suggestions for further improvements or contributions can be sent, preferentially, through the DiagnosisMed email list, or R-Forge website \url{https://r-forge.r-project.org/projects/diagnosismed/}.
}
\seealso{\code{\link{interact.ROC}},\code{\link{ROC}},\code{\link{diagnosis}},\link[ROCR]{performance},\link[epitools]{binom.conf.int},\code{ nonbinROC::performance}}
\examples{
# Loading a dataset.
data(tutorial)
# Attaching dataset
attach(tutorial)
# Running the analysis
TGROC(gold=Gold,test=Test_B)
rm(tutorial)
}
\keyword{iplot}
\keyword{univar}
\keyword{htest}
