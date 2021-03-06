% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TGROC.R, R/plot.TGROC.R, R/print.TGROC.R
\name{TGROC}
\alias{TGROC}
\alias{plot.TGROC}
\alias{print.TGROC}
\title{TGROC - Two Graphic Receiver Operating Characteristic}
\usage{
TGROC(ref, test, Cost = 1, CL = 0.95, Inconclusive = 0.95,
  Prevalence = NULL, pop.prevalence = NULL, t.max = NULL, t.min = NULL,
  precision = 0.05, reverse = "auto", n.neurons = c(1, 5, 1),
  learning.rate.global = 0.01, momentum.global = 0.3,
  error.criterium = "LMS", Stao = NA, hidden.layer = "sigmoid",
  output.layer = "sigmoid", method = "ADAPTgdwm", report = FALSE,
  show.step = 5000, n.shows = 1)

\method{plot}{TGROC}(x, ..., Plot.type = c("TGROC", "ROC", "None"),
  Plot = c("Binormal", "Non-parametric", "NN-parametric"),
  Plot.inc.area = TRUE, Plot.CL = FALSE, Plot.threshold = "None",
  threshold.arg = list(col = gray(0.5), lty = 6),
  ylab = "Sensitivity & Specificity", xlab = "Test scale", ylim = c(0, 1),
  xlim = "auto", auto.title = TRUE, title.arg = list(cex.sub = 0.65, line
  = 4.5), shade.args = list(col = gray(0.8), density = 45, border = NA),
  np.Se.args = list(type = "o", col = "blue", lty = 1, cex = 0.5),
  np.Se.ci.args = list(lty = 5, col = "blue"), np.Sp.args = list(type = "o",
  col = "red", lty = 2, cex = 0.5), np.Sp.ci.args = list(lty = 3, col =
  "red"), NN.Se.args = list(type = "l", col = "blue", lty = 1, cex = 0.5),
  NN.Se.ci.args = list(lty = 5, col = "blue"), NN.Sp.args = list(type = "l",
  col = "red", lty = 2, cex = 0.5), NN.Sp.ci.args = list(lty = 3, col =
  "red"), BN.Se.args = list(type = "l", col = "blue", lty = 1, cex = 0.5),
  BN.Se.ci.args = list(lty = 5, col = "blue"), BN.Sp.args = list(type = "l",
  col = "red", lty = 2, cex = 0.5), BN.Sp.ci.args = list(lty = 3, col =
  "red"), roc.np.line.args = list(col = "red", lwd = 1, type = "o", cex =
  0.5), roc.np.point.args = list(pch = 21, col = par()$bg, bg = "red", cex =
  2), roc.NN.line.args = list(col = "navyblue", lwd = 3),
  roc.NN.point.args = list(pch = 21, col = par()$bg, bg = "navyblue", cex =
  2), roc.BN.line.args = list(col = "darkgreen", lwd = 3),
  roc.BN.point.args = list(pch = 21, col = par()$bg, bg = "darkgreen", cex =
  2), roc.xlab = "1 - Specificity", roc.ylab = "Sensitivity", grid = TRUE,
  auto.legend = TRUE, legend.args = list(x = "top", border = NA, bty = "n",
  xpd = NA, inset = -0.18, ncol = 2, cex = 0.8))

\method{print}{TGROC}(x, ..., digits = 3, print.type = "all")
}
\arguments{
\item{ref}{The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)}

\item{test}{The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.}

\item{Cost}{Represents weights of wrong classifications as cost(FN)/cost(FP). See \code{\link{thresholds}}.}

\item{CL}{Confidence limit. The limits of the confidence interval. Must be coded as number in range from 0 to 1. Default value is 0.95.}

\item{Inconclusive}{Inconclusive is a value that ranges from 0 to 1 that will identify the test range where the performance of the test is not acceptable and thus considered inconclusive. It represents the researcher tolerance of how good the test should be. If it is set to 0.95 (which is the default value), test results that have less than 0.95 sensitivity and specificity will be in the inconclusive range.}

\item{Prevalence, pop.prevalence}{Population condition prevalence. The value of pop.prevalence will be passed to \code{\link{thresholds}} and The value of Prevalence will be passed to \code{\link{SS}}.}

\item{t.min, t.max, precision}{Test minimum, maximum and intervals to simulate the parametric estimation as \code{\link[base]{seq}}(t.min, t.max, precision). See \code{\link{SS}}.}

\item{reverse}{It accepts 'auto', TRUE or FALSE. ROC analysis assumes that higher values are from subjects with the condition and lower values are from subjects without the condition. If it occurs the other way around, the ROC analysis and its interpretation must be reversed. See \code{\link{np.auROCc}} and \code{\link{SS}}.}

\item{n.neurons}{Numeric vector containing the number of neurons of each layer. See \code{\link[AMORE]{newff}}.}

\item{learning.rate.global}{Learning rate at which every neuron is trained. See \code{\link[AMORE]{newff}}.}

\item{momentum.global}{Momentum for every neuron. See \code{\link[AMORE]{newff}}.}

\item{error.criterium}{Criteria used to measure to proximity of the neural network prediction to its target. See \code{\link[AMORE]{newff}}.}

\item{Stao}{Stao parameter for the TAO error criteria. See \code{\link[AMORE]{newff}}.}

\item{hidden.layer}{Activation function of the hidden layer neurons. See \code{\link[AMORE]{newff}}.}

\item{output.layer}{Activation function of the hidden layer neurons. See \code{\link[AMORE]{newff}}.}

\item{method}{Preferred training method. See \code{\link[AMORE]{newff}}.}

\item{report}{Logical value indicating whether the training function should keep quiet. See \code{\link[AMORE]{train}}.}

\item{show.step}{Number of epochs to train non-stop until the training function is allow to report. See \code{\link[AMORE]{train}}.}

\item{n.shows}{Number of times to report (if report is TRUE). See \code{\link[AMORE]{train}}.}

\item{x}{The output of \code{TGROC} function for printing and ploting.}

\item{...}{Additional arguments passed either to \code{\link[base]{print}} or \code{\link[graphics]{plot.default}}}

\item{Plot.type}{Valid options are "TGROC" and "ROC". This will define the type of plot to return.}

\item{Plot}{Valid options are "Binormal", "Non-parametric", and "NN-parametric" (for neural network). See \code{\link{SS}}. In the TGROC plot, this will define which sensitivity and specificity lines will be plotted. It may be one, two or all of them (although it become a polluted graph). The order they are set matters, as some options, e.g. Inconclusive argument, will pick the first option to show in the plot.}

\item{Plot.inc.area}{Acceptable options are TRUE (default) or FALSE. If TRUE and Plot.type = "TGROC", it will display a gray shade in the plot (using the \code{\link[graphics]{polygon}}) representing the inconclusive area defined by the first method in the Plot argument (e.g. "Binormal").}

\item{Plot.CL}{Acceptable options are TRUE or FALSE (default). If TRUE and Plot.type = "TGROC", it will pick the first method in the Plot argument (e.g. "Binormal"), and display lines in the plot representing a confidence band for this method.}

\item{Plot.threshold}{Method of the decision threshold to be represented as a vertical line in the plot. Acceptable values are "None" (default), "Max Youden", "Max Accuracy", "Max Accuracy area", "Max DOR", "Min Error rate", "Min ROC distance", "Se = Sp", "Min MCT", "Max Efficiency". See \code{\link{thresholds}}.}

\item{threshold.arg}{If Plot.type = "TGROC" and Plot.threshold is any value different from "None", this list of arguments with graphical parameters will be passed to \code{\link[graphics]{abline}} to plot a vertical line represented the chosen threshold. Internally, the \code{\link[graphics]{abline}} v argument will be replaced by the test value from Plot.threshold argument.}

\item{ylab, xlab, ylim, xlim, roc.xlab, roc.ylab}{Values to be passed to \code{\link[graphics]{plot.default}} for labels and limits of the axis. If Plot.type = "TGROC", then former will be used for axis labels, if Plot.type = "ROC", then the later will be use for axis labels. xlim also accepts the value "auto". If this is the case, the test range from the SS analysis of the first option of the Plot argument will be used as horizontal limits of the plot. In essence, it will be defined by t.min and t.max.}

\item{auto.title}{Logical. If TRUE, a subtitle will be automatically be placed at the bottom with the values of the inconclusive limits of the chosen threshold.}

\item{title.arg}{A list of arguments to be passed to \code{\link[graphics]{title}}. TGROC automatically generates a subtitle with values plotted. Internally the argument sub is replaced by values from other options. Sometimes tiny adjustments are required.}

\item{shade.args}{If Plot.type = "TGROC" and Plot.inc.area = TRUE, then this list of arguments will be passed to \code{\link[graphics]{polygon}} and plot a shade representing the inconclusive area. Internally, the coordinates will be extracted from in inconclusive limits from the first method defined ate the Plot argument.}

\item{np.Se.args}{If Plot.type = "TGROC" and Plot argument is "Non-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the sensitivity line. Internally, the x and y arguments of this list will be replaced by the non-parametric sensitivity values.}

\item{np.Se.ci.args}{If Plot.type = "TGROC" and Plot first argument is "Non-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the non-parametric sensitivity confidence band. Internally, the x and y arguments of this list will be replaced by the non-parametric sensitivity values.}

\item{np.Sp.args}{If Plot.type = "TGROC" and Plot argument is "Non-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the specificity line. Internally, the x and y arguments of this list will be replaced by the non-parametric specificity values.}

\item{np.Sp.ci.args}{If Plot.type = "TGROC" and Plot first argument is "Non-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the non-parametric specificity confidence band. Internally, the x and y arguments of this list will be replaced by the non-parametric specificity values.}

\item{NN.Se.args}{If Plot.type = "TGROC" and Plot argument is "NN-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the sensitivity line. Internally, the x and y arguments of this list will be replaced by the NN sensitivity values.}

\item{NN.Se.ci.args}{If Plot.type = "TGROC" and Plot first argument is "NN-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the NN-parametric sensitivity confidence band. Internally, the x and y arguments of this list will be replaced by the NN-parametric sensitivity values.}

\item{NN.Sp.args}{If Plot.type = "TGROC" and Plot argument is "NN-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the specificity line. Internally, the x and y arguments of this list will be replaced by the NN specificity values.}

\item{NN.Sp.ci.args}{If Plot.type = "TGROC" and Plot first argument is "NN-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the NN-parametric specificity confidence band. Internally, the x and y arguments of this list will be replaced by the NN-parametric specificity values.}

\item{BN.Se.args}{If Plot.type = "TGROC" and Plot argument is "Binormal", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the sensitivity line. Internally, the x and y arguments of this list will be replaced by the binormal sensitivity values.}

\item{BN.Se.ci.args}{If Plot.type = "TGROC" and Plot first argument is "Binormal" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the binormal sensitivity confidence band. Internally, the x and y arguments of this list will be replaced by the binormal sensitivity values.}

\item{BN.Sp.args}{If Plot.type = "TGROC" and Plot argument is "Binormal", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the specificity line. Internally, the x and y arguments of this list will be replaced by the binormal specificity values.}

\item{BN.Sp.ci.args}{If Plot.type = "TGROC" and Plot first argument is "Binormal" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the binormal specificity confidence band. Internally, the x and y arguments of this list will be replaced by the binormal specificity values.}

\item{roc.np.line.args}{If Plot.type = "ROC" and Plot argument is "Non-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the empirical ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical sensitivity and specificity.}

\item{roc.np.point.args}{If Plot.type = "ROC" and Plot argument is "Non-parametric" and Plot.threshold is different from "None", then this list of arguments will be passed to \code{\link[graphics]{points}} to plot a point representing the desired threshold on the empirical ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical point sensitivity and point specificity of the desired threshold.}

\item{roc.NN.line.args}{If Plot.type = "ROC" and Plot argument is "NN-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the NN ROC curve. Internally, the x and y arguments of this list will be replaced by the NN sensitivity and specificity.}

\item{roc.NN.point.args}{If Plot.type = "ROC" and Plot argument is "NN-parametric" and Plot.threshold is different from "None", then this list of arguments will be passed to \code{\link[graphics]{points}} to plot a point representing the desired threshold on the NN ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical point sensitivity and point specificity of the desired threshold.}

\item{roc.BN.line.args}{If Plot.type = "ROC" and Plot argument is "Binormal", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the binormal ROC curve. Internally, the x and y arguments of this list will be replaced by the binormal sensitivity and specificity.}

\item{roc.BN.point.args}{If Plot.type = "ROC" and Plot argument is "Binormal" and Plot.threshold is different from "None", then this list of arguments will be passed to \code{\link[graphics]{points}} to plot a point representing the desired threshold on the binormal ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical point sensitivity and point specificity of the desired threshold.}

\item{grid}{If Plot.type = "ROC" and grid = TRUE, then a grid will be plotted. Set to FALSE to avoid ploting the grid.}

\item{auto.legend}{Logical. If TRUE, then a legend will be automatically place in the plot. If Plot.type = "ROC", then the auto.legend will work only if Plot.threshold is different from "None".}

\item{legend.args}{If auto.legend = TRUE, then this list of arguments will be passed to \code{\link[graphics]{legend}}. Internally, several arguments from lines and shade lists of arguments will be passed to \code{\link[graphics]{legend}}. If Plot.type = "TGROC", then auto.legend will use the BN.Se.args and BN.Sp.args to extract the line type and line colors. If Plot.type = "ROC", the auto.legend is way less flexible, e.g. the positions are fixed for all methods.}

\item{digits}{Number of digits passed to \code{\link[base]{print}}.}

\item{print.type}{To print all the results (default) or just sections of the result. Options are "all", "NN-parametric", "Non-parametric", "Binormal". A combination of them is also possible.}
}
\value{
A list of class TGROC with the following:
\itemize{
  \item \code{sample.size} The sample size.
  \item \code{sample.prevalence} The sample prevalence.
  \item \code{pop.prevalence} The population prevalence.
  \item \code{cost} The informed cost.
  \item \code{test.summary} A data.frame statistics such as mean and SD from the test.
  \item \code{inc} The desired inconclusive limit.
  \item \code{conf.limit} The desired confidence limit.
  \item \code{AUC} The estimated area under the ROC curve, SE and confidence limits.
  \item \code{SS} A data.frame with the empirical sensitivity and specificity trade-off (ROC) analysis.
  \item \code{np.inconclusive} A data.frame with the non-parametric inconclusive limits and their validity measures.
  \item \code{np.best.threshold} A data.frame with the non-parametric thresholds and their validity measures.
  \item \code{NN.SS} A data.frame with the neural network  sensitivity and specificity trade-off (ROC) analysis.
  \item \code{NN.inconclusive} A data.frame with the inconclusive limits from neural network smoothed ROC analysis and their validity measures.
  \item \code{NN.best.threshold} A data.frame with the NN thresholds and their validity measures.
  \item \code{BN.SS} A data.frame with the binormal sensitivity and specificity trade-off (ROC) analysis.
  \item \code{BN.inconclusive} A data.frame with the inconclusive limits from binormal ROC analysis and their validity measures.
  \item \code{BN.best.threshold} A data.frame with the binormal thresholds and their validity measures.
}
}
\description{
\code{TGROC} is a wrapper of several functions to conduct performance analysis of continuous scale diagnostic tests, including good decision thresholds by a variety of methods and intermediate inconclusive range. See \code{\link{thresholds}}. The inconclusive range thresholds threechotomizes the test results into a range where the test is good to identify those with the target condition, a inconclusive range and a range where the test is good to identify those without the target condition according to a minimum required sensitivity and specificity. TGROC estimates non-parametric ROC analysis, bi-normal parametric ROC analysis, and uses the \pkg{AMORE} package to simulate a robust parametric curve values with a neural network. See \code{\link{SS}}. It draws a graph of sensitivity and specificity with the variations of a diagnostic test scale. Optionally, it may display a ROC plot. Automatically, the plot function will display subtitle with the values of the thresholds from the chosen methods, and optionally legend.
}
\details{
There are two main advantages of TG-ROC over ROC analysis: (1) for the uninitiated is much easier to understand how sensitivity and specificity changes with different thresholds; (2) and because of the graphical display is much easier to understand and estimate reasonable inconclusive test ranges. Occasionally, thresholds may be set outside the inconclusive range. This may happens with extreme values of Cost and population prevalence, or very skewed or unusual test values distributions. If this is the case, perhaps the inconclusive range may not be of interest or not applicable. Also, if the test is too accurate or inconclusive range tolerance is set too low, then there may be no inconclusive range at all, because sensitivity and specificity may not be below this tolerance at the same time. If this is the case, setting a higher inconclusive tolerance may work. Tests results matching the threshold values will be considered a positive test. \code{TGROC} assumes that subjects with higher values of the test are with the target condition and those with lower values are without the target condition. Tests that behave like glucose (middle values are supposed to be normal and extreme values are supposed to be abnormal) will not be correctly analyzed. If lower test values are from subjects with the condition and higher values are from subject without the condition, the analysis and its interpretation must reversed. The validity measures such as Sensitivity, Specificity and Likelihood ratios and its confidence limits are estimated as  in \code{\link{diagnosis}} function. Non-parametric confidence bands are estimated with \code{\link{binom.CI}} and parametric with normal confidence interval. The parametric curve and validity measures are estimated with a neural network strategy using the \pkg{AMORE} package. Usually neural networks, uses a subset of the data to estimate weights, a subset to calibrate/validate the weights and a third subset to simulate the function. TGROC uses only the estimate and simulate steps, therefore there is no early stopping rule for the neural network parametric estimation to prevent overfitting. The only way to check the fit of the neural network is to visually compare with the non-parametric curve. If the curve looks weird or not good enough, than progressive slight changes in the momentum, learning rate, number of layers or other parameters should work fine. The AUC (area under the ROC curve) is estimated by the trapezoidal method (also known as Mann-Whitney statistic), its confidence interval is estimated by DeLong method. See \code{\link{np.auROCc}}. The AUC confidence limits should be used only to compare the AUC with the null value for AUC which is 0.5 and not to compare the AUC from different tests.
}
\examples{

data(rocdata)
x <- TGROC(ref = rocdata$Gold, test = rocdata$test2, precision = .005)

# Printing just the NN part of the analysis
print(x, print.type = "NN-parametric")

# Ploting just the empirical ROC curve with Max Youden threshold.
plot(x, Plot.type = "ROC", Plot = "Non-parametric", Plot.threshold = "Max Youden")

# Ploting just the NN ROC curve
# The plot does not reach the corner because there is no 100\% Sensitivity at any threshold.
plot(x, Plot.type = "ROC", Plot = "NN-parametric", Plot.threshold = "Max Youden")

# Ploting just the Binormal ROC curve
plot(x, Plot.type = "ROC", Plot = "Binormal", Plot.threshold = "Max Youden")

# Ploting overplotted curves. The Binormal curve is way different from the remaining.
plot(x, Plot.type = "ROC",
     Plot = c("Non-parametric","NN-parametric","Binormal"),
     Plot.threshold = "Max Youden")

# Ploting the TGROC curve overploting the Non-parametric and the NN- parametric.
plot(x, Plot = c("Non-parametric","NN-parametric"))

# TGROC curve overploting the Non-parametric (and confidence band) and the NN- parametric.
# If legend is not good enough, scoot up a little bit decreasing the inset.
plot(x, Plot = c("NN-parametric", "Non-parametric"), Plot.CL = TRUE,
     legend.args = list(x = "top", border = NA, bty = "n", xpd = NA,
                             inset = -.25, ncol = 2, cex = .8))

# Changing the order of the Plot argument changes the inconclusive area shade.
plot(x, Plot = c("NN-parametric","Non-parametric"))

# The Binormal fit looks worst than the NN fit
plot(x, Plot = c("Binormal","Non-parametric"))

data(tutorial)
# Running TGROC with a very accurate test.
# warnings concerning ties and concerning the inconclusive limits.
x <- TGROC(ref = ifelse(tutorial$Gold == "pos", 1, 0), test = tutorial$Test_A, precision = .005)

# All lower inconclusive limits are higher than the upper inconclusive limit
plot(x, Plot = c("NN-parametric","Non-parametric"))
plot(x, Plot = c("Binormal","Non-parametric"))
plot(x, Plot = c("Non-parametric"))

# Increasing the inconclusive solves the issue.
x <- TGROC(ref = ifelse(tutorial$Gold == "pos", 1, 0),
     test = tutorial$Test_A, precision = .005, Inconclusive = .99)

# All lower inconclusive limits are higher than the upper inconclusive limit
# Nevertheless, check how the inconclisve area chages as method Changes.
plot(x, Plot = c("NN-parametric","Non-parametric"))
plot(x, Plot = c("Binormal","Non-parametric"))
plot(x, Plot = c("Non-parametric"))

# Ploting with different thrensholds without the inconclusive area.
x <- TGROC(ref = ifelse(tutorial$Gold == "pos", 1, 0), test = tutorial$Test_B,
     precision = .005, Inconclusive = .90)
plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.inc.area = FALSE,
     Plot.threshold = "Se = Sp")
plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.inc.area = FALSE,
     Plot.threshold = "Min ROC distance")
plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.inc.area = FALSE,
     Plot.threshold = "Max DOR")

# Ploting the inconclusive area and the threshold simultaneously.
# The subtitle may get overplotted, expanding the margins and adjusting the line may help.
par(mar = c(6.1, 4.1, 4.1, 2.1))
plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.threshold = "Max Efficiency",
     title.arg = list(cex.sub = 0.7, line = 5))

# Back to original margins.
par(mar = c(5.1, 4.1, 4.1, 2.1))

rm(rocdata, tutorial, x)

}
\references{
Greiner, M. (1996) Two-graph receiver operating characteristic (TG-ROC): update version supports optimization of cut-off values that minimize overall misclassification costs. J.Immunol.Methods 191:93-94.

M. Greiner (1995) Two-graph receiver operating characteristic (TG-ROC): a Microsoft-EXCEL template for the selection of cut-off values in diagnostic tests. Journal of Immunological Methods. 185(1):145-146.

M. Greiner, D. Sohr, P. Gobel (1995) A modified ROC analysis for the selection of cut-off values and the definition of intermediate results of serodiagnostic tests. Journal of immunological methods. 185(1):123-132.
}
\seealso{
\code{\link{diagnosis}}, \code{\link{binom.CI}}, \code{\link{SS}}, \code{\link{thresholds}}, \code{\link{np.auROCc}}
}
