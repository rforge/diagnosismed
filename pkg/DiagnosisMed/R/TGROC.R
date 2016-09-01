#' TGROC - Two Graphic Receiver Operating Characteristic
#'
#' @name TGROC
#'
#' @description \code{TGROC} is a wrapper of several functions to conduct performance analysis of concitnuous scale diagnostic tests, including good decision thresholds by a variety of methods and intermediate inconclusive range. See \code{\link{thresholds}}. The inconclusive range thresholds threechotomizes the test results into a range where the test is good to identify those with the target condition, a inconclusive range and a range where the test is good to identify those without the target condition according to a minimum required sensitivity and specificity. TGROC estimates non-parametric ROC analysis, bi-normal parametric ROC analysis, and uses the \pkg{AMORE} package to simulate a robust parametric curve values with a neural network. See \code{\link{SS}}. It draws a graph of sensitivity and specificity with the variations of a diagnostic test scale. Optioanlly, it may display a ROC plot. Automatically, the plot function will display subtitle with the values of the tresholds from the chosen methods, and optionally legend.
#'
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#'
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#'
#' @param CL Confidence limit. The limits of the confidence interval. Must be coded as number in range from 0 to 1. Default value is 0.95.
#'
#' @param Inconclusive Inconclusive is a value that ranges from 0 to 1 that will identify the test range where the performance of the test is not acceptable and thus considered inconclusive. It represents the researcher tolerance of how good the test should be. If it is set to 0.95 (which is the default value), test results that have less than 0.95 sensitivity and specificity will be in the inconclusive range.
#'
#' @param Prevalence, pop.prevalence Population condition prevalence. The value of pop.prevalence will be passed to \code{\link{thresholds}} and The value of Prevalence will be passed to \code{\link{SS}}.
#'
#' @param Cost Represents wiegths of wrong classifications as cost(FN)/cost(FP). See \code{\link{thresholds}}.
#'
#' @param t.min,t.max,precision Test minimum, maximum and intervals to simulate the parametric estimation as \code{\link[base]{seq}}(t.min, t.max, precision). See \code{\link{SS}}.
#'
#' @param reverse It accepts 'auto', TRUE or FALSE. ROC analysis assumes that higher values are from subjects with the condition and lower values are from subjects without the condition. If it occurs the other way around, the ROC analysis and its interpretation must be reversed. See \code{\link{np.auROCc}} and \code{\link{SS}}.
#'
#' @param n.neurons Numeric vector containing the number of neurons of each layer. See \code{\link[AMORE]{newff}}.
#'
#' @param learning.rate.global Learning rate at which every neuron is trained. See \code{\link[AMORE]{newff}}.
#'
#' @param momentum.global Momentum for every neuron. See \code{\link[AMORE]{newff}}.
#'
#' @param error.criterium Criteria used to measure to proximity of the neural network prediction to its target. See \code{\link[AMORE]{newff}}.
#'
#' @param Stao Stao parameter for the TAO error criteria. See \code{\link[AMORE]{newff}}.
#'
#' @param hidden.layer Activation function of the hidden layer neurons. See \code{\link[AMORE]{newff}}.
#'
#' @param output.layer Activation function of the hidden layer neurons. See \code{\link[AMORE]{newff}}.
#'
#' @param method Preferred training method. See \code{\link[AMORE]{newff}}.
#'
#' @param report Logical value indicating whether the training function should keep quiet. See \code{\link[AMORE]{train}}.
#'
#' @param show.step Number of epochs to train non-stop until the training function is allow to report. See \code{\link[AMORE]{train}}.
#'
#' @param n.shows Number of times to report (if report is TRUE). See \code{\link[AMORE]{train}}.
#'
#' @param x The output of \code{TGROC} function for printing and ploting.
#'
#' @param ylab,xlab,ylim,xlim,roc.xlab,roc.ylab Values to be passed to \code{\link[graphics]{plot.deafult}} for labels and limits of the axis. If Plot.type = "TGROC", then former will be used for axis labels, if Plot.type = "ROC", then the later will be use for axis labels. xlim also accepts the value "auto". If this is the case, the test range from the SS analysis of the first option of the Plot argument will be used as horizontal limits of the plot. In essence, it will be defined by t.min and t.max.
#'
#' @param Plot.type Valid options are "TGROC" and "ROC". This will define the type of plot to return.
#'
#' @param Plot Valid options are "Binormal", "Non-parametric", and "NN-parametric" (for neural network). See \code{\link{SS}}. In the TGROC plot, this will define which sensitivity and specificity lines will be plotted. It may be one, two or all of them (although it become a poluted graph). The order they are set matters, as some options, e.g. Inconclusive argument, will pick the first option to show in the plot.
#'
#' @param Plot.inc.area Acceptable options are TRUE (default) or FALSE. If TRUE and Plot.type = "TGROC", it will display a gray shade in the plot (using the \code{\link[graphics]{polygon}}) representing the inconclusive area defined by the first method in the Plot argument (e.g. "Binormal").
#'
#' @param shade.args If Plot.type = "TGROC" and Plot.inc.area = TRUE, then this list of arguments will be passed to \code{\link[graphics]{polygon}} and plot a shade representing the inconclusive area. Internally, the coordinates will be extracted from in inconclusive limits from the first method defined ate the Plot argument.
#'
#' @param Plot.CL Acceptable options are TRUE or FALSE (default). If TRUE and Plot.type = "TGROC", it will pick the first method in the Plot argument (e.g. "Binormal"), and display lines in the plot representing a confidence band for this method.
#'
#' @param Plot.threshold Method of the decision threshold to be represented as a vertical line in the plot. Acceptable values are "None" (default), "Max Youden", "Max Accuracy", "Max Accuracy area", "Max DOR", "Min Error rate", "Min ROC distance", "Se = Sp", "Min MCT", "Max Efficiency". See \code{\link{thresholds}}.
#'
#' @param threshold.arg If Plot.type = "TGROC" and Plot.threshold is any value different from "None", this list of aguments with graphical parameters will be passed to \code{\link[graphics]{abline}} to plot a vertical line represented the chosen trheshold. Internally, the \code{\link[graphics]{abline}} v argument will be replaced by the test value from Plot.threshold argument.
#'
#' @param auto.title Logical. If TRUE, a subtitle will be automatically be placed at the bottom with the values of the inconclusive limits of the chosen threshold.
#'
#' @param title.arg A list of arguments to be passed to \code{\link[graphics]{title}}. TGROC automatically generates a subtitle with values ploted. Internally the argument sub is replaced by values from other options. Sometimes tiny adjustements are required.
#'
#'@param np.Se.args If Plot.type = "TGROC" and Plot argument is "Non-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the sensitivity line. Internally, the x and y arguments of this list will be replaced by the non-parametric sensitivity values.
#'
#' @param np.Se.ci.args If Plot.type = "TGROC" and Plot first argument is "Non-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the non-parametric sensitivity confidence band. Internally, the x and y arguments of this list will be replaced by the non-parametric sensitivity values.
#'
#' @param np.Sp.args If Plot.type = "TGROC" and Plot argument is "Non-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the spcificity line. Internally, the x and y arguments of this list will be replaced by the non-parametric specificity values.
#'
#' @param np.Sp.ci.args If Plot.type = "TGROC" and Plot first argument is "Non-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the non-parametric specificty confidence band. Internally, the x and y arguments of this list will be replaced by the non-parametric specificity values.
#'
#'
#' @param NN.Se.args If Plot.type = "TGROC" and Plot argument is "NN-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the sensitivity line. Internally, the x and y arguments of this list will be replaced by the NN sensitivity values.
#'
#' @param NN.Se.ci.args If Plot.type = "TGROC" and Plot first argument is "NN-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the NN-parametric sensitivity confidence band. Internally, the x and y arguments of this list will be replaced by the NN-parametric sensitivity values.
#'
#' @param NN.Sp.args If Plot.type = "TGROC" and Plot argument is "NN-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the specificity line. Internally, the x and y arguments of this list will be replaced by the NN specificity values.
#'
#' @param NN.Sp.ci.args If Plot.type = "TGROC" and Plot first argument is "NN-parametric" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the NN-parametric specificty confidence band. Internally, the x and y arguments of this list will be replaced by the NN-parametric specificity values.
#'
#' @param BN.Se.args If Plot.type = "TGROC" and Plot argument is "Binormal", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the sensitivity line. Internally, the x and y arguments of this list will be replaced by the binormal sensitivity values.
#'
#' @param BN.Se.ci.args If Plot.type = "TGROC" and Plot first argument is "Binormal" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the binormal sensitivity confidence band. Internally, the x and y arguments of this list will be replaced by the binormal sensitivity values.
#'
#' @param BN.Sp.args If Plot.type = "TGROC" and Plot argument is "Binormal", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the specificity line. Internally, the x and y arguments of this list will be replaced by the binormal specificity values.
#'
#' @param BN.Sp.ci.args If Plot.type = "TGROC" and Plot first argument is "Binormal" and Plot.CL = TRUE, then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the binormal specificty confidence band. Internally, the x and y arguments of this list will be replaced by the binormal specificity values.
#'
#' @param roc.np.line.args If Plot.type = "ROC" and Plot argument is "Non-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the empirical ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical sensitivity and specificity.
#'
#' @param roc.np.point.args If Plot.type = "ROC" and Plot argument is "Non-parametric" and Plot.threshold is different from "None", then this list of arguments will be passed to \code{\link[graphics]{points}} to plot a point representing the desired threshold on the empirical ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical point sensitivity and point specificity of the desired threshold.
#'
#' @param roc.NN.line.args If Plot.type = "ROC" and Plot argument is "NN-parametric", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the NN ROC curve. Internally, the x and y arguments of this list will be replaced by the NN sensitivity and specificity.
#'
#' @param roc.NN.point.args If Plot.type = "ROC" and Plot argument is "NN-parametric" and Plot.threshold is different from "None", then this list of arguments will be passed to \code{\link[graphics]{points}} to plot a point representing the desired threshold on the NN ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical point sensitivity and point specificity of the desired threshold.
#'
#' @param roc.BN.line.args If Plot.type = "ROC" and Plot argument is "Binormal", then this list of arguments will be passed to \code{\link[graphics]{lines}} to plot the binormal ROC curve. Internally, the x and y arguments of this list will be replaced by the binormal sensitivity and specificity.
#'
#' @param roc.BN.point.args If Plot.type = "ROC" and Plot argument is "Binormal" and Plot.threshold is different from "None", then this list of arguments will be passed to \code{\link[graphics]{points}} to plot a point representing the desired threshold on the binormal ROC curve. Internally, the x and y arguments of this list will be replaced by the empirical point sensitivity and point specificity of the desired threshold.
#'
#' @param grid If Plot.type = "ROC" and grid = TRUE, then a grid will be ploted. Set to FALSE to avoid ploting the grid.
#'
#' @param auto.legend Logical. If TRUE, then a legend will be automatically place in the plot. If Plot.type = "ROC", then the auto.legend will work only if Plot.threshold is different from "None".
#'
#' @param legend.args If auto.legend = TRUE, then this list of arguments will be passed to \code{\link[graphics]{legend}}. Internnaly, several arguments from lines and shade lists of arguments will be passed to \code{\link[graphics]{legend}}. If Plot.type = "TGROC", then auto.legend will use the BN.Se.args and BN.Sp.args to extract the line type and line colors. If Plot.type = "ROC", the auto.legend is way less flexible, e.g. the positions are fixed for all methods.
#'
#' @param ... Additional arguments passed either to \code{\link[base]{print}} or \code{\link[graphics]{plot.default}}
#'
#' @details
#' There are two main advantages of TG-ROC over ROC analysis: (1) for the uninitiated is much easier to understand how sensitivity and specificity changes with different thresholds; (2) and because of the graphical display is much easier to understand and estimate reasonable inconclusive test ranges. Occasionally, thresholds may be set outside the inconclusive range. This may happens with extreme values of Cost and population prevalence, or very skewed or unsusual test values distributions. If this is the case, perhaps the inconclusive range may not be of interest or not applicable. Also, if the test is too accurate or inconclusive range tolerance is set too low, then there may be no inconclusive range at all, because sensitivity and specificity may not be below this tolerance at the same time. If this is the case, setting a higher inconclusive tolerance may work. Tests results matching the threshold values will be considered a positive test.TGROC assumes that subjects with higher values of the test are with the target condition and those with lower values are without the target condition. Tests that behave like glucose (middle values are supposed to be normal and extreme values are supposed to be abnormal) will not be correctly analyzed. If lower test values are from subjects with the condition and higher values are from subject without the condition, the analysis and its interpretation must reversed. The validity measures such as Sensitivity, Specificity and Likelihood ratios and its confidence limits are estimated as  in \code{\link{diagnosis}} function. Non-parametric confidence bands are estimated with \code{\link{binom.CI}} and parametric with normal confidence interval. The parametric curve and validity measures are estimated with a neural network strategy using the \pkg{AMORE} package. Usually neural networks, uses a subset of the data to estimate weights, a subset to calibrate/validate the weights and a third subset to simulate the function. TGROC uses only the estimate and simulate steps, therefore there is no early stopping rule for the neural network parametric estimation to prevent overfitting. The only way to check the fit of the neural network is to visually compare with the non-parametric curve. If the curve looks weird or not good enough, than progressive slight changes in the momentum, learning rate, number of layers or other parameters should work fine. The AUC (area under the ROC curve) is estimated by the trapezoidal method (also known as Mann-Whitney statistic), its confidence interval is estimated by DeLong method. See \code{\link{np.auROCc}}. The AUC confidence limits should be used only to compare the AUC with the null value for AUC which is 0.5 and not to compare the AUC from different tests.
#'
#' @return A list of class TGROC with the following:
#' \itemize{
#'   \item \code{sample.size} The sample size.
#'   \item \code{sample.prevalence} The sample prevalence.
#'   \item \code{pop.prevalence} The population prevalence.
#'   \item \code{cost} The informed cost.
#'   \item \code{test.summary} A data.frame statistics such as mean and sd from the test.
#'   \item \code{inc} The desired inconclusive limit.
#'   \item \code{conf.limit} The desired confidece limit.
#'   \item \code{AUC} The estimated area under the ROC curve, SE and confidence limits.
#'   \item \code{SS} A data.frame with the empirical sensitivity and specificity trade-off (ROC) analysis.
#'   \item \code{np.inconclusive} A data.frame with the non-parametric inconclusive limits and their validity measures.
#'   \item \code{np.best.threshold} A data.frame with the non-parametric thresholds and their validity measures.
#'   \item \code{NN.SS} A data.frame with the neural network  sensitivity and specificity trade-off (ROC) analysis.
#'   \item \code{NN.inconclusive} A data.frame with the inconclusive limits from neural network smoothed ROC analysis and their validity measures.
#'   \item \code{NN.best.threshold} A data.frame with the NN thresholds and their validity measures.
#'   \item \code{BN.SS} A data.frame with the binormal sensitivity and specificity trade-off (ROC) analysis.
#'   \item \code{BN.inconclusive} A data.frame with the inconclusive limits from binormal ROC analysis and their validity measures.
#'   \item \code{BN.best.threshold} A data.frame with the binormal thresholds and their validity measures.
#' }
#'
#' @references
#'Greiner, M. (1996) Two-graph receiver operating characteristic (TG-ROC): update version supports optimization of cut-off values that minimize overall misclassification costs. J.Immunol.Methods 191:93-94.
#'
#'M. Greiner (1995) Two-graph receiver operating characteristic (TG-ROC): a Microsoft-EXCEL template for the selection of cut-off values in diagnostic tests. Journal of Immunological Methods. 185(1):145-146.
#'
#'M. Greiner, D. Sohr, P. Gobel (1995) A modified ROC analysis for the selection of cut-off values and the definition of intermediate results of serodiagnostic tests. Journal of immunological methods. 185(1):123-132.
#'
#' @seealso \code{\link{diagnosis}}, \code{\link{binom.CI}}, \code{\link{SS}}, \code{\link{thresholds}}, \code{\link{np.auROCc}}
#'
#' @examples
#'
#'data(rocdata)
#'x <- TGROC(ref = rocdata$Gold, test = rocdata$test2, precision = .005)
#'
#' # Printing just the NN part of the analysis
#' print(x, print.type = "NN-parametric")
#'
#' # Ploting just the empirical ROC curve with Max Youden threshold.
#' plot(x, Plot.type = "ROC", Plot = "Non-parametric", Plot.threshold = "Max Youden")
#'
#' # Ploting just the NN ROC curve
#' # The plot does not reach the corner because there is no 100% Sensitivity at any threshold.
#' plot(x, Plot.type = "ROC", Plot = "NN-parametric", Plot.threshold = "Max Youden")
#'
#' # Ploting just the Binormal ROC curve
#' plot(x, Plot.type = "ROC", Plot = "Binormal", Plot.threshold = "Max Youden")
#'
#' # Ploting overploted curves. The Binormal curve is way different from the remaining.
#' plot(x, Plot.type = "ROC", Plot = c("Non-parametric","NN-parametric","Binormal"), Plot.threshold = "Max Youden")
#'
#' # Ploting the TGROC curve overploting the Non-parametric and the NN- parametric.
#' plot(x, Plot = c("Non-parametric","NN-parametric"))
#'
#' # TGROC curve overploting the Non-parametric (and confidence band) and the NN- parametric.
#' # If legend is not good enough, scoot up a little bit decreasing the inset.
#' plot(x, Plot = c("NN-parametric", "Non-parametric"), Plot.CL = TRUE,
#'      legend.args = list(x = "top", border = NA, bty = "n", xpd = NA,
#'                              inset = -.25, ncol = 2, cex = .8))
#'
#' # Changing the order of the Plot argument changes the inconclusive area shade.
#' plot(x, Plot = c("NN-parametric","Non-parametric"))
#'
#' # The Binormal fit looks worst than the NN fit
#' plot(x, Plot = c("Binormal","Non-parametric"))
#'
#' data(tutorial)
#' # Running TGROC with a very accurate test.
#' # warnings concerning ties and concerning the inconclusive limits.
#' x <- TGROC(ref = ifelse(tutorial$Gold == "pos", 1, 0), test = tutorial$Test_A, precision = .005)
#'
#' # All lower inconclusive limits are higher than the upper inconclusive limit
#' plot(x, Plot = c("NN-parametric","Non-parametric"))
#' plot(x, Plot = c("Binormal","Non-parametric"))
#' plot(x, Plot = c("Non-parametric"))
#'
#' # Increasing the inconclusive solves the issue.
#' x <- TGROC(ref = ifelse(tutorial$Gold == "pos", 1, 0), test = tutorial$Test_A, precision = .005, Inconclusive = .99)
#'
#' # All lower inconclusive limits are higher than the upper inconclusive limit
#' # Nevertheless, check how the inconclisve area chages as method Changes.
#' plot(x, Plot = c("NN-parametric","Non-parametric"))
#' plot(x, Plot = c("Binormal","Non-parametric"))
#' plot(x, Plot = c("Non-parametric"))
#'
#' # Ploting with different thrensholds without the inconclusive area.
#' x <- TGROC(ref = ifelse(tutorial$Gold == "pos", 1, 0), test = tutorial$Test_B, precision = .005, Inconclusive = .90)
#' plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.inc.area = FALSE, Plot.threshold = "Se = Sp")
#' plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.inc.area = FALSE, Plot.threshold = "Min ROC distance")
#' plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.inc.area = FALSE, Plot.threshold = "Max DOR")
#'
#' # Ploting the inconclusive area and the threshold simultaneously.
#' # The subtitle may get overploted, expanding the margins and adjusting the line may help.
#' par(mar = c(6.1, 4.1, 4.1, 2.1))
#' plot(x, Plot = c("NN-parametric","Non-parametric"), Plot.threshold = "Max Efficiency",
#'      title.arg = list(cex.sub = 0.7, line = 5))
#'
#' # Back to original margins.
#' par(mar = c(5.1, 4.1, 4.1, 2.1))
#'
#' rm(rocdata, tutorial, x)
#'
#' @export
TGROC <- function(ref,
                test,
                Cost = 1,
                CL = 0.95,
                Inconclusive = 0.95,
                Prevalence = NULL,
                pop.prevalence = NULL,
                t.max = NULL,
                t.min = NULL,
                precision = 0.05,
                reverse = "auto",
                n.neurons = c(1, 5, 1),
                learning.rate.global = 1e-2,
                momentum.global = 0.3,
                error.criterium = "LMS",
                Stao = NA,
                hidden.layer = "sigmoid",
                output.layer = "sigmoid",
                method = "ADAPTgdwm",
                report = FALSE,
                show.step = 5000,
                n.shows = 1){

  # Preventing wrong inputs and warnings ---------------------------------------
  if (any(levels(as.factor(ref)) != c(0, 1))) {
    stop("Your reference standard must be coded as 0 (absence) and 1 (presence). Check reference categories!")
  }
  if (is.null(precision) || !is.numeric(precision)) {
      stop("Precision must be set to a numeric value!")
  }

  # Making the non-parametric trade-off for Se and Sp (ROC analysis) -----------
  SeSp <- SS(ref, test, reverse = reverse, CL = CL, pop.prevalence = Prevalence)
  AUC <- np.auROCc(ref, test, reverse = reverse, CL = CL)

  # Setting requeired additional objects
  sample.prevalence <- SeSp$sample.prevalence
  sample.size <- SeSp$sample.size
  names(sample.prevalence) <- c("Condition prevalence in the sample")
  names(sample.size) <- c("Sample size")
  if (is.null(Prevalence)) { Prevalence <- sample.prevalence }
  names(Prevalence) <- c("Informed disease prevalence in the population.")
  cost <- Cost
  names(cost) <- c("Informed costs(FN)/costs(FP).")
  conf.limit <- CL
  inc <- Inconclusive
  names(inc) <- "Inconclusive tolerance level."

  # Overall test summary -------------------------------------------------------
  test.summary <- rbind(c(summary(test), sd(test)),
                        c(summary(test[which(ref == 1)]), sd(test[which(ref == 1)])),
                        c(summary(test[which(ref == 0)]), sd(test[which(ref == 0)])))
  colnames(test.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "SD")
  rownames(test.summary) <- c("Overall", "With the condition", "Without the condition")

  #Setting the decision Thresholds ---------------------------------------------
  np.best.threshold <- thresholds(SeSp, pop.prevalence = pop.prevalence, Cost = Cost)

  # Extracting the inconclusive ranges
  np.inconclusive <- inc.limits(SeSp, Inconclusive = Inconclusive)

  # Fitting the Neural Network for smoothed threshold estimation
  NN.SeSp <- NN.SS(SeSp, t.max = t.max, t.min = t.min, precision = precision,  n.neurons = n.neurons,  learning.rate.global = learning.rate.global,  momentum.global = momentum.global,  error.criterium = error.criterium,  Stao = Stao, hidden.layer = hidden.layer,  output.layer = output.layer, method = method,  report = report,  show.step = show.step,  n.shows = n.shows,  CL = CL)

  # Extracting the NN decision thresholds --------------------------------------------------
  NN.best.threshold <- thresholds(NN.SeSp, pop.prevalence = pop.prevalence, Cost = Cost)

  # Extracting the inconclusive ranges -----------------------------------------------------
  NN.inconclusive <- inc.limits(NN.SeSp, Inconclusive = Inconclusive)

  # Fitting the Binormal SS object
  BN.SeSp <- BN.SS(ref = ref, test = test, CL = CL, t.max = t.max, t.min = t.min, precision = precision, pop.prevalence = Prevalence)

  #Setting the decision Thresholds --------------------------------------------------
  BN.best.threshold <- thresholds(BN.SeSp, pop.prevalence = pop.prevalence, Cost = Cost)

  # Extracting the inconclusive ranges
  BN.inconclusive <- inc.limits(BN.SeSp, Inconclusive = Inconclusive)

  # Warnings Regarding the thresholds out of the inconclusive range
  if(np.inconclusive[1,1] > np.inconclusive[2,1]){
     warning("Non-parametric lower inconclusive limit is higher than upper inconclusive limit. \n Either ROC analysis should be reversed or Inconclusive argument should have higher value.")
  }
  if(NN.inconclusive[1,1] > NN.inconclusive[2,1]){
     warning("NN-parametric lower inconclusive limit is higher than upper inconclusive limit. \n Either ROC analysis should be reversed or Inconclusive argument should have higher value.")
  }
  if(BN.inconclusive[1,1] > BN.inconclusive[2,1]){
    warning("Binormal lower inconclusive limit is higher than upper inconclusive limit. \n Inconclusive argument should have higher value.")
  }

  # if(any(np.best.threshold[,1] > np.inconclusive[2,1])){
  #    warning("At least one of the non-parametric threshold is higher then upper inconclusive limit.")
  # }
  # if(any(np.best.threshold[,1] < np.inconclusive[1,1])){
  #    warning("At least one of the non-parametric threshold is lower then lower inconclusive limit.")
  # }
  # if(any(NN.best.threshold[,1] > NN.inconclusive[2,1])){
  #    warning("At least one of the NN-parametric best threshold is higher then upper inconclusive limit.")
  # }
  # if(any(NN.best.threshold[,1] < NN.inconclusive[1,1])){
  #    warning("At least one of the NN-parametric best threshold is lower then lower inconclusive limit.")
  # }
  # if(any(BN.best.threshold[,1] > BN.inconclusive[2,1])){
  #   warning("At least one of the Binormal best threshold is higher then upper inconclusive limit.")
  # }
  # if(any(BN.best.threshold[,1] < BN.inconclusive[1,1])){
  #   warning("At least one of the Binormal best threshold is lower then lower inconclusive limit.")
  # }
  output <- list(sample.size = sample.size,
                sample.prevalence = sample.prevalence,
                pop.prevalence = Prevalence,
                cost = cost,
                test.summary = test.summary,
                inc = inc,
                conf.limit = conf.limit,
                AUC = AUC,
                SS = SeSp$table,
                np.inconclusive = np.inconclusive,
                np.best.threshold = np.best.threshold,
                NN.SS = NN.SeSp$table,
                NN.inconclusive = NN.inconclusive,
                NN.best.threshold = NN.best.threshold,
                BN.SS = BN.SeSp$table,
                BN.inconclusive = BN.inconclusive,
                BN.best.threshold = BN.best.threshold
                )
  class(output) <- "TGROC"
  output
}
