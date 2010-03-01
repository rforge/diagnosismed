\name{ROC}
\alias{ROC}
\alias{print.ROC}
\alias{plot.ROC}
\title{Draw a ROC curve, estimate good cut-offs and compute validity measures for each cut-off}
\description{
  Draw a non-parametric (empirical) ROC - receiver operating characteristic - curve and compute test sensitivity, specificity, predictive values and likelihood ratios (and respective confidence limits) for each decision threshold. Estimate good decision threshold by a variety of methods.}
\usage{
ROC(gold,
    test,
    CL = 0.95,
    Cost = 1,
    Prevalence = 0,
    Plot = TRUE,
    Plot.point = "Min.ROC.Dist",
    p.cex=1,
    Full = FALSE,
    Print = TRUE)
\method{plot}{ROC}(x,
    Plot.point = "Min.ROC.Dist",
    cex.sub=.85,
    p.cex=1,
    ...)
\method{print}{ROC}(x,Full=FALSE,...)    
}
\arguments{
  \item{gold}{The reference standard. A column in a data frame indicating the classification by the reference test. The reference standard must have two levels: must be coded either as 0 - without target condition - or 1 - with the target condition; or could be coded \link[base]{as.factor} with the words "negative" - without target condition - and "positive"  - with the target condition.}
  \item{test}{The index test or test under evaluation. A column in a dataframe or vector indicating the test results in a continuous scale. Might also work with discrete or ordinal scale.}
  \item{CL}{Confidence limit. The limits of the confidence interval. Must be coded as number in a range from 0 to 1. Default value is 0.95}
  \item{Cost}{Cost = cost(FN)/cost(FP). MCT will be used to estimate a good cut-off. It is a value in a range from 0 to infinite. Could be financial cost or a health outcome with the perception that FN are more undesirable than FP (or the other way around). This item will run into MCT (misclassification cost term) - (1-prevalence)*(1-Sp)+Cost*prevalence(1-Se). Cost=1 means FN and FP have even cost. Cost = 0.9 means FP are 10 percent more costly. Cost = 0.769 means that FP are 30 percent more costly. Cost = 0.555 means that FP are 80 percent more costly. Cost = 0.3 means that FP are 3 times more costly. Cost = 0.2 means that FP are 5 times more costly. Also, it can be more easily inserted as any ratio such as 1/2.5 or 1/4.}
  \item{Prevalence}{Prevalence of the disease in the population who the test will be performed. It must be a value from 0 to 1. If left 0 (the default value), this will be replaced by the disease prevalence in the sample. This value will be used in the MCT and Efficiency formulas to estimate good  cut-offs.}
  \item{Plot}{If FALSE, the ROC curve plot will not be displayed. Default is TRUE. If default options of graphics parameters from \link[graphics]{par} are not satisfactory, then the sugestion is to assign the output into a object and call options from \command{plot} }
  \item{Plot.point}{The method of best cut-off estimation which will be displayed at ROC curve as a dot. Default is "Min.ROC.Dist". Possible options are: 
        
        "None" - only the AUC in the legend will appear;
        
        "Max.Accuracy" - the cut-off which maximize the accuracy;
        
        "Max.DOR" - the cut-off which maximize the diagnostic odds ratio;

        "Error.rate" - the cut-off which minimizes the error rate;

        "Max.Accuracy.area" - the cut-off which maximize the accuracy area;

        "Max.Sens+Spec" - the cut-off which maximize the sum of sensitivity with specificity;

        "Max.Youden" - the cut-off which maximize the Youden index;

        "Se=Sp" - the cut-off which Sensitivity is equal to Specificity;

        "Min.ROC.Dist" - the cut-off which minimize the distance between the curve and the upper left corner of the graph;

        "Max.Efficiency" - the cut-off which maximize the efficiency;

        "Min.MCT" - the cut-off which minimize the misclassification cost term.
}
  \item{cex.sub}{The magnification to be used for sub-titles relative to the current setting of cex. See \link[graphics]{par}}
  \item{p.cex}{Symbol expansion - a numerical vector - passed to points. See \link[graphics]{points}}
  \item{Print}{If FALSE, no results (detailed below in values section) will be displayed on the output window. Default is TRUE}
  \item{Full}{If TRUE, a table with sensitivity, specificity, predictive values and likelihood ratios (and respective confidence limits) for each decision threshold will be displayed. Default is FALSE.}
  \item{x}{For the \command{plot} and \command{print} functions, x is an object storing \command{ROC} function output.}    
  \item{...}{Other plot or print parameters form \link[graphics]{plot.default}}  
}
\details{
   Tests results matching the cut-off values will be considered a positive test. \command{ROC} assumes that subjects with higher values of the test are with the target condition and those with lower values are without the target condition. Tests that behave like glucose (middle values are supposed to be normal and extreme values are supposed to be abnormal) and immunefluorescence (lower values - higher dilutions - are suppose to be abnormal) will not be correctly analyzed. In the latter, multiplying the test results by -1 or other transformation before analysis could make it work. The AUC (area under the ROC curve) is estimated by the trapezoidal method (also known as Mann-Whitney statistic), its confidence interval is estimated by DeLong method. The AUC confidence limits should be used only to compare the AUC with the null value for AUC which is 0.5 and not to compare the AUC from different tests. The validity measures such as Sensitivity, Specificity and Likelihood ratios and its confidence limits are estimated as in \code{\link{diagnosis}} function. If \command{ROC} output is assign to an object (see example), tests results could be easily exported to a spreadsheet and other graphics that might be of interest could be easily done.
   
   Diagnostic odds ratio: \eqn{DOR = (TP*TN)/(FN*FP); the same as: DOR = PLR/NLR}

   Accuracy area: \eqn{AA = (TP*TN)/((TP+FN)*(FP+TN))}

   Youden index: \eqn{Y = Se+Sp-1; the same as: Y = Se-FPR}

   Minimum ROC distance: \eqn{m ROC Dis = (Sp-1)^2+(1-Se)^2}

   Efficiency: \eqn{Ef = Se*prevalence+(1-prevalence)*Sp}

   Misclassification Cost Term: \eqn{MCT = (1-prevalence)*(1-Sp)+(cost(FN)/cost(FP))*prevalence*(1-Se)}
}
\value{
  \item{pop.prevalence}{The disease prevalence informed by the user. If not informed, it will be the same as the sample prevalence.}
  \item{sample.prevalence}{The disease prevalence in the sample}
  \item{sample.size}{The number of subjects analyzed}
  \item{test.summary}{A table showing the quintiles, mean and standard deviation of overall test results, test results from those with the target condition and without the target condition}
  \item{AUC.summary}{A table showing the AUC estimated by DeLong method (trapezoidal)
        and its confidence limits.}
  \item{test.best.cutoff}{A table showing the best cut-offs estimated by methods described above, its corresponding sensitivity, specificity and positive likelihood ratio (and their confidence limits)}
}
\references{
Knotterus. The Evidence Based Clinical Diagnosis; BMJBooks, 2002.

Xiou-Hua Zhou, Nancy A Obuchowsky, Donna McClish. Statistical Methods  in diagnostic Medicine; Wiley, 2002.

Simel D, Samsa G, Matchar D (1991). Likelihood ratios with confidence: Sample size estimation for diagnostic test studies. Journal of Clinical Epidemiology 44: 763 - 770

S.B. Cantor, C.C. Sun, G. Tortolero-Luna, R. Richards-Kortum, and M. Follen. (1999) A comparison of C/B ratios from studies using receiver operating characteristic curve analysis. Journal of Clinical Epidemiology, 52(9):885-892.

Greiner, M. (1996) Two-graph receiver operating characteristic (TG-ROC): update version supports optimisation of cut-off values that minimize overall misclassification costs. J.Immunol.Methods 191:93-94.

Gengsheng Qin, Lejla Hotilovac. Comparison of non-parametric confidence intervals for the area under the ROC curve of a continuous-scale disagnostic test. Statistical Methods in Medical Research 2008; 17:207-221.

}

\author{Pedro Brasil; Beranrdo Rangel Tura - \email{diagnosismed-list@lists.r-forge.r-project.org}}
\note{Bug reports, malfunctioning, or suggestions for further improvements or contributions can be sent, preferentially, through the DiagnosisMed email list, or R-Forge website \url{https://r-forge.r-project.org/projects/diagnosismed/}.
}
\seealso{\link[epitools]{binom.conf.int},\code{\link{diagnosis}},\code{\link{interact.ROC}},\code{\link{TGROC}},\link[ROCR]{performance},\link[nonbinROC]{contROC}}
\examples{
# loading a dataset
data(rocdata)
# Attaching the data set.
attach(rocdata)
# A little description of the data set to check if it is ok!
str(rocdata)
# Running ROC analysis with the full table option 
# and storing ROC objects into 'x' from which there are tables to draw the graphs below.
x<-DiagnosisMed::ROC(Gold,test2,Full=TRUE) 
# There is no need to stick the package name before the function if it is loaded as first in search path!
# Adding a title to the graph.
title(main="ROC graph")
# Some graphs that may be of interest. Validity measures at each test value. 
# Setting the plot window to get nine graphs
# Some graphs showing some validity measures and some indexes variations used to choose good cut-offs
par(mfrow=c(3,3))
plot(x$test.diag.table$test.values,x$test.diag.table$DOR,type="l",ylab="DOR",xlab="Test values")
title(main="Test values x DOR")
plot(x$test.diag.table$test.values,x$test.diag.table$MCT,type="l",ylab="MCT",xlab="Test values")
title(main="Test values x MCT")
plot(x$test.diag.table$test.values,x$test.diag.table$Efficiency,type="l",ylab="Efficiency",xlab="Test values")
title(main="Test values x Efficiency")
plot(x$test.diag.table$test.values,x$test.diag.table$Youden,type="l",ylab="Youden index",xlab="Test values")
title(main="Test values x Youden index")
plot(x$test.diag.table$test.values,x$test.diag.table$PLR,type="l",ylim=c(0,49),ylab="Likelihood ratios",xlab="Test values")
lines(x$test.diag.table$test.values,x$test.diag.table$NLR,type="l",lty=2)
legend("right",lty=c(1,2),legend=c("PLR","NLR"),bty = 'n')
title(main="Test values x Likelihood ratios")
plot(x$test.diag.table$test.values,x$test.diag.table$PPV,type="l",ylab="Predictive values",xlab="Test values")
lines(x$test.diag.table$test.values,x$test.diag.table$NPV,type="l",lty=2)
legend("bottomright",lty=c(1,2),legend=c("PPV","NPV"),bty = 'n')
title(main="Test values x Predictive values")
plot(x$test.diag.table$test.values,x$test.diag.table$Accuracy.area,type="l",,ylab="Accuracy area",xlab="Test values")
title(main="Test values x Accuracy area")
plot(x$test.diag.table$test.values,x$test.diag.table$MinRocDist,type="l",,ylab="ROC distance",xlab="Test values")
title(main="Test values x ROC distance")
plot(x$test.diag.table$test.values,x$test.diag.table$Accuracy,type="l",ylab="Error rate & Accuracy",xlab="Test values",ylim=c(0,1))
lines(x$test.diag.table$test.values,x$test.diag.table$Error.rate,type="l",lty=2)#,xlim=c(0,2.5))
legend("bottomright",lty=c(1,2),legend=c("Accuracy","Error rate"),bty = 'n')
par(mfrow=c(1,1))

# Also, results from ROC analysis could easily exported to a spreadsheet file or to a odt file by OdfWeave.
# Exporting the full table:
# write.csv(x$test.diag.table[,-c(2:5,24:34)],'MytestFulltable.csv')
# Exporting AUC summary:
# write.csv(x$AUC.summary,'MytestAUC.csv')      
# Exporting Test summary:
# write.csv(x$test.summary,'MytestSummary.csv')
# Exporting Test best-cut-offs table:
# write.csv(x$test.best.cutoff,'MytestBestcutof.csv')                                                
rm(rocdata,x)
}
\keyword{iplot}
\keyword{univar}
\keyword{htest}
