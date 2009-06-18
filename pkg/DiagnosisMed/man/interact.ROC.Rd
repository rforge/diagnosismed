\name{interact.ROC}
\alias{interact.ROC}
\title{Interactively draw a ROC curve with your data}
\description{
Draw a ROC curve with the user data, interactively, sliding a button, and watch how changes at the cutoff, correlate with the changes of Sensitivity and Specificity while building the ROC curve itself.
}
\usage{
interact.ROC(gold, test)
}
\arguments{
  \item{gold}{A column in a data frame or a vector indicating the classification by the reference test. Must be coded either as 0 - without target condition - or 1 - with the target condition.}
  \item{test}{A column in a data frame or a vector indicating the test under study (index test) results. Must be \link[base]{numeric}.}
}
\details{
interact.ROC is a call from roc.demo function in TeachingDemo package. The difference is that interact.ROC allow the input data be displayed as usually is in diagnostic studies, a column with the test result and the other with the reference standard results. Inside this function the "ask" option - which controls the "next" button to see the next graph - is turned off - par(ask=FALSE). Also, other options must be set to this function to work fine. Type "options()". The option ask and device.ask.default should be set as FALSE, as they usually are by default. To turn it on again later, type par(ask=TRUE). The test must have a rationale that higher values of the index test belong to those with the target disease and those with lower values belong to those without the target disease. If this is not the case, the suggestion is to transform the tests results by multiplying it by -1 before running interact.ROC.
}
\value{
interact.ROC generates two graphs in the same window: the upper graph is a ROC graph (Sensitivity on the vertical axis and 1-Specificity on the horizontal axis); the lower graph is a density plot (the density on the vertical axis and the test cut-off (or threshold) on the horizontal axis). With a sliding button is possible to interact and see how  the Sensitivity and Specificity changes while changing the cut-off. In the upper graph 
 the cut-off is represented by the different dots and the purple line represents the distance to the "optimal" threshold. At the lower graph, the red line and dashes represent the density and the test result from those without the target condition respectively. While the blue ones represent those with the target condition. If the dashes are at the bottom of the lower graph then they are classified as without the target condition, if at the top, with the target condition. The green vertical line represents the cut-off and changes with the sliding button. The cut-off itself can be seen right above the sliding button and the respective sensitivity and specificity at the bottom of the graph window. 

}
\references{

JA Knotterus, The Evidence Based Clinical Diagnosis; BMJBooks, 2002

}
\author{Pedro Brasil - \email{diagnosismed-list@lists.r-forge.r-project.org}}
\note{Bug reports, malfunctioning, or suggestions for further improvements or contributions can be sent, preferentially, through the DiagnosisMed email list, or R-Forge website \url{https://r-forge.r-project.org/projects/diagnosismed/}.
}
\seealso{\code{\link{diagnosis}},\code{\link{ROC}},\code{\link{TGROC}},\link[ROCR]{performance},\link[TeachingDemos]{roc.demo},\link[nonbinROC]{contROC}}
\examples{
data(rocdata)
attach(rocdata)
interact.ROC(Gold,test2)
rm(rocdata)
}
\keyword{dynamic}
\keyword{iplot}

