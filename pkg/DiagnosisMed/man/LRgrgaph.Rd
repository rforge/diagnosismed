\name{LRgraph}
\alias{LRgraph}
\title{Comparing diagnositic tests: a simple graphic using likelihood ratios.}
\description{
LRgraph graphicaly compares two or more (all of them with the first test) diagnostic tests with binary results through their likelihood ratios, based on the rationale that the predictive hability of a test is a more interresting charachteristic than sensitivity and/or specificity. It is possible to see thruogh the graph that if the tests with smaller sensitivity or specifity may have superior predictive habilty, that is, increases the prediction hability with small sensitivity/specificity trade-off.
}
\usage{
LRgraph(a = cbind(t1, t2), lwd = 2, lty = 1, cex = 1, leg.cex = 1.5, pt.cex = 2,...)
}
\arguments{
  \item{a}{a is a matrix composed by tow or more tests. Therefore, the user should edit only what is between the parethesis in a=cbind().The user may insert as many tests as one wishes. See t1, t2 below.}
  \item{t1,t2}{ t1, t2... tn are objects created by the \code{\link{diagnosis}} function (see example below). Therefore, one should first analyze data for each single test and store it in a object. Later, stick these objects in the LRgraph command. The objects names will appear in the graph legend.}
  \item{lwd}{Line width. See \link[graphics]{par},\link[graphics]{points},\link[graphics]{legend}}
  \item{lty}{Line type. See \link[graphics]{par}}
  \item{cex}{Symbols and text size. See \link[graphics]{par},\link[graphics]{points}}
  \item{leg.cex}{Legend text size, this will replace the cex option in the legend. See \link[graphics]{legend}}
  \item{pt.cex}{Size of the symbols in the legend. See \link[graphics]{legend}}
  \item{...}{Other graph parameters. See \link[graphics]{plot.default}}
}
\details{When a diagnostic test has both sensitivity and specificity higher than a competing test is easy to see that the former is superior than the later. However, sometimes a test may have superior sensitivity and inferior specificity (or the other way around). In this case, a good decision may be toward the test that have a better prediction ability. The graph visually helps the user to see and compare these abilities. The graph is very similar to the ROC graph. The vertical and horizontal axis have the same length as the ROC graph. However, the diagnostic tests are represented as dots instead of curves. The solid line passing through (0,0) is the likelihood ratio positive-line and the solid line passing through (1,1) is the likelihood ratio negative-line. Both negative and positive likelihood are numerically equivalent to the slopes of the solid lines. The solid lines split the graph into four areas (run the example). Also, there are dashed lines representing the sensitivity and specificity of the first test plotted. One may see that there are areas that a test may have superior sensitivity (or specificity) and yet the dot may be below the likelihood solid line. That is because the sensitivity / specificity trade-off is not reasonable, making the test with less predictive ability. 

}
\value{
Returns only a graph which is divided in four areas by the black solid lines (the likelihood ratios of the firts test). The interpretation of the comparisons will depend on which area the competing tests will fall in. See and run the example to have the idea on how interpretation must be done.
}
\references{Biggerstaff, B.J. Comparing diagnostic tests: a simple graphic using likelihood ratios. Statistics in Medicine. 2000; 19(5):649-663}
\author{Pedro Brasil - \email{diagnosismed-list@lists.r-forge.r-project.org}}
\note{Bug reports, malfunctioning, or suggestions for further improvements or contributions can be sent, preferentially, through the DiagnosisMed email list, or R-Forge website \url{https://r-forge.r-project.org/projects/diagnosismed/}.
}
\seealso{\code{\link{diagnosis}}}
\examples{
# Making tests with diagnosis function with different performances for comparison.
# mytest5 is the one which all others will be compared with.
mytest5<-diagnosisI(80,20,20,80,print=FALSE)

# mytest1 has higher sensitivity and specificity.
# mytest1 is overall superior compared to mytest5.
mytest1<-diagnosisI(90,10,10,90,print=FALSE)

# mytest2 has lower sensitivity but higher specificity.
# mytest2 is better to identify the absence of the target condition compared to mytest5.
mytest2<-diagnosisI(72,28,3,97,print=FALSE)

# mytest3 has higher sensitivity but lower specificity.
# mytest3 is better to identify the presence of the target condition compared to mytest5.
mytest3<-diagnosisI(92,8,37,63,print=FALSE)

# mytest41 has lower sensitivity and specificity.
# mytest41 is overall inferior compared to mytest5.
mytest41<-diagnosisI(72,28,35,65,print=FALSE)

# mytest42 has lower sensitivity but higher specificity.
# Yet, mytest42 is overall inferior compared to mytest5.
mytest42<-diagnosisI(82,18,42,58,print=FALSE)

# But that becomes clear only after ploting the tests.
LRgraph(a=cbind(mytest5,mytest1,mytest2,mytest3,mytest41,mytest42),cex=2.5)

# The texts below are not part of the function but helps to understand the areas
text(x=.5, y =.5, labels ="Area 4: Overall inferior", col="lightgray",cex=.8)
text(x=.5, y =1, labels ="Area 2: Absence", col="lightgray",cex=.8)
text(x=.07, y =.73, labels ="Area 3: Presence", col="lightgray",cex=.8)
text(x=.1, y =1, labels ="Area 1: Overall superior", col="lightgray",cex=.8)

rm(mytest1)
rm(mytest2)
rm(mytest3)
rm(mytest41)
rm(mytest42)
rm(mytest5)
}
\keyword{iplot}
