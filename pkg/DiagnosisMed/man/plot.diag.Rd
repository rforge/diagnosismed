\name{plot.diag}
\alias{plot.diag}
\title{Nomogram from 2x2 table diagnostic test}
\description{
This function draw a simple nomogram from an object created by diagnosis or
diagnosisI
}
\usage{
\method{plot}{diag}(x, print = FALSE, ...)
}
\arguments{
  \item{x}{An object from diagnosis or diagnosisI}
  \item{print}{If TRUE, plot.diag will print in the output window the table
  with all pre-test and its corresponding post-test probabilities. Default is
  FALSE}
  \item{\dots}{Further generic options for plot function}
}
\details{
plot.diag will draw a very simple nomogram as many examples from wikipedia
 \url{http://en.wikipedia.org/wiki/Nomogram}. This is not a generic nomogram
 as shown in many evidenced based medicine texts, because this one shows only
 pre-test and post-test variations with a fixed positive likelihood ratio. This
 likelihood is a statistic from an object created by diagnosis or diagnosisI function.
 Its usage is the same as applying the Bayes theorem where the pre-test odds times
 positive likelihood ratio equals the pos-test odd (transforming the odds to
 probabilities). To use it, draw, with a rule, a vertical line from a desired
 pre-test  probability, and to find the corresponding post-test probability, draw a
 horizontal line from the intersection of the curve and the vertical line toward the vertical axis.
}

\references{

JA Knotterus, The Evidence Based Clinical Diagnosis; BMJBooks, 2002

}
\author{Pedro Brasil - \email{diagnosismed-list@lists.r-forge.r-project.org}}
\note{Bug reports, malfunctioning, or suggestions for further improvements or contributions can be sent, preferentially, through the DiagnosisMed email list, or R-Forge website \url{https://r-forge.r-project.org/projects/diagnosismed/}.
}
\seealso{\code{\link{diagnosis}},\link[Design]{nomogram},\link[epiR]{epi.nomogram},\code{\link{ROC}},\code{\link{LRgraph}}}
\examples{
#This will draw a ROC curve
mytest<-diagnosisI(364,22,17,211, print=TRUE, plot=TRUE)
#This is different, and will draw a nomogram
mytest<-diagnosisI(364,22,17,211)
plot(mytest)
#Draw a nomogram and print it's table - printing is optional
mytest<-diagnosisI(364,22,17,211)
plot(mytest,print=TRUE)
# The same as
plot(diagnosisI(364,22,17,211),print=TRUE)

}
\keyword{univar}
\keyword{htest}
\keyword{aplot}
