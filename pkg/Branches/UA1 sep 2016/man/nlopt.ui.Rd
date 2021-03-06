% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlopt.ui.R
\name{nlopt.ui}
\alias{nlopt.ui}
\title{Function for the determination of the population thresholds an inconclusive interval for bi-normal distributed test scores.}
\usage{
nlopt.ui(Se = 0.55, Sp = 0.55, mu0 = 0, sd0 = 1, mu1 = 1, sd1 = 1,
  intersection = NULL, start = NULL, print.level = 0)
}
\arguments{
\item{Se}{(default = .55). Desired sensitivity of the test scores within the uncertain interval. A value <= .5 is not allowed.}

\item{Sp}{(default = .55). Desired specificity of the test scores within the uncertain interval. A value <= .5 is not allowed.}

\item{mu0}{Population value or estimate of the mean of the test scores of the persons without the targeted condition.}

\item{sd0}{Population value or estimate of the standard deviation of the test scores of the persons without the targeted condition.}

\item{mu1}{Population value or estimate of the mean of the test scores of the persons with the targeted condition.}

\item{sd1}{Population value or estimate of the standard deviation of the test scores of the persons with the targeted condition.}

\item{intersection}{Default NULL. If not null, the supplied value is used as the estimate of the intersection of the two bi-normal distributions. Otherwise, it is calculated.}

\item{start}{Default NULL. If not null, the first two values of the supplied vector are used as the starting values for the \code{nloptr} optimization function.}

\item{print.level}{Default is 0. The option print_level controls how much output is shown during the optimization process. Possible values: 0 (default)    no output; 1    show iteration number and value of objective function; 2    1 + show value of (in)equalities; 3    2 + show value of controls.}
}
\value{
List of values:
\describe{
  \item{$status: }{Integer value with the status of the optimization (0 is success).}
  \item{$message: }{More informative message with the status of the optimization}
  \item{$results: }{Vector with the following values:}
     \itemize{
      \item{exp.Sp.ui: }{The population value of the specificity in the Uncertain Interval, given mu0, sd0, mu1 and sd1. This value should be very near the supplied value of Sp.}
      \item{exp.Sp.ui: }{The population value of the sensitivity in the Uncertain Interval, given mu0, sd0, mu1 and sd1. This value should be very near the supplied value of Se.}
      \item{mu0: }{The value that has been supplied for mu0.}
      \item{sd0: }{The value that has been supplied for sd0.}
      \item{mu1: }{The value that has been supplied for mu1.}
      \item{sd1: }{The value that has been supplied for sd1.}
      }
  \item{$solution: }{Vector with the following values:}
     \itemize{
      \item{L: }{The population value of the lower threshold of the Uncertain Interval.}
      \item{U: }{The population value of the upper threshold of the Uncertain Interval.}
      }
}
}
\description{
Function for the determination of the population thresholds an inconclusive interval for bi-normal distributed test scores.
}
\details{
The function can be used to determinate the uncertain interval of two bi-normal distributions.
The Uncertain Interval is defined as an interval below and above the intersection of the two distributions, with a sensitivity and specificity below a desired value (default .55).

Only a single intersection is assumed (or an second intersection where the overlap is negligible).
}
\examples{
# A simple test model:
nlopt.ui()
# Using another bi-normal distribution:
nlopt.ui(mu0=0, sd0=1, mu1=1.6, sd1=2)
}
