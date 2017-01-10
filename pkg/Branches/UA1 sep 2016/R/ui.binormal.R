#' Function for the determination of the sample thresholds an inconclusive interval for bi-normal distributed test scores.
#'
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition).
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param Se (default = .55). Desired sensitivity of the test scores within the uncertain interval. A value below .5 is not allowed, while a value larger than .6 is not recommended.
#' @param Sp (default = .55). Desired specificity of the test scores within the uncertain interval. A value below .5 is not allowed, while a value larger than .6 is not recommended.
#' @param intersection Default NULL. If not null, the supplied value is used as the estimate of the intersection of the two bi-normal distributions. Otherwise, it is calculated.
#' @param start Default NULL. If not null, the first two values of the supplied vector are used as the starting values for the \code{nloptr} optimization function.
#' @param print.level Default is 0. The option print_level controls how much output is shown during the optimization process. Possible values: 0 (default)	no output; 1	show iteration number and value of objective function; 2	1 + show value of (in)equalities; 3	2 + show value of controls.
#'
#' @details{
#' This function can be used for a test with bi-normal distributed scores.
#' The Uncertain Interval is defined as an interval below and above the intersection, with a sensitivity and specificity below a desired value (default .55).
#'
#' Only a single intersection is assumed (or an second intersection where the overlap is negligible). If another intersection exists and the overlap around this intersection is considerable, the test with such a non-negligible overlap is problematic and difficult to apply and interpret.
#'
#' In general, when estimating decision thresholds, a sample of sufficient size should be used. It is recommended to use at least a sample of 100 patients with the targeted condition, and a 'healthy' sample (without the targeted condition) of the same size or larger.
#'
#' Lastly, it should be noted that the Uncertain interval method has been developed recently, and future research may provide more satisfactory answers.
#' }
#' @return List of values:
#' \describe{
#'   \item{$status: }{Integer value with the status of the optimization (0 is success).}
#'   \item{$message: }{More informative message with the status of the optimization}
#'   \item{$results: }{Vector with the following values:}
#'      \itemize{
#'       \item{exp.Sp.ui: }{The population value of the specificity in the Uncertain Interval, given mu0, sd0, mu1 and sd1. This value should be very near the supplied value of Sp.}
#'       \item{exp.Sp.ui: }{The population value of the sensitivity in the Uncertain Interval, given mu0, sd0, mu1 and sd1. This value should be very near the supplied value of Se.}
#'       \item{mu0: }{The value that has been supplied for mu0.}
#'       \item{sd0: }{The value that has been supplied for sd0.}
#'       \item{mu1: }{The value that has been supplied for mu1.}
#'       \item{sd1: }{The value that has been supplied for sd1.}
#'       }
#'   \item{$solution: }{Vector with the following values:}
#'      \itemize{
#'       \item{L: }{The population value of the lower threshold of the Uncertain Interval.}
#'       \item{U: }{The population value of the upper threshold of the Uncertain Interval.}
#'       }
#' }
#' @export
#' @examples
#' # A simple test model
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1))
#' ui.binormal(ref, test)
#'

ui.binormal <- function(ref, test, Se = .55, Sp = .55,
                        intersection = NULL, start=NULL, print.level=0){
  df = check.data(ref, test)
  n0 = df$test[ref==0]
  n1 = df$test[ref==1]
  m0 = mean(n0)
  sd0 = sd(n0)
  m1 = mean(n1)
  sd1 = sd(n1)
  # res0=fitdistr(n0,"normal") # MASS
  # m0=res0$estimate[1]; sd0=res0$estimate[2]
  # res1=fitdistr(n1,"normal") # fitdistrplus
  # m1=res1$estimate[1]; sd1=res1$estimate[2]



  # Se=.55; Sp=.55;
  # mu0 = m0; sd0 = sd0;
  # mu1 = m1; sd1 = sd1;
  # intersection = NULL; start=NULL; print.level=0
  res=nlopt.ui(Se = Se, Sp = Sp,
               mu0 = m0, sd0 = sd0,
               mu1 = m1, sd1 = sd1,
               intersection = intersection,
               start=start, print.level=print.level)
  return(res)
}
