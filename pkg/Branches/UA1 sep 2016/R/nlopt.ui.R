#' Function for the determination of the population thresholds an inconclusive interval for bi-normal distributed test scores.
#'
#' @param Se (default = .55). Desired sensitivity of the test scores within the uncertain interval. A value below .5 is not allowed, while a value larger than .6 is not recommended.
#' @param Sp (default = .55). Desired specificity of the test scores within the uncertain interval. A value below .5 is not allowed, while a value larger than .6 is not recommended.
#' @param mu0 Population value or estimate of the mean of the test scores of the persons without the targeted condition.
#' @param sd0 Population value or estimate of the standard deviation of the test scores of the persons without the targeted condition.
#' @param mu1 Population value or estimate of the mean of the test scores of the persons with the targeted condition.
#' @param sd1 Population value or estimate of the standard deviation of the test scores of the persons with the targeted condition.
#' @param intersection Default NULL. If not null, the supplied value is used as the estimate of the intersection of the two bi-normal distributions. Otherwise, it is calculated.
#' @param start Default NULL. If not null, the first two values of the supplied vector are used as the starting values for the \code{nloptr} optimization function.
#' @param print.level Default is 0. The option print_level controls how much output is shown during the optimization process. Possible values: 0 (default)	no output; 1	show iteration number and value of objective function; 2	1 + show value of (in)equalities; 3	2 + show value of controls.
#'
#' @details The function can be used to determinate the uncertain interval of two bi-normal distributions.
#' The Uncertain Interval is defined as an interval below and above the intersection of the two distributions, with a sensitivity and specificity below a desired value (default .55).
#'
#' Only a single intersection is assumed (or an second intersection where the overlap is negligible).
#'
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
#' @importFrom rootSolve uniroot.all
#' @importFrom nloptr nloptr
#' @importFrom stats dnorm pnorm sd
#'
#' @examples
#' # A simple test model:
#' nlopt.ui()
#' # Using another bi-normal distribution:
#' nlopt.ui(mu0=0, sd0=1, mu1=1.6, sd1=2)
nlopt.ui <- function(Se = .55, Sp = .55,
                     mu0 = 0, sd0 = 1,
                     mu1 = 1, sd1 = 1,
                     intersection = NULL,
                     start=NULL, print.level=0) {
  c01 = Sp / (1 - Sp)
  c11 = Se / (1 - Se)
  if (is.null(intersection)) {
    f <- function(x)
      dnorm(x, mean = mu0, sd = sd0) - dnorm(x, mean = mu1, sd = sd1)

    I = uniroot.all(f, interval = c(max(mu0 - 3 * sd0, mu1 - 3 * sd1),
                                    min(mu0 + 3 * sd0, mu1 + 3 * sd1)))
    if (length(I) > 1) {
      d = dnorm(I, mu0, sd0)
      I = I[which.max(d)]
      warning('More than one point of intersection. Point with highest density used.')
    }
  } else I = intersection

  # objective function -(H - L)^2
  eval_f0 <- function(x) {
    return(-(x[2] - x[1]) ^ 2)
    # return(-(x[2] - x[1]) )
  }

  eval_grad_f0 <- function(x) {
    return(c(2 * (x[2] - x[1]),-2 * (x[2] - x[1])))
    # return(c(1,-1))
  }

  # constraint function
  eval_g0 <- function(x) {
    a0 = pnorm(I, mu0, sd0) - pnorm(x[1], mu0, sd0) # TN
    b0 = pnorm(x[2], mu0, sd0) - pnorm(I, mu0, sd0) # FP
    a1 = pnorm(x[2], mu1, sd1) - pnorm(I, mu1, sd1) # TP
    b1 = pnorm(I, mu1, sd1) - pnorm(x[1], mu1, sd1) # FN

    return(c(a0 - c01 * b0,
             a1 - c11 * b1)) # vector with two constraint values
  }

  # jacobian of constraints
  eval_jac_g0 <- function(x) {
    return(rbind(c(
      -dnorm(x[1], mu0, sd0) , -c01 * dnorm(x[2], mu0, sd0)
    ),
    c(
      c11 * dnorm(x[1], mu1, sd1), dnorm(x[2], mu1, sd1)
    )))
  }

  if (is.null(start)) par0=c(L=I-.5*sd1,
                             H=I+.5*sd0) else par0=c(L=start[1], H=start[2])

  res0 <- nloptr( x0=par0,
                  eval_f=eval_f0,
                  lb = c(I-2*sd1, I),
                  ub = c(I, I+2*sd0),
                  eval_g_ineq = eval_g0,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-8,
                              print_level=print.level)
  )

  # res0 <- nloptr(
  #   x0 = c(I - sd1, I + sd0),
  #   eval_f = eval_f0,
  #   eval_grad_f = eval_grad_f0,
  #   eval_g_ineq = eval_g0,
  #   eval_jac_g_ineq = eval_jac_g0,
  #   lb = c(I-2*sd1, I),
  #   ub = c(I, I+2*sd0),
  #   opts = list(
  #     "algorithm" = "NLOPT_LD_MMA",
  #     "xtol_rel" = 1.0e-8,
  #     "print_level" = 2 #,
  #     # "check_derivatives" = TRUE,
  #     # "check_derivatives_print" = "all"
  #     )
  #   )

  TN = pnorm(I, mu0, sd0) - pnorm(res0$solution[1], mu0, sd0)  # area check Sp: lower area / upper area
  FP = pnorm(res0$solution[2], mu0, sd0) - pnorm(I, mu0, sd0)
  TP = pnorm(res0$solution[2], mu1, sd1) - pnorm(I, mu1, sd1)  # area check Se: upper area / lower area
  FN = pnorm(I, mu1, sd1) - pnorm(res0$solution[1], mu1, sd1)
  res = list()
  res$status = res0$status
  res$message = res0$message
  res$intersection = I
  res$results = c(exp.Sp.ui = ifelse((TN > 1e-4),TN / (FP + TN), NA),
                  exp.Se.ui = ifelse(TP > 1e-4, TP / (FN + TP), NA),
                  mu0=unname(mu0), sd0=unname(sd0), mu1=unname(mu1), sd1=unname(sd1))
  res$solution = c(L = res0$solution[1], U = res0$solution[2])

  # TN.FP=unname((pnorm(I, mu0, sd0)-pnorm(res$solution['L'], mu0, sd0)) / # area check Sp: lower area / upper area
  #                (pnorm(res$solution['H'], mu0, sd0)-pnorm(I, mu0, sd0)))
  # TP.FN=unname((pnorm(res$solution['H'], mu1, sd1)-pnorm(I, mu1, sd1)) / # area check Se: upper area / lower area
  #                (pnorm(I, mu1, sd1)-pnorm(res$solution['L'], mu1, sd1)))
  # res$results2=c(Sp=TN.FP/(1+TN.FP), Se=TP.FN/(1+TP.FN))
  return(res)
}

