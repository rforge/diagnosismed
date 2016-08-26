#' Confidence intervals for binomial counts or proportions
#'
#' @description Calculates confidence intervals for binomial counts or proportions
#'
#' @param x Number of successes in n trials, can be a vector.
#' @param n Number of Bernoulli trials, can be a vector.
#' @param conf.level Confidence level (default = 0.95), can be a vector.
#' @param type Possible values are "wilson", "exact" and "approximate". See datails.
#'
#' @details \code{binom.CI} with type "exact", calculates exact confidence intervals for binomial counts or proportions. This function uses R's binom.test function; however, the arguments to this function can be numeric vectors of any length.
#' \code{binom.CI} with type "wilson", calculates confidence intervals for binomial counts or proportions using Wilson's formula which approximate the exact method. The arguments to this function can be numeric vectors of any length.
#' \code{binom.CI} with type "approximate" calculates confidence intervals for binomial counts or proportions using a normal approximation to the binomial distribution. The arguments to this function can be numeric vectors of any length.
#' This function is a clone from epitools binomial confidence intervals functions, but with a wrapper to include the method as an option instead of separate functions.
#'
#' @return This function returns a n x 6 matrix with the following colnames:
#' x number of successes in n trials
#' n number of Bernoulli trials
#' prop proportion = x / n
#' lower lower confidence interval limit
#' upper upper confidence interval limit
#' conf.level confidence level
#'
#' @references
#' Tomas Aragon, et al. Applied Epidemiology Using R.
#' Kenneth Rothman (2002), Epidemiology: An Introduction, Oxford University Press, 1st Edition.
#'
#' @seealso pois.exact, binom.test

# Examples
#
# binom.exact(1:10, seq(10, 100, 10))
# binom.wilson(1:10, seq(10, 100, 10))
# binom.approx(1:10, seq(10, 100, 10))

binom.CI <- function(x, n, conf.level = 0.95, type = c("wilson","exact","approximate")) {
# warning section
  if (conf.level > 1 || conf.level < 0) {
    stop("conf.level argument must between o and 1")
  }
  if (!(any(type %in% c("wilson","exact","approximate")))){
    stop("type argument must be either 'wilson','exact', or 'approximate'")
  }
  binom.wilson <- function(x, n, conf.level = conf.level) {
    Z <- qnorm(0.5 * (1 + conf.level))
    Zinsert <- Z * sqrt(((x * (n - x)) / n ^ 3) + Z ^ 2 / (4 * n ^ 2))
    R.lower <- (n / (n + Z ^ 2)) * (x / n + Z ^ 2 / (2 * n) - Zinsert)
    R.upper <- (n/(n + Z ^ 2)) * (x / n + Z ^ 2 / (2 * n) + Zinsert)
    output <- data.frame(x = x, n = n, proportion = x / n, lower = R.lower, upper = R.upper, conf.level = conf.level)
    output
  }
  binom.exact <- function(x, n, conf.level = conf.level) {
    xnc <- cbind(x, n, conf.level)
    lower <- numeric(nrow(xnc))
    upper <- numeric(nrow(xnc))
    for (i in 1:nrow(xnc)) {
      ci <- binom.test(x = xnc[i, 1], n = xnc[i, 2], conf.level = xnc[i, 3])$conf.int
      lower[i] <- ci[1]
      upper[i] <- ci[2]
    }
    output <- data.frame(x = x, n = n, proportion = x/n, lower = lower, upper = upper, conf.level = conf.level)
    output
  }
  binom.approx <- function(x, n, conf.level = conf.level) {
    Z <- qnorm(0.5 * (1 + conf.level))
    SE.R <- sqrt(x * (n - x) / (n ^ 3))
    R.lci <- x / n - Z * SE.R
    R.uci <- x / n + Z * SE.R
    output <- data.frame(x = x, n = n, proportion = x/n, lower = R.lci, upper = R.uci, conf.level = conf.level)
    output
  }
  if (type[1] == "wilson") {
    output <- binom.wilson(x = x, n = n, conf.level = conf.level)
  }
  if (type[1] == "exact") {
    output <- binom.exact(x = x, n = n, conf.level = conf.level)
  }
  if (type[1] == "approximate") {
    output <- binom.approximate(x = x, n = n, conf.level = conf.level)
  }
  output
}

# binom.CI(100, 500)
