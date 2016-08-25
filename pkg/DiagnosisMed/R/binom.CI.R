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
