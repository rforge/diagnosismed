#' Gwet's AC1 inter-rater statistic for 2 raters special case.
#'
#' @description Compute inter-rater or intra-rater agreement. Kappa, a common agreement statisitcs, assumes that agreement is at random and it's index express the agreement beyond the one observed at random. Gwet's AC1 statistic assumes that are agreement between observers are not totally at random such as Kappa. There will be some cases easy to agree in the condition absence, and some cases easy to agree in the condition presence and some that will be difficult to agree. Gwet also discuss that AC1 is way less prompt to a paradox that occurs with Kappa with exteme prevalences of the condition, where low Kappa are observed even with high crude agreement. TIn simulations studies, Gwet estimated bias for this statistic, therefore, for apropriate interpretation of AC1, it is necessary to subtract AC1 from  a critical value on table 6.16 on Gwet's book, and then use Fleis table as bencmarking. See \url{http://agreestat.com/}.
#'
#' @param tab k x k table which represents \code{\link[base]{table}}(rater1,rater2), must have equal number of rows and columns and have the same levels in the same order.
#'
#' @param conflev Confidence Level associated with the confidence interval (0.95 is the default value).
#'
#' @param N population size which will be stick in standard error correction, N = Inf is no correction.
#'
#' @param print Logical. Should results be printed on the console?
#'
#' @author Marcel Quintana and Pedro Brasil. Special thanks to Mr. Gwet for reviewing the code.
#'
#' @references Gwet. Computing inter-rater reliability and its variance in the presence of high agreement. The British Journal of Mathematical and Statistical Psychology. 61. 29-48. May 2008.
#'
#' @seealso package irr for many agreement statistics. epiR::epi.kappa for Kappa with confidence intervals. Also take a look at epiDisplay::kap for multi-rater multi-reader Kappa.
#'
#' @examples
#' table3 <- matrix(c(118, 2, 5, 0), nrow = 2, ncol = 2)
#'
#' AC1(table3)
#'
#' # The same thing as...
#' x <- AC1(table3, print = FALSE)
#' x
#'
#' @import stats
#' @export
AC1 <- function(tab, conflev = 0.95, N = Inf, print = TRUE){
  if(nrow(tab) != ncol(tab)){
    stop('The table should have the same number of rows and columns!')
  }
  if (conflev > 1 || conflev < 0) {
    stop('conflev argument should be between 0 and 1!')
  }
  if (!is.numeric(N) || N < 0) {
    stop('N must be a positive number!')
  }
  if (!is.logical(print)) {
    stop('print argument must be either TRUE or FALSE!')
  }

  # sifnificant part
  nc <- ncol(tab)
  nr <- nrow(tab)

  n <- sum(tab)
  f <- n / N
  pa <- sum(diag(tab)) / n # formula 18
  q <- ncol(tab) # number of categories
  pkk <- diag(tab) / n
  pak <- sapply(1:q, function(i) sum(tab[i , ])) / n
  pbk <- sapply(1:q, function(i) sum(tab[ , i])) / n
  pik <- (pak + pbk) / 2
  pegama <- (sum(pik * (1 - pik))) / (q - 1)
  gama <- (pa - pegama) / (1 - pegama) # AC1 statistics
  # 2 raters special case variance
  pkl <- tab / n
  soma <- 0
  for(k in 1:q){
	  for(l in 1:q){
		  soma <- soma + (pkl[k , l] * ((1 - (pik[k] + pik[l]) / 2) ^ 2))
	  }
  }
  vgama <- ((1 - f) / (n * (1 - pegama) ^ 2)) * (pa*(1 - pa) - 4 * (1 - gama) * ((1 / (q - 1)) * sum(pkk * (1 - pik)) - pa * pegama) + 4 * ((1 - gama) ^ 2) * ((1 / ((q - 1) ^ 2)) * soma - pegama ^ 2))
  epgama <- sqrt(vgama) # AC1 standard error
  lcb <- max(0, gama - epgama * qnorm(1 - (1 - conflev) / 2, 0, 1)) # lower confidence bound
  ucb <- min(1, gama + epgama * qnorm(1 - (1 - conflev) / 2, 0, 1)) # upper confidence bound
  names(pa) <- "Raw agreement"
  names(pegama) <- "Chance-independent agreement"
  names(gama) <- "Agreement coeficient (AC1)"
  names(epgama) <- "Standard error"
  names(lcb) <- "Lower CI"
  names(ucb) <- "Upper CI"
  if (print) {
    cat('Raw agreement:',pa,'; Chance-independent agreement:',pegama,'\n')
    cat('Agreement coefficient (AC1):',gama,'; AC1 standard error:',epgama,'\n')
    cat(conflev * 100,'% Confidence Interval: (',lcb,',',ucb,')\n')
  }
  invisible(c(pa,pegama,gama,epgama,lcb,ucb))
}
# For apropriate interpretation of AC is necessary to subtract AC1
# from  a critical value on table 6.16 on gwet's book and then use Fleis table as bencmarking
#nc <- ncol(table)
#nr <- nrow(table)
#cm <- sapply(1:q,function(i)sum(table[,i]))# colum marginals
#rm <- sapply(1:q,function(i)sum(table[i,]))# row marginals
#rt <- matrix(NA,nc,nr)# random table from data
#for(c in 1:nc){
#  for(r in 1:nr){
#    rt[r,c] <- (cm[c] + rm[r])/2
#  }
#}

