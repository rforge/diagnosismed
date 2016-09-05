#' get.intersection Obtain the intersection of two distributions using the kernel method
#'
#' @name get.intersection
#' @description \code{get.intersection} Obtain the intersection of two distributions using the kernel method. Warning: This function does not check the parameters ref and test.
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results on a continuous scale.
#' @param ... passing arguments to the kernel density function, other than kernel='gaussian' (default).
#' @return A vector of points of intersection, ordered on their density. The tail has the highest density.
#' @seealso \code{\link{density}}
#'
#' @examples
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,2))
#' get.intersection(ref, test)
#' @export

get.intersection <- function(ref, test, ...){
  
  # determine points of intersection 
  x0 = test[ref == 0]
  x1 = test[ref == 1]
  lower = min(test)
  upper = max(test)
  # generate kernel densities
  d0 <- stats::density(x0, from = lower, to = upper, n = 2048, ...)
  d1 <- stats::density(x1, from = lower, to = upper, n = 2048, ...)
  dif = as.logical(abs(diff(d0$y < d1$y)))
  d.h0 = d0$y[dif]
  intersection <- d1$x[dif]
  o=order(d.h0)
  intersection=intersection[o] # order on density; tail has highest density
}


