#' get.intersection Obtain the intersection of two distributions using the kernel method
#'
#' @name get.intersection
#' @description \code{get.intersection} Obtain the intersection of two distributions using the kernel method. Warning: This function does not check the parameters ref and test.
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results on a continuous scale.
#' @param model The model used for estimating the intersection(s). Default = 'kernel'.
#' @param ... passing arguments to the kernel density function, other than kernel='gaussian' (default).
#' @return A vector of points of intersection, ordered on their density. The tail has the highest density.
#' @references Landsheer, J. A. (2016). Interval of Uncertainty: An Alternative Approach for the Determination of Decision Thresholds, with an Illustrative Application for the Prediction of Prostate Cancer. PloS One, 11(11), e0166007.
#' @seealso \code{\link{density}}
#'
#' @examples
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,2))
#' (get.intersection(ref, test)) # two intersections! Generates warning in other functions!
#' @export

get.intersection <-
  function(ref,
           test,
           model = c('kernel', 'binormal', 'ordinal'),
           ...) {

    model <- match.arg(model)
    df = check.data(ref, test, model=model)

    intersect.binormal1 <-
      function(mu0, sd0, mu1, sd1) {
        if (sd0==sd1){
          is <- (mu1+mu0)/2
        }else{
          B <- (mu0 / sd0 ^ 2 - mu1 / sd1 ^ 2)
          A <- 0.5 * (1 / sd1 ^ 2 - 1 / sd0 ^ 2)
          C <-
            0.5 * (mu1 ^ 2 / sd1 ^ 2 - mu0 ^ 2 / sd0 ^ 2) - log(sd0 / sd1)
          
          is = (-B + c(1, -1) * sqrt(B ^ 2 - 4 * A * C)) / (2 * A)
        }
        d = dnorm(is, mu0, sd0) + dnorm(is, mu1, sd1)
        is[order(d)] # tail has highest density
      }
    # determine points of intersection
    x0 = test[ref == 0]
    x1 = test[ref == 1]

    if (model == 'kernel') {
      lower = min(test)
      upper = max(test)
      # generate kernel densities
      d0 <- stats::density(x0,
                           from = lower,
                           to = upper,
                           n = 2048,
                           ...)
      d1 <- stats::density(x1,
                           from = lower,
                           to = upper,
                           n = 2048,
                           ...)
      dif = as.logical(abs(diff(d0$y < d1$y)))
      d.h0 = d0$y[dif] + d1$y[dif]
      intersection <- d1$x[dif]
      o = order(d.h0)
      return(intersection[o]) # order on density; tail has highest density
    } else if (model == 'ordinal') {
        lower = min(test)
        upper = max(test)
        # generate kernel densities
        if (any(names(list(...))== 'adjust')) {
          d0 <- stats::density(x0,
                               from = lower,
                               to = upper,
                               n = 2048,
                               ...)
          d1 <- stats::density(x1,
                               from = lower,
                               to = upper,
                               n = 2048,
                               ...)
        } else {
          d0 <- stats::density(x0,
                               from = lower,
                               to = upper,
                               n = 2048,
                               adjust=2,
                               ...)
          d1 <- stats::density(x1,
                               from = lower,
                               to = upper,
                               n = 2048,
                               adjust=2,
                               ...)

        }
        dif = as.logical(abs(diff(d0$y < d1$y)))
        d.h0 = d0$y[dif] + d1$y[dif]
        intersection <- d1$x[dif]
        o = order(d.h0)
        return(intersection[o]) # order on density; tail has highest density
      } else if (model == 'binormal') {

        mu0 = mean(df$test[df$ref == 0])
        sd0 = sd(df$test[df$ref == 0])
        mu1 = mean(df$test[df$ref == 1])
        sd1 = sd(df$test[df$ref == 1])
        return(intersect.binormal1(mu0, sd0, mu1, sd1))

      }
    }



