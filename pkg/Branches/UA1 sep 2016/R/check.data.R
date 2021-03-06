#' Function to check the dataset of individuals with (1) and without (0) the targeted condition.
#'
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param model The model used for estimation. Default = 'kernel'. When model != ordinal, the dataset is checked whether the test has a sufficient number of different values (>= 20).
#' @details The first check is whether ref and test have equal length. If not, checkdata is aborted with an error message.
#' The second check is whether ref is coded solely with 0 and 1. If not, check.data is aborted and an error message is shown.
#' The third check is whether ref and test have missing values. If true, list wise deletion is applied and a warning message is shown.
#' The fourth check is whether test is continuous or not. If test has less than 20 different values, a warning message is shown. This test is omitted when ordinal = TRUE.
#'
#' This function is called from every function that requires data. A call is only useful to check warnings and errors.
#'
#' @return Either a valid dataset as data.frame with two variables ref and test or an error message.
#' @export
#'
#' @examples
#' #' set.seed(1)
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1.2))
#' check.data(ref, test) # model = 'kernel'

check.data <- function(ref, test, model = c('kernel', 'binormal', 'ordinal')){

  model <- match.arg(model)

   if (length(ref) != length(test)) stop('parameters ref and test have unequal length')

   if (any(levels(as.factor(ref)) != c(0, 1))) {
     stop("Your reference standard must be coded as 0 (absence) and 1 (presence). Check reference categories!")
   }
  sel = stats::complete.cases(ref, test) # sum(sel)
  ref=ref[sel]
  test=test[sel]
  if (!all(sel)) warning('The data has missing values. List wise deletion has been applied, but multiple imputation may provide better results.')
   t=table(test)
   if (model != 'ordinal') {
     if (length(names(t)) < 20) warning('The ui.nonpar and ui.binormal functions have been developed for continuous data. Your test has less than 20 different values.')
   }
   return(data.frame(ref=ref, test=test))
}

# ref[1000]=NA
# test=round(test)
