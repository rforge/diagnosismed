#' @rdname diagnosis
#' @export
print.diag <- function(x, digits = 3, ...){
  cat("Reference standard:",x$reference.name,"\n")
  cat("Index test        :",x$index.name,"\n")
  cat("---------------------------------------------------------------\n")
  print(x$tabmarg)
  cat("\n")
  cat("The test has the following parameters [",x$Conf.limit * 100,"% confidence interval]\n", sep = "")
  cat("---------------------------------------------------------------\n")
  print(x$results, digits = digits, ...)
  cat("---------------------------------------------------------------\n")
}
