#' @rdname display
#' @export
diag.stack <- function(ref, test, data, CL = .95, CL.type = c("wilson", "exact", "approximate"),
                        var.labels = attr(data, "var.labels")[match(test, names(data))]){
  if (!is.character(ref)) {
    stop("'ref' should be a character vector.")
  }
  if (length(ref) != 1) {
    stop("'ref' should have length 1.")
  }
  if (!is.character(test)) {
    stop("'test' should be a character vector.")
  }
  if (length(test) < 2) {
    stop("'test' should have at least length 2.")
  }
  if (!is.data.frame(data)) {
    stop("'data' should be a data.frame.")
  }
  if (!all(test %in% names(data))) {
    stop(paste0("The following variables of 'tests' are not in 'data': ",test[!(test %in% names(data))], collapse = ","))
  }
  if (!is.null(var.labels)) {
    if (length(var.labels) != length(test)) {
      warning("'test' and 'var.labels' do not have the same lenth. Check carefully if labels match correctly.")
    }
  }
  output <- data.frame(Tests = test, n = NA, prevalence = NA, Sensitivity = NA,
                       Se.inf.cl = NA, Se.sup.cl = NA, Specificity = NA, Sp.inf.cl = NA, Sp.sup.cl = NA,
                       PLR = NA, PLR.inf.cl = NA, PLR.sup.cl = NA)
  for (i in test) {
    # i = "BAAR1"
    tab <- table(data[, i], data[, ref])
    if (all(dim(tab) == 2)) {
      output[output$Tests == i, 2:12] <- diagnosis(tab = tab, CL = CL, CL.type = CL.type[1])[c(2:3,6:14)]
    } else {
      warning(paste0("Not able to make a 2 x 2 table with the variable ", i))
    }
  }
  output$Tests <- as.character(output$Tests)
  for (i in output$Tests) {
    if (all(is.na(output[output$Tests == i, 2:12]))) {output <- output[-which(output$Tests == i),]}
  }
  if (!is.null(var.labels)) {
    output$Tests <- var.labels
  }
  class(output) <- c("data.frame", "diag.stack")
  output
}

