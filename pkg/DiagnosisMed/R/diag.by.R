#' @rdname display
#' @export
diag.by <- function(ref, test, group.var, data, CL = .95, CL.type = c("wilson", "exact", "approximate"), group.var.labels = attr(data, "var.labels")[match(group.var, names(data))], test.var.labels = attr(data, "var.labels")[match(test, names(data))]){
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
  if (!is.character(group.var)) {
    stop("'group.var' should be a character vector.")
  }
  if (length(group.var) != 1) {
    stop("'group.var' should have length 1.")
  }
  if (!is.data.frame(data)) {
    stop("'data' should be a data.frame.")
  }
  if (!all(test %in% names(data))) {
    stop(paste0("The following variables of 'test' are not in 'data': ",test[!(test %in% names(data))], collapse = ","))
  }
  if (!is.null(group.var.labels)) {
    if (length(group.var.labels) != length(group.var)) {
      warning("'group.var' and 'group.var.labels' do not have the same lenth. Check carefully if labels match correctly.")
    }
  }
  if (!is.null(test.var.labels)) {
    if (length(test.var.labels) != length(test)) {
      warning("'test' and 'test.var.labels' do not have the same lenth. Check carefully if labels match correctly.")
    }
  }

  Group <- c(unlist(sapply(1:nlevels(data[,group.var]), function(i) rep(levels(data[,group.var])[i], length(test)))))
  Tests <- rep(test, nlevels(data[,group.var]))

  ################## COntinuar daqui
  output <- data.frame(Group, Tests, n = NA, prevalence = NA, Sensitivity = NA,
                       Se.inf.cl = NA, Se.sup.cl = NA, Specificity = NA, Sp.inf.cl = NA, Sp.sup.cl = NA,
                       PLR = NA, PLR.inf.cl = NA, PLR.sup.cl = NA)
  # i = 1
  for (i in 1:nrow(output)) {
    tab <- table(data[data[,group.var] == output$Group[i], as.character(output$Tests[i])], data[data[,group.var] == output$Group[i], as.character(ref)])
    if (all(dim(tab) == 2)) {
      output[i, 3:13] <- diagnosis(tab = tab, CL = CL, CL.type = CL.type[1])[c(2:3,6:14)]
    } else {
      warning(paste0("Not able to make a 2x2 table in the subset in the group '", as.character(output[i, "Group"]),"' with in test ", as.character(output[i, "Tests"])))
    }
  }
  output$Group <- as.character(output$Group)
  for (i in length(output$Group):2) {
    # i = 30
    if (output$Group[i - 1] == output$Group[i]) {output$Group[i] <- NA}
  }
  for (i in 1:nrow(output)) {
    if (all(is.na(output[i,3:13]))) {output <- output[-i,]}
  }
  if (!is.null(test.var.labels)) {
    output$Tests <- test.var.labels
  }
  if (!is.null(group.var.labels)) {
    output$Group[!is.na(output$Group)] <- paste0(group.var.labels, " = ", output$Group[!is.na(output$Group)])
  }
  class(output) <- c("data.frame", "diag.by")
  output
}

