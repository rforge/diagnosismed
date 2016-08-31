#' @rdname TGROC
#' @export
print.TGROC <- function(x, ..., digits = 3, print.type = "all"){
  if(print.type != "all" && print.type != "NN-parametric" && print.type != "Non-parametric" && type != "Binormal"){
    stop("print.type argument in print.TGROC is not valid.")
  }
  cat("                Sample size:",x$sample.size,"\n")
  cat("          Sample prevalence:",round(x$sample.prevalence, digits = digits),"\n")
  cat("      Population prevalence:",round(x$pop.prevalence, digits = digits)," - same as sample prevalence if not informed\n")
  cat("      Informed cost - FP/FN:",round(x$cost, digits = digits),"\n")
  cat("Informed inconclusive level:",round(x$inc, digits = digits),"\n")
  cat("\n\n")
  cat("Test accuracy ------------------------------------------------------------\n")
  print(x$AUC, ... , digits = digits)
  cat("\n\n")
  cat("Test summary ------------------------------------------------------------\n")
  print(x$test.summary, ... , digits = digits)
  cat("\n\n")
  if(any(print.type %in% c("all" , "Non-parametric"))){
    cat("Non-paramentric inconclusive threshold limits with",x$inc,"inconclusive tolerance -----------------------------------------------\n")
    print(x$np.inconclusive, ... , digits = digits)
    cat("\n\n")
    cat("Non-paramentric best threshold estimations with",x$conf.limit,"confidence limits ------------------------------------------------------\n")
    print(x$np.best.threshold, ... , digits = digits)
    cat("\n\n")
  }
  if(any(print.type %in% c("all" , "NN-parametric"))){
    cat("Paramentric (neural network) inconclusive threshold limits with",x$inc,"inconclusive tolerance ----------------------------------\n")
    print(x$NN.inconclusive, ... , digits = digits)
    cat("\n\n")
    cat("Paramentric (neural network) best threshold estimations with",x$conf.limit,"confidence limits -----------------------------------------\n")
    print(x$NN.best.threshold, ... , digits = digits)
    cat("\n\n")
  }
  if(any(print.type %in% c("all" , "Binormal"))){
    cat("Paramentric (Binormal) inconclusive threshold limits with",x$inc,"inconclusive tolerance ----------------------------------\n")
    print(x$BN.inconclusive, ... , digits = digits)
    cat("\n\n")
    cat("Paramentric (Binormal) best threshold estimations with",x$conf.limit,"confidence limits -----------------------------------------\n")
    print(x$BN.best.threshold, ... , digits = digits)
    cat("\n\n")
  }
}
