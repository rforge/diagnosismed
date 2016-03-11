#' @export
print.ROC<-function(x,Full=FALSE,...){
  ## View cannot be used in examples !!!
  if (Full==TRUE){ View(x$test.diag.table[,-c(2:5,24:34)])}
  cat("          Sample size:",x$sample.size,"\n")
  cat("    Sample prevalence:",round(x$sample.prevalence,digits = 4),"\n")
  cat("Population prevalence:",round(x$pop.prevalence,digits = 4)," - same as sample prevalence if not informed\n")
  cat("Informed Cost - cost(FN)/cost(FP):",x$cost,"\n")  
  cat("\n\n")
  cat("Non-parametric AUC (trapezoidal method) and its",x$CL,"confidence limits (DeLong method)\n")  
  cat(" Area under ROC curve:",paste(round(x$AUC.summary[2],digits = 4)),"[",paste(round(x$AUC.summary[1],digits = 4)),"-",
        paste(round(x$AUC.summary[3],digits = 4)),"]\n")
  cat("\n\n")
  cat("Test summary-----------------------------------------------------\n")  
  print(x$test.summary)
  cat("\n\n")
  cat("Best cut-off estimations with",x$CL,"confidence limits -----------\n")  
  print(x$test.best.cutoff)
}  
