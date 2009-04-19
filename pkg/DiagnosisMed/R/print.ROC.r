print.ROC<-function(x,Full=FALSE,...){
  if (Full==TRUE){ page(x$test.diag.table,method="print")}
  cat("          Sample size:",paste(x$sample.size),"\n")
  cat("    Sample prevalence:",paste(round(x$sample.prevalence,digits = 4)),"\n")
  cat("Population prevalence:",paste(round(x$pop.prevalence,digits = 4))," - same as sample prevalence if not informed\n")
  cat("\n\n")
  cat("Non-parametric AUC (trapezoidal method) and its confidence",x$CL," limits (DeLong method)\n")  
  cat(" Area under ROC curve:",paste(round(x$AUC.summary[2],digits = 4)),"[",paste(round(x$AUC.summary[1],digits = 4)),"-",
        paste(round(x$AUC.summary[3],digits = 4)),"]\n")
  cat("\n\n")
  cat("Test summary-----------------------------------------------------\n")  
  print(x$test.summary)
  cat("\n\n")
  cat("Best cut-off estimations with",x$CL,"confidence limits -----------\n")  
  print(x$test.best.cutoff)
}  
