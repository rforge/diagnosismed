#' @export
print.diag <- function(x,...){
  cat("Reference standard:",x$reference.name,"\n")
  cat("Index test        :",x$index.name,"\n")  
  cat("---------------------------------------------------------------\n")
  print(x$tabmarg)
  cat("\n")
  cat("The test has the following parameters [",x$Conf.limit*100,"% confidence interval]\n",sep="")
  cat("---------------------------------------------------------------\n")
  cat("Sample size:                  ",x$n,"\n")
  cat("Prevalence considered(%):     ",formatC(x$p*100,digits=2,format="f"),"\n")
  cat("Sensitivity(%):               ",formatC(x$Se*100,digits=2,format="f")," [",formatC(x$Se.cl[1]*100,digits=2,format="f")," - ",formatC(x$Se.cl[2]*100,digits=2,format="f"),"]\n")
  cat("Specificity(%):               ",formatC(x$Sp*100,digits=2,format="f")," [",formatC(x$Sp.cl[1]*100,digits=2,format="f")," - ",formatC(x$Sp.cl[2]*100,digits=2,format="f"),"]\n")
  cat("Positive predictive value(%): ",formatC(x$PPV*100,digits=2,format="f")," [",formatC(x$PPV.cl[1]*100,digits=2,format="f")," - ",formatC(x$PPV.cl[2]*100,digits=2,format="f"),"]\n")
  cat("Negative predictive value:(%):",formatC(x$NPV*100,digits=2,format="f")," [",formatC(x$NPV.cl[1]*100,digits=2,format="f")," - ",formatC(x$NPV.cl[2]*100,digits=2,format="f"),"]\n")
  cat("Positive likelihood ratio:    ",formatC(x$PLR,digits=2,format="f"),"  [",formatC(x$PLR.inf.cl,digits=2,format="f")," - ",formatC(x$PLR.sup.cl,digits=2,format="f"),"]\n")
  cat("Negative likelihood ratio:    ",formatC(x$NLR,digits=2,format="f"),"  [",formatC(x$NLR.inf.cl,digits=2,format="f")," - ",formatC(x$NLR.sup.cl,digits=2,format="f"),"]\n")
  cat("Diagnostic odds ratio:        ",formatC(x$DOR,digits=2,format="f")," [",formatC(x$DOR.inf.cl,digits=2,format="f")," - ",formatC(x$DOR.sup.cl,digits=2,format="f"),"]\n")
  cat("Error trade off (FN : FP)      ",round(x$ET,digits=2)," : 1 \n",sep='')
  cat("Error rate(%):                ",formatC(x$ER*100,digits=2,format="f")," [",formatC(x$ER.cl[1]*100,digits=2,format="f")," - ",formatC(x$ER.cl[2]*100,digits=2,format="f"),"]\n")
  cat("Accuracy(%):                  ",formatC(x$accu*100,digits=2,format="f")," [",formatC(x$accu.cl[1]*100,digits=2,format="f")," - ",formatC(x$accu.cl[2]*100,digits=2,format="f"),"]\n")
  cat("Youden index:                 ",formatC(x$Youden,digits=4,format="f")," [",formatC(x$Youden.inf.cl,digits=4,format="f")," - ",formatC(x$Youden.sup.cl,digits=4,format="f"),"]\n")
  cat("Area under ROC curve:         ",round(x$AUC,digits=4),"\n")
  cat("---------------------------------------------------------------\n")
}