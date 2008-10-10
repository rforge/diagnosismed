print.diag <- function(x,...){
  #options(digits=4)
  cat("\n\n")
  print(x$tabmarg)
  cat("\n\n")
  cat("The test has the following parameters [",paste(x$Conf.limit*100),"% confidence interval]\n")
  cat("---------------------------------------------------------------\n")
  cat("Sample size:                ",x$n,"\n")
  cat("Prevalence considered:      ",round(x$p*100,digits=2),"%\n")
  cat("Sensitivity:                ",round(x$Se*100,digits=2),"% [",round(x$Se.cl[1]*100,digits=2)," - ",round(x$Se.cl[2]*100,digits=2),"]\n")
  cat("Specificity:                ",round(x$Sp*100,digits=2),"% [",round(x$Sp.cl[1]*100,digits=2)," - ",round(x$Sp.cl[2]*100,digits=2),"]\n")
  cat("Positive predictive value:  ",round(x$PPV*100,digits=2),"% [",round(x$PPV.cl[1]*100,digits=2)," - ",round(x$PPV.cl[2]*100,digits=2),"]\n")
  cat("Negative predictive value:  ",round(x$NPV*100,digits=2),"% [",round(x$NPV.cl[1]*100,digits=2)," - ",round(x$NPV.cl[2]*100,digits=2),"]\n")
  cat("Positive likelihood ratio:  ",round(x$PLR,digits=2),"  [",round(x$PLR.inf.cl,digits=2)," - ",round(x$PLR.sup.cl,digits=2),"]\n")
  cat("Negative likelihood ratio:  ",round(x$NLR,digits=2),"   [",round(x$NLR.inf.cl,digits=2)," - ",round(x$NLR.sup.cl,digits=2),"]\n")
  cat("Diagnostic odds ratio:      ",round(x$DOR,digits=2)," [",round(x$DOR.inf.cl,digits=2)," - ",round(x$DOR.sup.cl,digits=2),"]\n")
  cat("Error trade off (FN : FP)   ",round(x$ET,digits=2)," : 1 \n")
  cat("Error rate:                 ",round(x$ER*100,digits=2),"% [",round(x$ER.cl[1]*100,digits=2)," - ",round(x$ER.cl[2]*100,digits=2),"]\n")
  cat("Accuracy:                   ",round(x$accu*100,digits=2),"% [",round(x$accu.cl[1]*100,digits=2)," - ",round(x$accu.cl[2]*100,digits=2),"]\n")
  cat("Youden index:               ",round(x$Youden,digits=4)," [",round(x$Youden.inf.cl,digits=4)," - ",round(x$Youden.sup.cl,digits=4),"]\n")
  cat("Area under ROC curve:       ",round(x$AUC,digits=4),"\n")
  cat("---------------------------------------------------------------\n")
  #options(digits=7)
}
