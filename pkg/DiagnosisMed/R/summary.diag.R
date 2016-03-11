#' @export
summary.diag <- function(object,...){
  diag.tab <- matrix(
        c(object$n,NA,paste(formatC(object$p*100,digits=2,format="f")),NA,
          formatC(object$Se*100,digits=2,format="f"),
                     paste('[',formatC(object$Se.cl[1]*100,digits=2,format="f"),'-',formatC(object$Se.cl[2]*100,digits=2,format="f"),']'),
          formatC(object$Sp*100,digits=2,format="f"),
                     paste('[',formatC(object$Sp.cl[1]*100,digits=2,format="f"),'-',formatC(object$Sp.cl[2]*100,digits=2,format="f"),']'),
          formatC(object$PPV*100,digits=2,format="f"),
                     paste('[',formatC(object$PPV.cl[1]*100,digits=2,format="f"),'-',formatC(object$PPV.cl[2]*100,digits=2,format="f"),']'),
          formatC(object$NPV*100,digits=2,format="f"),
                     paste('[',formatC(object$NPV.cl[1]*100,digits=2,format="f"),'-',formatC(object$NPV.cl[2]*100,digits=2,format="f"),']'),
          formatC(object$PLR,digits=2,format="f"),paste('[',formatC(object$PLR.inf.cl,digits=2,format="f"),'-',formatC(object$PLR.sup.cl,digits=2,format="f"),']'),
          formatC(object$NLR,digits=2,format="f"),paste('[',formatC(object$NLR.inf.cl,digits=2,format="f"),'-',formatC(object$NLR.sup.cl,digits=2,format="f"),']'),
          formatC(object$DOR,digits=2,format="f"),paste('[',formatC(object$DOR.inf.cl,digits=2,format="f"),'-',formatC(object$DOR.sup.cl,digits=2,format="f"),']'),
          paste(round(object$ET,digits=2),':1',sep=''),NA,
          formatC(object$ER*100,digits=2,format="f"),
                     paste('[',formatC(object$ER.cl[1]*100,digits=2,format="f"),'-',formatC(object$ER.cl[2]*100,digits=2,format="f"),']'),
          formatC(object$accu*100,digits=2,format="f"),
                     paste('[',formatC(object$accu.cl[1]*100,digits=2,format="f"),'-',formatC(object$accu.cl[2]*100,digits=2,format="f"),']'),
          formatC(object$Youden,digits=4,format="f"),
                     paste('[',formatC(object$Youden.inf.cl,digits=4,format="f"),'-',formatC(object$Youden.sup.cl,digits=4,format="f"),']'),
          round(object$AUC,digits=4),NA  
        ),nrow = 14, ncol=2, byrow=TRUE,
        dimnames = list(c('Sample size:','Prevalence(%):','Sensitivity(%):','Specificity(%):',
                          'Postive predictive value(%):','Negative predictive value(%):',
                          'Positive likelihood ratio:','Negative likelihood ratio:',
                          'Diagnostic Odds Ratio:','Error trade off (FN : FP):','Error rate(%):',
                          'Accuracy(%):','Youden index:','Area under ROC curve:'),
                        c('Estimate', paste(object$Conf.limit,'Confidence limits')))  
  )
  diag.tab <- format(diag.tab,justify = "centre",na.encode = FALSE,trim = TRUE)
  colnames(diag.tab) <-format(colnames(diag.tab),justify = "centre",trim = TRUE)
  cat("--------------------------------------------------------------------------\n")
  print(diag.tab,quote=FALSE,na.print='')
  cat("--------------------------------------------------------------------------\n")
  invisible(diag.tab)
}
