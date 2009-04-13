diagnosis <- function(gold,test,CL=0.95,print=TRUE,plot=FALSE){
  #require(epitools)
  #require(epicalc)
  # to do ...
  # by option (multicenter)
  # test with 3 categories (indeterminate results)
  #   testef <- as.factor(teste)
  #   if(nlevels(teste)  2)
  #                  {
                     # analysis this way
  #                   }
  tab<-table(test,gold,dnn=c(deparse(substitute(test)),deparse(substitute(gold))))
  dimnames(tab) <- list(test = c("negative" , "positive"), gold.standard = c("negative","positive"))
  tabmarg<-addmargins(tab)
  Conf.limit<-CL
  #dnn is option of table command - specifies the names of row and column
  TN<-tab[1,1]
  FN<-tab[1,2]
  FP<-tab[2,1]
  TP<-tab[2,2]
  # sample size
  n<-sum(tab)
  # prevalence
  p<-(TP+FN)/n
  # sensitivity and confidence limits
  Se<-TP/(TP+FN)
  Se.cl<-as.numeric(binom.wilson(TP, TP+FN, conf.level = CL)[4:5])
  # especificity and confidence limits
  Sp<-TN/(FP+TN)
  Sp.cl<-as.numeric(binom.wilson(TN, FP+TN, conf.level = CL)[4:5])
  # positive and negative likelyhood ratios and confidence limits
  PLR<-Se/(1-Sp)
  # LR confidence limists inspired in epi.tests{epiR}
  PLR.inf.cl<-exp(log(PLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-Se)/(
    (TP+FN)*Sp)+(Sp)/((FP+TN)*(1-Sp))))
  PLR.sup.cl<-exp(log(PLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-Se)/(
    (TP+FN)*Sp)+(Sp)/((FP+TN)*(1-Sp))))
  NLR<-(1-Se)/Sp
  NLR.inf.cl<-exp(log(NLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((Se)/((TP+
    FN)*(1-Se))+(1-Sp)/((FP+TN)*(Sp))))
  NLR.sup.cl<-exp(log(NLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((Se)/((TP+
    FN)*(1-Se))+(1-Sp)/((FP+TN)*(Sp))))
  #accuracy and confidence limits
  accu<-(TP+TN)/n
  accu.cl<-as.numeric(binom.wilson(TP+TN, n, conf.level = CL)[4:5])
  # positive and negative predictive values and confidence limits
  PPV<-TP/(TP+FP)
  PPV.cl<-as.numeric(binom.wilson(TP, TP+FP, conf.level = CL)[4:5])
  NPV<-TN/(TN+FN)
  NPV.cl<-as.numeric(binom.wilson(TN, TN+FN, conf.level = CL)[4:5])
  # diagnostic odds ratio and confidence limits
  OR<-oddsratio(tab,conf.level = CL)
  DOR<-OR$measure[2,1]
  #DOR<-(TP*TN)/(FP*FN)
  DOR.inf.cl<-OR$measure[2,2]
  DOR.sup.cl<-OR$measure[2,3]
  rm(OR)
  # error rate and error trade
  #ER<-((FN/(FN+TN))*p)+(((FP/(FP+TP))*(TN+FP))
  ER<-(FN+FP)/n
  ER.cl<-as.numeric(binom.wilson(FN+FP, n, conf.level = CL)[4:5])
  ET<-(FN/FP)
  # pre-test and pos-test odds (to do)
  # area under ROC curve
  if(plot==TRUE)
    {ROC<-roc.from.table(tab, graph = TRUE)}
  if(plot==FALSE)
    {ROC<-roc.from.table(tab, graph = FALSE)}
  AUC<-ROC$auc
  # gives same results as AUC<-(Se+Sp)/2
  Youden<-Se+Sp-1
  Youden.inf.cl<-Youden-qnorm(CL/2)*sqrt(((Se * (1 - Se))/(TP+FN) +
           ((Sp * (1 - Sp))/(TN+FP))))
  Youden.sup.cl<-Youden+qnorm(CL/2)*sqrt(((Se * (1 - Se))/(TP+FN) +
           ((Sp * (1 - Sp))/(FP+TN))))
  rm(ROC)
  rm(tab)
  # results evaluations
  reteval <- list(tabmarg=tabmarg,n=n,p=p,Se=Se,Se.cl=Se.cl,Sp=Sp,Sp.cl=Sp.cl,PLR=PLR,
    PLR.inf.cl=PLR.inf.cl,PLR.sup.cl=PLR.sup.cl,NLR=NLR,NLR.inf.cl=NLR.inf.cl,
    NLR.sup.cl=NLR.sup.cl,accu=accu,accu.cl=accu.cl,PPV=PPV,PPV.cl=PPV.cl,NPV=NPV,NPV.cl=NPV.cl,
    DOR=DOR,DOR.inf.cl=DOR.inf.cl,DOR.sup.cl=DOR.sup.cl,ET=ET,ER=ER,ER.cl=ER.cl,
    Youden=Youden,Youden.inf.cl=Youden.inf.cl,Youden.sup.cl=Youden.sup.cl,AUC=AUC,
    Conf.limit=Conf.limit)
  class(reteval) <- "diag"
  if(print==TRUE)  {print(reteval)}
  invisible(reteval)
}