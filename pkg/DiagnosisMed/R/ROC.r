ROC<-function(gold,
              test,
              CL=0.95,
              Cost=1,
              Prevalence=0,
              Plot=TRUE,
              Plot.point="Min.ROC.Dist",
              p.cex=1,
              Full=FALSE,
              Print=TRUE
              ){
  # A simple warning ...
  if(any(is.na(test) | is.na(gold))){
     stop('It seems there are NAs either in the index test or in the reference test. Consider imputing or removing NAs!')
  }
  test.table<-table(test,gold)
  if (dim(test.table)[2] != 2){
      stop("It seems that your gold standard has more than 2 categories")
  }
  CL<-CL
  cost<-Cost
  # Sample size
  sample.size<-sum(test.table)
  # Sample prevalence, replace by pop prevalence if adequate
  sample.prevalence<-(sum(test.table[,2])/sample.size)
  if (Prevalence==0){
    pop.prevalence<-sample.prevalence
  }
  if (Prevalence>0){
    (pop.prevalence<-Prevalence)
  }

  if (is.numeric(gold)==TRUE){
  X<-sort(test[gold==0]) 
  Y<-sort(test[gold==1]) 
  #X<-test[gold==0] 
  #Y<-test[gold==1] 
  AUC <- ((as.double(length(test[gold == 0]))) * (as.double(length(test[gold ==1]))) + ((as.double(length(test[gold == 0]))) * ((as.double(length(test[gold == 0]))) + 1))/2 - sum(rank(test,ties.method = "average")[gold == 0]))/((as.double(length(test[gold == 0]))) * (as.double(length(test[gold == 1]))))
  AUC[AUC < 0.5] <- 1 - AUC
  }

  if (is.factor(gold)==TRUE){
  #X<-test[gold=="negative"] 
  #Y<-test[gold=="positive"]
  X<-sort(test[gold=="negative"]) 
  Y<-sort(test[gold=="positive"]) 
  AUC <- ((as.double(length(test[gold == "negative"]))) * (as.double(length(test[gold == "positive"]))) + ((as.double(length(test[gold == "negative"]))) * ((as.double(length(test[gold == "negative"]))) + 1))/2 - sum(rank(test,ties.method = "average")[gold == "negative"]))/((as.double(length(test[gold == "negative"]))) * (as.double(length(test[gold == "positive"]))))
  AUC[AUC < 0.5] <- 1 - AUC
  }
  m<-as.double(length(X)) 
  n<-as.double(length(Y))
   
  test.summary<-round(c(summary(test),sd(test)),digits=5)
  test.summary<-rbind(test.summary,round(c(summary(X),sd(X)),digits=5))
  test.summary<-rbind(test.summary,round(c(summary(Y),sd(Y)),digits=5))
  colnames(test.summary)<-c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","SD")
  rownames(test.summary)<-c("Overall summary","Without disease", "With disease")
    
  #D10X<-function(Xi){(1/n)*sum(Y>=Xi)} 
  #D01Y<-function(Yi){(1/m)*sum(Yi>=X)} 
  D10X <- function(Xi) {(1/n) * sum(Y >= Xi[1])}
  D01Y <- function(Yi) {(1/m) * sum(Yi[1] >= X)}
  VAR.AUC<-sum((tapply(X,X,"D10X")-AUC)^2)/(m*(m-1))+sum((tapply(Y,Y,"D01Y")-AUC)^2)/(n*(n-1))
  SD.AUC<-sqrt(VAR.AUC)
  alpha<-1-CL
  AUC.summary<-c(AUC- qnorm(1-alpha/2)*SD.AUC,AUC,AUC+ qnorm(1-alpha/2)*SD.AUC)
  #names(AUC.summary)<-c("AUC inf conf limit", "AUC","AUC sup conf limit")
  
  #TP sum(test.table[i:nrow(test.table),2])
  #FP sum(test.table[i:nrow(test.table),1])
  #TN sum(test.table[1:i-1,1])
  #FN sum(test.table[1:i-1,2])
  D<-sum(test.table[,2])
  ND<-sum(test.table[,1])

  # Taking the rownames of the test.table to be results first column
  test.values<-(as.numeric(rownames(unclass(test.table))))
  test.diag.table<-as.data.frame(test.values)
  # Making a table with Se Sp PLR NLR PPV NPV and its confidence limits for each cut-off
  for (i in 1:nrow(test.diag.table)) {
    test.diag.table$TP[i] <- sum(test.table[i:nrow(test.table),2])
    test.diag.table$FN[i] <- sum(test.table[1:i-1,2])
    test.diag.table$FP[i] <- sum(test.table[i:nrow(test.table),1])
    test.diag.table$TN[i] <- sum(test.table[1:i-1,1])
  }  
  test.diag.table$Sensitivity <- round(test.diag.table$TP/D,digits=4)
  test.diag.table$Se.inf.cl <- round(binom.wilson(test.diag.table$TP,D,conf.level=CL)[4]$lower,digits=4)
  test.diag.table$Se.sup.cl <- round(binom.wilson(test.diag.table$TP,D,conf.level=CL)[5]$upper,digits=4)
  test.diag.table$Specificity <- round(test.diag.table$TN/ND,digits=4)
  test.diag.table$Sp.inf.cl <- round(binom.wilson(test.diag.table$TN,ND,conf.level=CL)[4]$lower,digits=4)
  test.diag.table$Sp.sup.cl <- round(binom.wilson(test.diag.table$TN,ND,conf.level=CL)[5]$upper,digits=4)
  test.diag.table$PPV <- round(test.diag.table$TP/(test.diag.table$TP + test.diag.table$FP),digits=4)
  test.diag.table$PPV.inf.cl <- round(binom.wilson(test.diag.table$TP,(test.diag.table$TP + test.diag.table$TP),conf.level=CL)[4]$lower,digits=4)
  test.diag.table$PPV.sup.cl <- round(binom.wilson(test.diag.table$TP,(test.diag.table$TP + test.diag.table$FN),conf.level=CL)[5]$upper,digits=4)
  test.diag.table$NPV <- round(test.diag.table$TN/(test.diag.table$TN + test.diag.table$FN),digits=4)
  test.diag.table$NPV.inf.cl <- round(binom.wilson(test.diag.table$TN,(test.diag.table$TN + test.diag.table$FN),conf.level=CL)[4]$lower,digits=4)
  test.diag.table$NPV.sup.cl <- round(binom.wilson(test.diag.table$TN,(test.diag.table$TN + test.diag.table$FN),conf.level=CL)[5]$upper,digits=4)
  test.diag.table$PLR<-round(test.diag.table$Sensitivity/(1-test.diag.table$Specificity),digits=2)
  test.diag.table$PLR.inf.cl<-round(exp(log(test.diag.table$PLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-test.diag.table$Sensitivity)/(
    (D)*test.diag.table$Specificity)+(test.diag.table$Specificity)/((ND)*(1-test.diag.table$Specificity)))),digits=2)
  test.diag.table$PLR.sup.cl<-round(exp(log(test.diag.table$PLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-test.diag.table$Sensitivity)/(
    (D)*test.diag.table$Specificity)+(test.diag.table$Specificity)/((ND)*(1-test.diag.table$Specificity)))),digits=2)
  test.diag.table$NLR<-round((1-test.diag.table$Sensitivity)/test.diag.table$Specificity,digits=2)
  test.diag.table$NLR.inf.cl<-round(exp(log(test.diag.table$NLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((test.diag.table$Sensitivity)/((D)*(1-test.diag.table$Sensitivity))+(1-test.diag.table$Specificity)/((ND)*(test.diag.table$Specificity)))),digits=2)
  test.diag.table$NLR.sup.cl<-round(exp(log(test.diag.table$NLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((test.diag.table$Sensitivity)/((D)*(1-test.diag.table$Sensitivity))+(1-test.diag.table$Specificity)/((ND)*(test.diag.table$Specificity)))),digits=2)
  test.diag.table$Accuracy <- (test.diag.table$TN + test.diag.table$TP)/sample.size
  test.diag.table$DOR <- ((test.diag.table$TN)*(test.diag.table$TP))/((test.diag.table$FP)*(test.diag.table$FN))
  test.diag.table$DOR<-ifelse(test.diag.table$DOR==Inf,NA,test.diag.table$DOR)  
  test.diag.table$Error.rate <- ((test.diag.table$FP)+(test.diag.table$FN))/sample.size
  test.diag.table$Accuracy.area <- ((test.diag.table$TP)*(test.diag.table$TN))/(D*ND)
  test.diag.table$Max.Se.Sp <- test.diag.table$Sensitivity + test.diag.table$Specificity
  test.diag.table$Youden <- test.diag.table$Sensitivity + test.diag.table$Specificity - 1
  test.diag.table$Se.equals.Sp <- abs(test.diag.table$Specificity-test.diag.table$Sensitivity)
  test.diag.table$MinRocDist <- (test.diag.table$Specificity-1)^2+(1-test.diag.table$Sensitivity)^2
  test.diag.table$Efficiency<-(test.diag.table$Sensitivity*(pop.prevalence))+((1-(pop.prevalence))*test.diag.table$Specificity)
  test.diag.table$MCT<-(1-(pop.prevalence))*(1-test.diag.table$Specificity)+(cost*(pop.prevalence))*(1-test.diag.table$Sensitivity)

  # Making a table with the test result for each best cut-off and attaching validity measures
  test.best.cutoff <- as.data.frame(rbind(
       test.diag.table[which.max(test.diag.table$Accuracy),c(1,6:11,18:20)],
       test.diag.table[which.max(test.diag.table$DOR),c(1,6:11,18:20)],
       test.diag.table[which.min(test.diag.table$Error.rate),c(1,6:11,18:20)],
       test.diag.table[which.max(test.diag.table$Accuracy.area),c(1,6:11,18:20)],
       test.diag.table[which.max(test.diag.table$Max.Se.Sp),c(1,6:11,18:20)],
       test.diag.table[which.max(test.diag.table$Youden),c(1,6:11,18:20)],
       test.diag.table[which.min(test.diag.table$Se.equals.Sp),c(1,6:11,18:20)],
       test.diag.table[which.min(test.diag.table$MinRocDist),c(1,6:11,18:20)],
       test.diag.table[which.max(test.diag.table$Efficiency),c(1,6:11,18:20)],
       test.diag.table[which.min(test.diag.table$MCT),c(1,6:11,18:20)]
       ))
  rownames(test.best.cutoff)<- c("Max. Accuracy", "Max. DOR","Min. Error rate",
   "Max. Accuracy area","Max. Sens+Spec","Max. Youden","Se=Sp","Min. ROC distance",
   "Max. Efficiency", "Min. MCT")
   
  #names(pop.prevalence)<-c("Informed disease prevalence - same as sample prevalence if not informed")
  #names(sample.prevalence)<-c("Observed prevalence by gold standard")
   reteval<-list(pop.prevalence=pop.prevalence,
                 sample.size=sample.size,
                 sample.prevalence=sample.prevalence,
                 test.summary=test.summary,
                 AUC.summary=AUC.summary,
                 test.table=test.table,
                 test.best.cutoff=test.best.cutoff,
                 test.diag.table=test.diag.table,
                 CL=CL,
                 cost=cost)
  
  class(reteval)<-"ROC"
  if(Print==TRUE){
     if(Full==TRUE){ print(reteval,Full=TRUE) }
     else{ print(reteval) }
  }
  # the plot commands
  if(Plot==TRUE){
  plot(reteval,Plot.point=Plot.point,p.cex=p.cex)
  }
  invisible(reteval)
}