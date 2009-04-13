ROC<-function(gold,
              test,
              CL=0.95,
              Cost=1,
              Prevalence=0,
              Plot=TRUE,
              Plot.point="Min.ROC.Dist",
              Print.full=FALSE,
              Print=TRUE
              ){
  # A simple warning ...
  test.table<-table(test,gold)
  if (dim(test.table)[2] != 2){
      stop("It seems that your gold standard has more than 2 categories")
  }
  CL<-CL
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
  # Making test summary, overall, disease only, without disease only
  test.summary<-round(c(summary(test),sd(test)),digits=5)
  test.summary<-rbind(test.summary,round(c(summary(test[gold==0]),sd(test[gold==0])),digits=5))
  test.summary<-rbind(test.summary,round(c(summary(test[gold==1]),sd(test[gold==1])),digits=5))
  colnames(test.summary)<-c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","SD")
  rownames(test.summary)<-c("Overall summary","Without disease", "With disease")
  # Estimating the AUC and confidence limits, inspired in auc{PresenceAbsence}
  # m length(test[gold == 0])
  # n length(test[gold == 1])
  AUC <- ((as.double(length(test[gold == 0]))) * (as.double(length(test[gold ==1]))) + ((as.double(length(test[gold == 0]))) * ((as.double(length(test[gold == 0]))) + 1))/2 - sum(rank(test,ties.method = "average")[gold == 0]))/((as.double(length(test[gold == 0]))) * (as.double(length(test[gold == 1]))))
  AUC[AUC < 0.5] <- 1 - AUC
  X<-test[gold==0] # pag 209
  Y<-test[gold==1] # pag 209
  }
  if (is.factor(gold)==TRUE){
  # The same tests summary but with different reference standard codes
  test.summary<-round(c(summary(test),sd(test)),digits=5)
  test.summary<-rbind(test.summary,round(c(summary(test[gold=="negative"]),sd(test[gold=="negative"])),digits=5))
  test.summary<-rbind(test.summary,round(c(summary(test[gold=="positive"]),sd(test[gold=="positive"])),digits=5))
  colnames(test.summary)<-c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","SD")
  rownames(test.summary)<-c("Overall summary","Without disease", "With disease")
  AUC <- ((as.double(length(test[gold =="negative"]))) * (as.double(length(test[gold =="positive"]))) + ((as.double(length(test[gold =="negative"]))) * ((as.double(length(test[gold =="negative"]))) + 1))/2 - sum(rank(test,ties.method = "average")[gold =="negative"]))/((as.double(length(test[gold =="negative"]))) * (as.double(length(test[gold == "positive"]))))
  AUC[AUC < 0.5] <- 1 - AUC
  X<-test[gold=="negative"] # pag 209
  Y<-test[gold=="positive"] # pag 209
  }
  m<-length(X) # pag 209  
  n<-length(Y) # pag 209  
  D10X<-function(Xi){(1/n)*sum(Y>=Xi)} # pag 211
  D01Y<-function(Yi){(1/m)*sum(Yi>=X)} # pag 211
  VAR.AUC<-sum((tapply(X,X,"D10X")-AUC)^2)/(m*(m-1))+sum((tapply(Y,Y,"D01Y")-AUC)^2)/(n*(n-1)) # pag 211
  SD.AUC<-sqrt(VAR.AUC)
  alpha<-1-CL
  AUC.summary<-c(AUC- qnorm(1-alpha/2)*SD.AUC,AUC,AUC+ qnorm(1-alpha/2)*SD.AUC)
  #names(AUC.summary)<-c("AUC inf conf limit", "AUC","AUC sup conf limit")
  
  #TP sum(test.table[i:nrow(test.table),2])
  #FP sum(test.table[i:nrow(test.table),1])
  #TN sum(test.table[1:i-1,1])
  #FN sum(test.table[1:i-1,2])

  # Taking the rownames of the test.table to be results first column
  test.values<-(as.numeric(rownames(unclass(test.table))))
  test.diag.table<-as.data.frame(test.values)
  # Making a table with Se Sp PLR NLR PPV NPV and its confidence limits for each cut-off
  for (i in 1:nrow(test.diag.table)) {
    test.diag.table$Sensitivity[i] <- round(sum(test.table[(i:nrow(test.table)),2])/
      sum(test.table[,2]),digits=4)
    test.diag.table$Se.inf.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[i:nrow(
      test.table),2]),sum(test.table[,2]),conf.level=CL)[4]),digits=4)
    test.diag.table$Se.sup.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[i:nrow(
      test.table),2]),sum(test.table[,2]),conf.level=CL)[5]),digits=4)
    test.diag.table$Specificity[i] <- round((sum(test.table[(1:i-1),1]))/(sum(
      test.table[,1])),digits=4)
    test.diag.table$Sp.inf.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[(
      1:i-1),1]),sum(test.table[,1]),conf.level=CL)[4]),digits=4)
    test.diag.table$Sp.sup.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[
      (1:i-1),1]),sum(test.table[,1]),conf.level=CL)[5]),digits=4)
    test.diag.table$PPV[i]<-round((sum(test.table[i:nrow(test.table),2]))/(
      (sum(test.table[i:nrow(test.table),2]))+
      (sum(test.table[i:nrow(test.table),1]))),digits=4)
    test.diag.table$PPV.inf.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[
      i:nrow(test.table),2]),((sum(test.table[i:nrow(test.table),2]))+(sum(
      test.table[i:nrow(test.table),1]))),conf.level=CL)[4]),digits=4)
    test.diag.table$PPV.sup.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[
      i:nrow(test.table),2]),((sum(test.table[i:nrow(test.table),2]))+(sum(
      test.table[i:nrow(test.table),1]))),conf.level=CL)[5]),digits=4)
    test.diag.table$NPV[i]<-round((sum(test.table[1:i-1,1]))/((sum(test.table[
      1:i-1,1]))+sum(test.table[1:i-1,2])),digits=4)
    test.diag.table$NPV.inf.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[
      1:i-1,1]),((sum(test.table[1:i-1,1]))+(sum(test.table[1:i-1,])))
      ,conf.level=CL)[4]),digits=4)
    test.diag.table$NPV.sup.cl[i]<-round(as.numeric(binom.wilson(sum(test.table[
      1:i-1,1]),((sum(test.table[1:i-1,1]))+(sum(test.table[1:i-1,])))
      ,conf.level=CL)[5]),digits=4)
  }
  test.diag.table$PLR<-round(test.diag.table$Sensitivity/(1-test.diag.table$Specificity),digits=2)
  test.diag.table$PLR.inf.cl<-round(exp(log(test.diag.table$PLR)-
    (qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-test.diag.table$Sensitivity)/
    ((sum(test.diag.table[,2]))*test.diag.table$Specificity)+
    (test.diag.table$Specificity)/((sum(test.diag.table[,1]))*
    (1-test.diag.table$Specificity)))),digits=2)
  test.diag.table$PLR.sup.cl<-round(exp(log(test.diag.table$PLR)+(qnorm(1-((1-CL)/
    2),mean=0,sd=1))*sqrt((1-test.diag.table$Sensitivity)/((sum(
    test.diag.table[,2]))*test.diag.table$Specificity)+(
    test.diag.table$Specificity)/((sum(test.diag.table[,1]))*(1-
    test.diag.table$Specificity)))),digits=2)
  test.diag.table$NLR<-round((1-test.diag.table$Sensitivity)/test.diag.table$Specificity,digits=2)
  test.diag.table$NLR.inf.cl<-round(exp(log(test.diag.table$NLR)-(qnorm(1-((1-CL)/2)
    ,mean=0,sd=1))*sqrt((test.diag.table$Sensitivity)/((sum(test.diag.table[,
    2]))*(1-test.diag.table$Sensitivity))+(1-test.diag.table$Specificity)/((
    sum(test.diag.table[,1]))*(test.diag.table$Specificity)))),digits=2)
  test.diag.table$NLR.sup.cl<-round(exp(log(test.diag.table$NLR)+(qnorm(1-((1-CL)/2)
    ,mean=0,sd=1))*sqrt((test.diag.table$Sensitivity)/((sum(test.diag.table[,2])
    )*(1-test.diag.table$Sensitivity))+(1-test.diag.table$Specificity)/((sum(
    test.diag.table[,1]))*(test.diag.table$Specificity)))),digits=2)

  # The following is not result but will help to find the best cut-offs
  test.cutoff.table<-as.data.frame(test.values)
  for(i in 1:nrow(test.cutoff.table)) {
    #Accuracy is (TP+TN)/sample.size
    test.cutoff.table$Accuracy[i]<-round((sum(test.table[i:nrow(test.table),2])+
      sum(test.table[1:i-1,1]))/sample.size,digits=4)
    test.cutoff.table$DOR[i]<-round(((sum(test.table[i:nrow(test.table),2])*(sum(test.table[1:i-1,1])))/
    ((sum(test.table[i:nrow(test.table),1]))*(sum(test.table[1:i-1,2])))),digits=2)
    test.cutoff.table$Error.rate[i]<-round(((sum(test.table[1:i-1,2]))+(sum(
      test.table[i:nrow(test.table),1])))/sample.size,digits=4)
    #Error.trade.off - error out of limits, only zeros can be used with negative subscripts
    #Error.trade.off[i]<-(sum(test.table[1:i-1,2])-sum(test.table[1:i-2,2]))
        #    (sum(test.table[i:nrow(test.table),1])-sum(test.table[i+
        #    1:nrow(test.table),1]))
    test.cutoff.table$Accuracy.area[i]<-round((sum(test.table[i:nrow(test.table)
      ,2])*sum(test.table[1:i-1,1]))/((sample.size-pop.prevalence)*pop.prevalence),digits=4)
  }
  test.cutoff.table$DOR<-ifelse(test.cutoff.table$DOR==Inf,NA,test.cutoff.table$DOR)
  test.cutoff.table$Max.Se.Sp<-round((test.diag.table$Sensitivity+
    test.diag.table$Specificity),digits=4)
  test.cutoff.table$Youden<-test.diag.table$Sensitivity+test.diag.table$Specificity-1
  test.cutoff.table$Se.equals.Sp<-abs(test.diag.table$Specificity-test.diag.table$Sensitivity)
  test.cutoff.table$MinRocDist<-(test.diag.table$Specificity-1)^2+(1-test.diag.table$Sensitivity)^2
  # Efficiency= Se*prevalence+(1-prevalence)*Sp
  test.cutoff.table$Efficiency<-(test.diag.table$Sensitivity*(pop.prevalence)
    )+((1-(pop.prevalence))*test.diag.table$Specificity)
  # MissclassificatioCostTerm(MCT)=(1-prevalence)(1-Sp)+r*prevalence(1-Se) - r=cost(FN)/cost(FP)
  test.cutoff.table$MCT<-(1-(pop.prevalence))*(1-test.diag.table$Specificity)+(Cost*
    (pop.prevalence))*(1-test.diag.table$Sensitivity)
  # Making a table with the test result for each best cut-off and attaching validity measures
  best.cutoff<-test.cutoff.table$test.values[which.max(test.cutoff.table$Accuracy)]
  test.best.cutoff<-cbind(best.cutoff,test.diag.table[which.max(test.cutoff.table$Accuracy)
           ,2:7],test.diag.table[which.max(test.cutoff.table$Accuracy),14:16])
  best.cutoff<-test.cutoff.table$test.values[which.max(test.cutoff.table$DOR)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.max(test.cutoff.table$DOR),2:7],test.diag.table[which.max
           (test.cutoff.table$DOR),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.min(test.cutoff.table$Error.rate)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.min(test.cutoff.table$Error.rate),2:7],test.diag.table[which.min
           (test.cutoff.table$Error.rate),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.max(test.cutoff.table$Accuracy.area)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.max(test.cutoff.table$Accuracy.area),2:7],test.diag.table[which.max
           (test.cutoff.table$Accuracy.area),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.max(test.cutoff.table$Max.Se.Sp)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.max(test.cutoff.table$Max.Se.Sp),2:7],test.diag.table[which.max
           (test.cutoff.table$Max.Se.Sp),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.max(test.cutoff.table$Youden)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.max(test.cutoff.table$Youden),2:7],test.diag.table[which.max
           (test.cutoff.table$Youden),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.min(test.cutoff.table$Se.equals.Sp)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.min(test.cutoff.table$Se.equals.Sp),2:7],test.diag.table[which.min
           (test.cutoff.table$Se.equals.Sp),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.min(test.cutoff.table$MinRocDist)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.min(test.cutoff.table$MinRocDist),2:7],test.diag.table[which.min
           (test.cutoff.table$MinRocDist),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.max(test.cutoff.table$Efficiency)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.max(test.cutoff.table$Efficiency),2:7],test.diag.table[which.max
           (test.cutoff.table$Efficiency),14:16]))
  best.cutoff<-test.cutoff.table$test.values[which.min(test.cutoff.table$MCT)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,test.diag.table[
           which.min(test.cutoff.table$MCT),2:7],test.diag.table[which.min
           (test.cutoff.table$MCT),14:16]))
  rownames(test.best.cutoff)<- c("Max. Accuracy", "Max. DOR","Min. Error rate",
   "Max. Accuracy area","Max. Sens+Spec","Max. Youden","Se=Sp","Min. ROC distance",
   "Max. Efficiency", "Min. MCT")
  rm(best.cutoff)
  # the plot commands
  if(Plot==TRUE){
    plot(1-test.diag.table$Specificity,test.diag.table$Sensitivity,type="l",
      col=1,xlab="1-Specificity",ylab="Sensitivity",xlim=c(0,1),ylim=c(0,1))
    grid()
    segments(0,0,1,1,col="lightgray")
    if(Plot.point=="Max.Accuracy")
      {points(1-test.diag.table$Specificity[which.max(test.cutoff.table$Accuracy)],
         test.diag.table$Sensitivity[which.max(test.cutoff.table$Accuracy)],
         col=1,pch=19)
      title(sub="Cut-off estimated by maximazing accuracy",cex.sub=0.85)   
      legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.max
            (test.cutoff.table$Accuracy)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.max
            (test.cutoff.table$Accuracy)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.max
            (test.cutoff.table$Accuracy)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Max.DOR")
      {points(1-test.diag.table$Specificity[which.max(test.cutoff.table$DOR)],
         test.diag.table$Sensitivity[which.max(test.cutoff.table$DOR)],
         col=1,pch=19)
       title(sub="Cut-off estimated by maximazing diagnostic odds ratio",cex.sub=0.85)  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.max
            (test.cutoff.table$DOR)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.max
            (test.cutoff.table$DOR)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.max
            (test.cutoff.table$DOR)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Error.rate")
      {points(1-test.diag.table$Specificity[which.min(test.cutoff.table$Error.rate)],
         test.diag.table$Sensitivity[which.min(test.cutoff.table$Error.rate)],
         col=1,pch=19)
       title(sub="Cut-off estimated by minimizing error rate",cex.sub=0.85)  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.min
            (test.cutoff.table$Error.rate)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.min
            (test.cutoff.table$Error.rate)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.min
            (test.cutoff.table$Error.rate)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Max.Accuracy.area")
      {points(1-test.diag.table$Specificity[which.max(test.cutoff.table$Accuracy.area)],
         test.diag.table$Sensitivity[which.max(test.cutoff.table$Accuracy.area)],
         col=1,pch=19)
       title(sub="Cut-off estimated by maximazing the area related to accuracy",cex.sub=0.85)  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.max
            (test.cutoff.table$Accuracy.area)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.max
            (test.cutoff.table$Accuracy.area)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.max
            (test.cutoff.table$Accuracy.area)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Max.Sens+Spec")
      {points(1-test.diag.table$Specificity[which.max(test.cutoff.table$Max.Se.Sp)],
         test.diag.table$Sensitivity[which.max(test.cutoff.table$Max.Se.Sp)],
         col=1,pch=19)
      title(sub="Cut-off value where the sum Se + Sp is maximized",cex.sub=0.85)         
      legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.max
            (test.cutoff.table$Max.Se.Sp)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.max
            (test.cutoff.table$Max.Se.Sp)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.max
            (test.cutoff.table$Max.Se.Sp)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Max.Youden")
      {points(1-test.diag.table$Specificity[which.max(test.cutoff.table$Youden)],
         test.diag.table$Sensitivity[which.max(test.cutoff.table$Youden)],
         col=1,pch=19)
       title(sub="Cut-off estimated by maximazing Youden Index",cex.sub=0.85)         
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.max
            (test.cutoff.table$Youden)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.max
            (test.cutoff.table$Youden)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.max
            (test.cutoff.table$Youden)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Se=Sp")
      {points(1-test.diag.table$Specificity[which.min(test.cutoff.table$Se.equals.Sp)],
         test.diag.table$Sensitivity[which.min(test.cutoff.table$Se.equals.Sp)],
         col=1,pch=19)
       title(sub="Cut-off value where Se is the closest to Sp",cex.sub=0.85)
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.min
            (test.cutoff.table$Se.equals.SP)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.min
            (test.cutoff.table$Se.equals.SP)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.min
            (test.cutoff.table$Se.equals.SP)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Min.ROC.Dist")
      {points(1-test.diag.table$Specificity[which.min(test.cutoff.table$MinRocDist)],
         test.diag.table$Sensitivity[which.min(test.cutoff.table$MinRocDist)],
         col=1,pch=19)
       title(sub="Cut-off that minimizes the distance between the curve and upper left corner",cex.sub=0.85)         
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.min
            (test.cutoff.table$MinRocDist)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.min
            (test.cutoff.table$MinRocDist)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.min
            (test.cutoff.table$MinRocDist)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Max.Efficiency")
      {points(1-test.diag.table$Specificity[which.max(test.cutoff.table$Efficiency)],
         test.diag.table$Sensitivity[which.max(test.cutoff.table$Efficiency)],
         col=1,pch=19)
       title(sub=paste("Cut-off maximizing efficiency: population prevalence =",
             formatC(pop.prevalence,digits=2)),cex.sub=0.85)         
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.max
            (test.cutoff.table$Efficiency)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.max
            (test.cutoff.table$Efficiency)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.max
            (test.cutoff.table$Efficiency)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
    if(Plot.point=="Min.MCT")
      {points(1-test.diag.table$Specificity[which.min(test.cutoff.table$MCT)],
         test.diag.table$Sensitivity[which.min(test.cutoff.table$MCT)],
         col=1,pch=19)
       title(sub=paste("Cut-off minimazing MCT: population prevalence =",
             formatC(pop.prevalence,digits=2),"; cost(FN)/cost(FP)=",
             formatC(Cost,digits=2)),cex.sub=0.85)                  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(test.cutoff.table$test.values[which.min
            (test.cutoff.table$MCT)],digits=4)),
         paste("Sensitivity:",formatC(test.diag.table$Sensitivity[which.min
            (test.cutoff.table$MCT)],digits=4)),
         paste("Specificity:",formatC(test.diag.table$Specificity[which.min
            (test.cutoff.table$MCT)],digits=4)),
         paste("AUC:",formatC(AUC,digits=4))
       )),bty="n")}
  }
  #names(pop.prevalence)<-c("Informed disease prevalence - same as sample prevalence if not informed")
  #names(sample.prevalence)<-c("Observed prevalence by gold standard")
   reteval<-list(pop.prevalence=pop.prevalence,sample.size=sample.size,sample.prevalence=sample.prevalence,test.summary=test.summary,AUC.summary=AUC.summary,test.best.cutoff=test.best.cutoff,test.diag.table=test.diag.table,CL=CL,test.cutoff.table=test.cutoff.table)
  invisible(reteval)
  class(reteval)<-"ROC"
  if(Print==TRUE){
     if(Print.full==TRUE){ print(reteval,Full=TRUE) }
     else{ print(reteval) }
  }
}