#' @export
#' @import epitools
#' @import AMORE

TGROC<-function(gold,
                test,
                Cost=1,
                CL=0.95,
                Inconclusive=0.95,
                Prevalence=0,
                t.max=NULL,
                t.min=NULL,
                precision=.0001,  
                n.neurons=c(1,5,1),
                learning.rate.global=1e-2,
                momentum.global=0.3,
                error.criterium="LMS",
                Stao=NA,
                hidden.layer="sigmoid",
                output.layer="sigmoid",
                method="ADAPTgdwm",
                report=FALSE,
                show.step=5000,
                n.shows=1,
                Plot="Both",
                Plot.inc.range=TRUE,
                Plot.Cl=FALSE,
                Plot.cutoff="None",
                cex=0.5,
                cex.sub=0.85,
                Print=TRUE){
                
  #require(epitools)
  #TP sum(test.table[i:nrow(test.table),2])
  #FP sum(test.table[i:nrow(test.table),1])
  #TN sum(test.table[1:i-1,1])
  #FN sum(test.table[1:i-1,2])
  
  test.table<-table(test,gold)
  if (dim(test.table)[2] != 2){
      stop("It seems that your gold standard has more than 2 categories!")
  }
  if(is.null(precision)||!is.numeric(precision)){
      stop("Precision must be set to a numeric value!")
  }
  sample.prevalence<-(sum(test.table[,2]))/(sum(test.table))
  if (Prevalence==0){
    pop.prevalence<-sample.prevalence
  }
  if (Prevalence>0){
    (pop.prevalence<-Prevalence)
  }

  names(sample.prevalence)<-c("Disease prevalence in the sample")
  names(pop.prevalence)<-c("Informed disease prevalence in the population")
  sample.size<-sum(test.table)
  names(sample.size)<-c("Sample size")
  test.summary<-round(c(summary(test),sd(test)),digits=5)
  names(test.summary)<-c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","SD")
  cost<-Cost
  names(cost)<-c("Informed costs(FN)/costs(FP)")
  conf.limit<-CL
  inc<-Inconclusive
  names(inc)<-"Inconclusive tolerance level"

  D<-sum(test.table[,2])
  ND<-sum(test.table[,1])

  # Taking the rownames of the test.table to be results first column
  test.values<-(as.numeric(rownames(unclass(test.table))))
  non.parametric<-as.data.frame(test.values)
  # Making a table with Se Sp PLR NLR PPV NPV and its confidence limits for each cut-off
  for (i in 1:nrow(non.parametric)) {
    non.parametric$TP[i] <- sum(test.table[i:nrow(test.table),2])
    non.parametric$FN[i] <- sum(test.table[1:i-1,2])
    non.parametric$FP[i] <- sum(test.table[i:nrow(test.table),1])
    non.parametric$TN[i] <- sum(test.table[1:i-1,1])
  }  

  non.parametric$Sensitivity <- round(non.parametric$TP/D,digits=4)
  non.parametric$Se.inf.cl <- round(binom.wilson(non.parametric$TP,D,conf.level=CL)[4]$lower,digits=4)
  non.parametric$Se.sup.cl <- round(binom.wilson(non.parametric$TP,D,conf.level=CL)[5]$upper,digits=4)
  non.parametric$Specificity <- round(non.parametric$TN/ND,digits=4)
  non.parametric$Sp.inf.cl <- round(binom.wilson(non.parametric$TN,ND,conf.level=CL)[4]$lower,digits=4)
  non.parametric$Sp.sup.cl <- round(binom.wilson(non.parametric$TN,ND,conf.level=CL)[5]$upper,digits=4)

  non.parametric$PLR<-round(non.parametric$Sensitivity/(1-non.parametric$Specificity),digits=2)
  non.parametric$PLR.inf.cl<-round(exp(log(non.parametric$PLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-non.parametric$Sensitivity)/(
    (D)*non.parametric$Specificity)+(non.parametric$Specificity)/((ND)*(1-non.parametric$Specificity)))),digits=2)
  non.parametric$PLR.sup.cl<-round(exp(log(non.parametric$PLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-non.parametric$Sensitivity)/(
    (D)*non.parametric$Specificity)+(non.parametric$Specificity)/((ND)*(1-non.parametric$Specificity)))),digits=2)
    
  # Se=Sp cut-off
  non.parametric$Se.equals.Sp<-abs(non.parametric$Specificity-non.parametric$Sensitivity)
  # Efficiency= Se*prevalence+(1-prevalence)*Se
  non.parametric$Efficiency<-(non.parametric$Sensitivity*(pop.prevalence))+((1-(pop.prevalence))*non.parametric$Specificity)
  # MissClassificatioCost(MCT)=(1-prevalence)(1-Sp)+r*prevalence(1-Se) - r=cost(FN)/cost(FP)
  non.parametric$MCT<-(1-(pop.prevalence))*(1-non.parametric$Specificity)+(Cost*(pop.prevalence))*(1-non.parametric$Sensitivity)

#  np.test.best.cutoff<-subset(non.parametric,subset=c(
#       non.parametric[which.min(non.parametric$Se.equals.Sp)],
#       non.parametric[which.max(non.parametric$Efficiency)],
#       non.parametric[which.min(non.parametric$MCT)],
#       select=test.values:PLR.sup.cl       
#       )) Does Not work... somethig worng with the lines or columns selection

  np.test.best.cutoff<-as.data.frame(rbind(
       non.parametric[which.min(non.parametric$Se.equals.Sp),c(1,6:14)],
       non.parametric[which.max(non.parametric$Efficiency),c(1,6:14)],
       non.parametric[which.min(non.parametric$MCT),c(1,6:14)]))
  
  #best.cutoff,non.parametric[which.min(non.parametric$Se.equals.Sp),1:10])
  #best.cutoff<-non.parametric$test.values[which.max(non.parametric$Efficiency)]
  #np.test.best.cutoff<-rbind(np.test.best.cutoff,cbind(best.cutoff,non.parametric[
  #         which.max(non.parametric$Efficiency),2:10]))
  #best.cutoff<-non.parametric$test.values[which.min(non.parametric$MCT)]
  #np.test.best.cutoff<-rbind(np.test.best.cutoff,cbind(best.cutoff,non.parametric[
  #         which.min(non.parametric$MCT),2:10]))
  rownames(np.test.best.cutoff)<- c("Se=Sp","Max. Efficiency","Min. MCT")

  non.parametric.inconclusive<-as.data.frame(rbind(
      non.parametric[which.min(abs(Inconclusive-non.parametric$Sensitivity)),c(1,6:14)],
      non.parametric[which.min(abs(Inconclusive-non.parametric$Specificity)),c(1,6:14)]))
  rownames(non.parametric.inconclusive)<-c("Lower inconclusive","Upper inconclusive")

  if(is.null(t.max)){
  t.max<-max(non.parametric$test.values)
  }
  if(is.null(t.min)){  
  t.min<-min(non.parametric$test.values)
  }
  net <- newff(n.neurons=n.neurons, learning.rate.global=learning.rate.global, momentum.global=momentum.global,
        error.criterium=error.criterium, Stao=Stao, hidden.layer=hidden.layer, 
        output.layer=output.layer, method=method)
  net.Se <- train(net,P=non.parametric$test.values,T=non.parametric$Sensitivity,error.criterium=error.criterium,
        report=report,show.step=show.step,n.shows=n.shows)
  net.Sp <- train(net,P=non.parametric$test.values,T=non.parametric$Specificity,error.criterium=error.criterium,
        report=report,show.step=show.step,n.shows=n.shows)
  test.values<-seq(t.min,t.max,precision)
  parametric<-as.data.frame(test.values)
  parametric$Sensitivity <- as.numeric(sim(net.Se$net, test.values))
  parametric$Se.inf.cl<-parametric$Sensitivity - qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Sensitivity * 
        (1 - parametric$Sensitivity))/(sample.size * sample.prevalence)))
  parametric$Se.sup.cl<-parametric$Sensitivity + qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Sensitivity * 
        (1 - parametric$Sensitivity))/(sample.size * sample.prevalence)))
  parametric$Specificity <- as.numeric(sim(net.Sp$net, test.values))
  parametric$Sp.inf.cl <- parametric$Specificity - qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Specificity * 
        (1 - parametric$Specificity))/(sample.size * (1-sample.prevalence))))
  parametric$Sp.sup.cl <- parametric$Specificity + qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Specificity * 
        (1 - parametric$Specificity))/(sample.size * (1-sample.prevalence))))
  parametric$PLR <- parametric$Sensitivity/(1-parametric$Specificity) 
  parametric$PLR.inf.cl<-exp(log(parametric$PLR)-(qnorm(1-((1-conf.limit)/2),mean=0,sd=1))*sqrt((1-parametric$Sensitivity)/
        ((sample.size * sample.prevalence)*parametric$Specificity)+(parametric$Specificity)/((sample.size * 
        (1-sample.prevalence))*(1-parametric$Specificity))))
  parametric$PLR.sup.cl<-exp(log(parametric$PLR)+(qnorm(1-((1-conf.limit)/2),mean=0,sd=1))*sqrt((1-parametric$Sensitivity)/
        ((sample.size * sample.prevalence)*parametric$Specificity)+(parametric$Specificity)/((sample.size * 
        (1-sample.prevalence))*(1-parametric$Specificity))))            
  #parametric$NLR <- (1-parametric$Specificity)/parametric$Sensitivity
  #parametric$NLR.inf.cl <- exp(log(parametric$NLR)-(qnorm(1-((1-(1-conf.limit))/2),mean=0,sd=1))*
  #     sqrt((parametric$Sensitivity)/((sample.size * sample.prevalence)*(1-parametric$Sensitivity))+(1-parametric$Specificity)/
  #     ((sample.size * (1-sample.prevalence))*(parametric$Specificity))))
  #parametric$NLR.sup.cl <- exp(log(parametric$NLR)+(qnorm(1-((1-conf.limit)/2),mean=0,sd=1))*sqrt((parametric$Sensitivity)/
  #     ((sample.size * sample.prevalence)*(1-parametric$Sensitivity))+(1-parametric$Specificity)/((sample.size * 
  #     (1-sample.prevalence))*(parametric$Specificity))))
   parametric$Se.equals.Sp<-abs(parametric$Sensitivity-parametric$Specificity)
   parametric$Efficiency<-parametric$Sensitivity*pop.prevalence+(1-pop.prevalence)*parametric$Specificity
   parametric$MCT<-(1-(pop.prevalence))*(1-parametric$Specificity)+(cost*(pop.prevalence))*(1-parametric$Sensitivity)
   
   parametric.inconclusive<-as.data.frame(rbind(
             parametric[which.min(abs(inc-parametric$Sensitivity)),1:10],
             parametric[which.min(abs(inc-parametric$Specificity)),1:10]
             ))
   rownames(parametric.inconclusive)<-c("Lower inconclusive","Upper inconclusive")

   par.test.best.cutoff<-as.data.frame(rbind(
             parametric[which.min(parametric$Se.equals.Sp),1:10],
             parametric[which.max(parametric$Efficiency),1:10],
             parametric[which.min(parametric$MCT),1:10]
             ))
   rownames(par.test.best.cutoff)<- c("Se=Sp","Max. Efficiency","Min. MCT")  

#  rm(test.values)
  if(non.parametric.inconclusive[1,1]>non.parametric.inconclusive[2,1]){
     warning("Non-parametric lower inconclusive limit is higher than upper inconclusive limit.")
  }
  if(parametric.inconclusive[1,1]>parametric.inconclusive[2,1]){
     warning("Parametric lower inconclusive limit is higher than upper inconclusive limit.")
  }
  if(np.test.best.cutoff[1,1]>non.parametric.inconclusive[2,1]|
     np.test.best.cutoff[2,1]>non.parametric.inconclusive[2,1]|
     np.test.best.cutoff[3,1]>non.parametric.inconclusive[2,1]){
     warning("At least one of the non-parametric best cut-off is higher then upper inconclusive limit.")
  }        
  if(np.test.best.cutoff[1,1]<non.parametric.inconclusive[1,1]|
     np.test.best.cutoff[2,1]<non.parametric.inconclusive[1,1]|
     np.test.best.cutoff[3,1]<non.parametric.inconclusive[1,1]){
     warning("At least one of the non-parametric best cut-off is lower then lower inconclusive limit.")
  }
  if(par.test.best.cutoff[1,1]>parametric.inconclusive[2,1]|
     par.test.best.cutoff[2,1]>parametric.inconclusive[2,1]|
     par.test.best.cutoff[3,1]>parametric.inconclusive[2,1]){
     warning("At least one of the parametric best cut-off is higher then upper inconclusive limit.")
  }        
  if(par.test.best.cutoff[1,1]<parametric.inconclusive[1,1]|
     par.test.best.cutoff[2,1]<parametric.inconclusive[1,1]|
     par.test.best.cutoff[3,1]<parametric.inconclusive[1,1]){
     warning("At least one of the parametric best cut-off is lower then lower inconclusive limit.")
  }                  
  reteval<-list(sample.size = sample.size,
                sample.prevalence = sample.prevalence,
                pop.prevalence = pop.prevalence,
                cost = cost,
                test.summary = test.summary,
                inc = inc,
                conf.limit = conf.limit,
                non.parametric = non.parametric,
                non.parametric.inconclusive = non.parametric.inconclusive,
                np.test.best.cutoff = np.test.best.cutoff,
                parametric.inconclusive=parametric.inconclusive,
                par.test.best.cutoff=par.test.best.cutoff,
                parametric=parametric)
  class(reteval)<-"TGROC"
  if(Print==TRUE){
     print(reteval)
  }
  if(Plot!="None"){
  plot(reteval,
       Plot=Plot,
       Plot.inc.range=Plot.inc.range,
       Plot.Cl=Plot.Cl,
       Plot.cutoff=Plot.cutoff,
       cex=cex,
       cex.sub=cex.sub)
   }    
  invisible(reteval)
}