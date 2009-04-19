TGROC<-function(gold,
                test,
                Cost=1,
                CL=0.95,
                Inconclusive=0.95,
                Prevalence=0,
                Plot=TRUE,
                Plot.inc.range=TRUE,
                Plot.Cl=FALSE,
                Plot.cutoff="None",
                cex.sub=0.85,
                Print=TRUE){
                
  #require(epitools)
  #TP sum(test.table[i:nrow(test.table),2])
  #FP sum(test.table[i:nrow(test.table),1])
  #TN sum(test.table[1:i-1,1])
  #FN sum(test.table[1:i-1,2])
  
  test.table<-table(test,gold)
  if (dim(test.table)[2] != 2){
      stop("It seems that your gold standard has more than 2 categories")
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
  test.values<-(as.numeric(rownames(unclass(test.table))))
  non.parametric<-as.data.frame(test.values)
  for (i in 1:nrow(non.parametric)) {
    non.parametric$Sensitivity[i] <- round(sum(test.table[(i:nrow(test.table)),2])/
      sum(test.table[,2]),digits=4)
    non.parametric$Se.inf.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[i:nrow(
      test.table),2]),sum(test.table[,2]),conf.level=CL)[4]),digits=4)
    non.parametric$Se.sup.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[i:nrow(
      test.table),2]),sum(test.table[,2]),conf.level=CL)[5]),digits=4)
    non.parametric$Specificity[i] <- round((sum(test.table[(1:i-1),1]))/(sum(
      test.table[,1])),digits=4)
    non.parametric$Sp.inf.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[(
      1:i-1),1]),sum(test.table[,1]),conf.level=CL)[4]),digits=4)
    non.parametric$Sp.sup.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[
      (1:i-1),1]),sum(test.table[,1]),conf.level=CL)[5]),digits=4)
    }
  non.parametric$PLR<-round(non.parametric$Sensitivity/(1-non.parametric$Specificity),digits=2)
  non.parametric$PLR.inf.cl<-round(exp(log(non.parametric$PLR)-
    (qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-non.parametric$Sensitivity)/
    ((sum(non.parametric[,2]))*non.parametric$Specificity)+
    (non.parametric$Specificity)/((sum(non.parametric[,1]))*
    (1-non.parametric$Specificity)))),digits=2)
  non.parametric$PLR.sup.cl<-round(exp(log(non.parametric$PLR)+(qnorm(1-((1-CL)/
    2),mean=0,sd=1))*sqrt((1-non.parametric$Sensitivity)/((sum(
    non.parametric[,2]))*non.parametric$Specificity)+(
    non.parametric$Specificity)/((sum(non.parametric[,1]))*(1-
    non.parametric$Specificity)))),digits=2)

  cutoff.table<-as.data.frame(test.values)
  cutoff.table$Se.equals.Sp<-abs(non.parametric$Specificity-non.parametric$Sensitivity)
  # Efficiency= Se*prevalence+(1-prevalence)*Se
  cutoff.table$Efficiency<-(non.parametric$Sensitivity*(pop.prevalence)
    )+((1-(pop.prevalence))*non.parametric$Specificity)
  # MissClassificatioCost(MCT)=(1-prevalence)(1-Sp)+r*prevalence(1-Se) - r=cost(FN)/cost(FP)
  cutoff.table$MCT<-(1-(pop.prevalence))*(1-non.parametric$Specificity)+(Cost*
    (pop.prevalence))*(1-non.parametric$Sensitivity)


  best.cutoff<-cutoff.table$test.values[which.min(cutoff.table$Se.equals.Sp)]
  test.best.cutoff<-cbind(best.cutoff,non.parametric[
           which.min(cutoff.table$Se.equals.Sp),2:10])
  best.cutoff<-cutoff.table$test.values[which.max(cutoff.table$Efficiency)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,non.parametric[
           which.max(cutoff.table$Efficiency),2:10]))
  best.cutoff<-cutoff.table$test.values[which.min(cutoff.table$MCT)]
  test.best.cutoff<-rbind(test.best.cutoff,cbind(best.cutoff,non.parametric[
           which.min(cutoff.table$MCT),2:10]))
  rownames(test.best.cutoff)<- c("Se=Sp","Max. Efficiency","Min. MCT")
  best.cutoff<-test.best.cutoff
  rm(test.best.cutoff)
  non.parametric.inconclusive<-non.parametric[which.min(abs(
      Inconclusive-non.parametric$Sensitivity)),1:10]
  non.parametric.inconclusive<-rbind(non.parametric.inconclusive,non.parametric[
      which.min(abs(Inconclusive-non.parametric$Specificity)),1:10])
  rownames(non.parametric.inconclusive)<-c("Lower inconclusive","Upper inconclusive")
  if(Plot==TRUE){
    plot(non.parametric$test.values,non.parametric$Sensitivity,
         ylim=c(0,1),type="l",col=2, xlab="test scale",ylab="Sensitivity & Specificity",
         lty=1)
    lines(non.parametric$test.values,non.parametric$Specificity,col=4,lty=2)
    leg.txt<-c("Se", "Sp")      
    fill.col<-c(2,4)
    line.type<-c(1,2)
    if(Plot.inc.range==TRUE){
      abline(h=Inconclusive,col="lightgray",lty=4)
      abline(v=(non.parametric.inconclusive[1,1]),col="lightgray",lty=4)
      abline(v=(non.parametric.inconclusive[2,1]),col="lightgray",lty=4)
      leg.txt<-c(leg.txt,c("Inc limits"))      
      fill.col<-c(fill.col,c("lightgray"))
      line.type<-c(line.type,4)
      subtitle<-paste("Non-parametric inconclusive limits at",formatC(inc),"level:",formatC(non.parametric.inconclusive[1,1]),"-",formatC(non.parametric.inconclusive[2,1]),".")
      }
    if(Plot.Cl==TRUE){
    lines(non.parametric$test.values,non.parametric$Se.inf.Cl,lty=5,col=2)
    lines(non.parametric$test.values,non.parametric$Se.sup.Cl,lty=5,col=2)
    lines(non.parametric$test.values,non.parametric$Sp.inf.Cl,lty=3,col=4)
    lines(non.parametric$test.values,non.parametric$Sp.sup.Cl,lty=3,col=4)
    leg.txt<-c(leg.txt,c("Se conf band","Sp conf band"))      
    fill.col<-c(fill.col,c(2,4))
    line.type<-c(line.type,5,3)
    }
    if(Plot.cutoff=="Se=Sp"){
          abline(v=(best.cutoff[1,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Se=Sp:",formatC(best.cutoff[1,1])))
          }
    if(Plot.cutoff=="Max.Efficiency"){
          abline(v=(best.cutoff[2,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Max. Efficiency:",formatC(best.cutoff[2,1]),"."))
          #"Pop. prevalence:",formatC(pop.prevalence)))) Does not fit in the graph
          }
    if(Plot.cutoff=="Min.MCT"){
          abline(v=(best.cutoff[3,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Minimizing MCT:",formatC(best.cutoff[3,1],".")))
          #,"Pop. prevalence:",formatC(pop.prevalence),"Cost FN/FP:",formatC(cost))) Does not fit in the graph
          }
    legend("right",legend=leg.txt,col=fill.col, lty=line.type, bty="n")
    title(sub=subtitle,cex.sub=cex.sub)         
  }          
  rm(test.values)
  if(non.parametric.inconclusive[1,1]>=non.parametric.inconclusive[2,1])
     {
     warning("Lower inconclusive limit is higher than upper inconclusive limit.")
     }
  reteval<-list(sample.size=sample.size,sample.prevalence=sample.prevalence,pop.prevalence=pop.prevalence,cost=cost,test.summary=test.summary,inc=inc,conf.limit=conf.limit,non.parametric=non.parametric,non.parametric.inconclusive=non.parametric.inconclusive,best.cutoff=best.cutoff,cutoff.table=cutoff.table)
  class(reteval)<-"TGROC"
  if(Print==TRUE)  {print(reteval)}
  invisible(reteval)
}