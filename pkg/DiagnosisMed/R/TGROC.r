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
  if(Plot!="Both" & Plot!="Non-parametric" & Plot!="Parametric" & Plot!="None"){
      stop("Plot must be set either to 'None','Both','Non-parametric' or 'Parametric'!")
  }
  if(Plot.cutoff!="Min.MCT" & Plot.cutoff!="Se=Sp" & Plot.cutoff!="Max.Efficiency" & Plot.cutoff!="None"){
      stop("Plot.cutoff must be set either to 'None','Max.Efficiency','Min.MCT' or 'Se=Sp'!")
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
    non.parametric$Sensitivity[i] <- round(sum(test.table[(i:nrow(test.table)),2])/sum(test.table[,2]),digits=4)
    non.parametric$Se.inf.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[i:nrow(test.table),2]),sum(test.table[,2]),
        conf.level=CL)[4]),digits=4)
    non.parametric$Se.sup.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[i:nrow(test.table),2]),sum(test.table[,2]),
        conf.level=CL)[5]),digits=4)
    non.parametric$Specificity[i] <- round((sum(test.table[(1:i-1),1]))/(sum(test.table[,1])),digits=4)
    non.parametric$Sp.inf.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[(1:i-1),1]),sum(test.table[,1]),
        conf.level=CL)[4]),digits=4)
    non.parametric$Sp.sup.Cl[i]<-round(as.numeric(binom.wilson(sum(test.table[(1:i-1),1]),sum(test.table[,1]),
        conf.level=CL)[5]),digits=4)
  }
  non.parametric$PLR<-round(non.parametric$Sensitivity/(1-non.parametric$Specificity),digits=2)
  non.parametric$PLR.inf.cl<-round(exp(log(non.parametric$PLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-non.parametric$Sensitivity)/
    ((sum(non.parametric[,2]))*non.parametric$Specificity)+(non.parametric$Specificity)/((sum(non.parametric[,1]))*
    (1-non.parametric$Specificity)))),digits=2)
  non.parametric$PLR.sup.cl<-round(exp(log(non.parametric$PLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt(
    (1-non.parametric$Sensitivity)/((sum(non.parametric[,2]))*non.parametric$Specificity)+(
    non.parametric$Specificity)/((sum(non.parametric[,1]))*(1-non.parametric$Specificity)))),digits=2)
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
       non.parametric[which.min(non.parametric$Se.equals.Sp),1:10],
       non.parametric[which.max(non.parametric$Efficiency),1:10],
       non.parametric[which.min(non.parametric$MCT),1:10]))
  
  #best.cutoff,non.parametric[which.min(non.parametric$Se.equals.Sp),1:10])
  #best.cutoff<-non.parametric$test.values[which.max(non.parametric$Efficiency)]
  #np.test.best.cutoff<-rbind(np.test.best.cutoff,cbind(best.cutoff,non.parametric[
  #         which.max(non.parametric$Efficiency),2:10]))
  #best.cutoff<-non.parametric$test.values[which.min(non.parametric$MCT)]
  #np.test.best.cutoff<-rbind(np.test.best.cutoff,cbind(best.cutoff,non.parametric[
  #         which.min(non.parametric$MCT),2:10]))
  rownames(np.test.best.cutoff)<- c("Se=Sp","Max. Efficiency","Min. MCT")

  non.parametric.inconclusive<-as.data.frame(rbind(
      non.parametric[which.min(abs(Inconclusive-non.parametric$Sensitivity)),1:10],
      non.parametric[which.min(abs(Inconclusive-non.parametric$Specificity)),1:10]))
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
  parametric$Se.inf.Cl<-parametric$Sensitivity - qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Sensitivity * 
        (1 - parametric$Sensitivity))/(sample.size * sample.prevalence)))
  parametric$Se.sup.Cl<-parametric$Sensitivity + qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Sensitivity * 
        (1 - parametric$Sensitivity))/(sample.size * sample.prevalence)))
  parametric$Specificity <- as.numeric(sim(net.Sp$net, test.values))
  parametric$Sp.inf.Cl <- parametric$Specificity - qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Specificity * 
        (1 - parametric$Specificity))/(sample.size * (1-sample.prevalence))))
  parametric$Sp.sup.Cl <- parametric$Specificity + qnorm(1 - (1-conf.limit)/2) * sqrt(((parametric$Specificity * 
        (1 - parametric$Specificity))/(sample.size * (1-sample.prevalence))))
  parametric$PLR <- parametric$Sensitivity/(1-parametric$Specificity) 
  parametric$PLR.inf.Cl<-exp(log(parametric$PLR)-(qnorm(1-((1-conf.limit)/2),mean=0,sd=1))*sqrt((1-parametric$Sensitivity)/
        ((sample.size * sample.prevalence)*parametric$Specificity)+(parametric$Specificity)/((sample.size * 
        (1-sample.prevalence))*(1-parametric$Specificity))))
  parametric$PLR.sup.Cl<-exp(log(parametric$PLR)+(qnorm(1-((1-conf.limit)/2),mean=0,sd=1))*sqrt((1-parametric$Sensitivity)/
        ((sample.size * sample.prevalence)*parametric$Specificity)+(parametric$Specificity)/((sample.size * 
        (1-sample.prevalence))*(1-parametric$Specificity))))            
  #parametric$NLR <- (1-parametric$Specificity)/parametric$Sensitivity
  #parametric$NLR.inf.Cl <- exp(log(parametric$NLR)-(qnorm(1-((1-(1-conf.limit))/2),mean=0,sd=1))*
  #     sqrt((parametric$Sensitivity)/((sample.size * sample.prevalence)*(1-parametric$Sensitivity))+(1-parametric$Specificity)/
  #     ((sample.size * (1-sample.prevalence))*(parametric$Specificity))))
  #parametric$NLR.sup.Cl <- exp(log(parametric$NLR)+(qnorm(1-((1-conf.limit)/2),mean=0,sd=1))*sqrt((parametric$Sensitivity)/
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
       
# Plot Commands
    if(Plot!="None"){
    if(Plot=="Parametric"|Plot=="Both"){
      plot(parametric$test.values,parametric$Sensitivity,ylim=c(0,1),type="l",col=2, xlab="Test scale",
              ylab="Sensitivity & Specificity",cex=cex,lty=1)
      lines(parametric$test.values,parametric$Specificity,col=4,type="l",lty=2,cex=cex)
        if(Plot=="Both"){
             lines(non.parametric$test.values,non.parametric$Sensitivity,col=2,type="o",lty=1,cex=cex)
             lines(non.parametric$test.values,non.parametric$Specificity,col=4,type="o",lty=2,cex=cex)
        }          
      leg.txt<-c("Se", "Sp")      
      fill.col<-c(2,4)
      line.type<-c(1,2)
      subtitle<-""
    }
    if(Plot=="Non-parametric"){
       plot(non.parametric$test.values,non.parametric$Sensitivity,
         ylim=c(0,1),type="o",col=2, xlab="test scale",ylab="Sensitivity & Specificity",
         lty=1,cex=cex)
       lines(non.parametric$test.values,non.parametric$Specificity,col=4,type="o",lty=2,cex=cex)
      leg.txt<-c("Se", "Sp")      
      fill.col<-c(2,4)
      line.type<-c(1,2)
      subtitle<-""
    }  
    if(Plot.inc.range==TRUE){
      abline(h=Inconclusive,col="lightgray",lty=4)
      if(Plot=="Parametric"|Plot=="Both"){
            abline(v=(parametric.inconclusive[1,1]),col="lightgray",lty=4)
            abline(v=(parametric.inconclusive[2,1]),col="lightgray",lty=4)
            subtitle<-paste("Parametric inconclusive limits at",formatC(inc),"level:",formatC(parametric.inconclusive[1,1]),
                      "-",formatC(parametric.inconclusive[2,1]),".")
      }
      if(Plot=="Non-parametric"){
            abline(v=(non.parametric.inconclusive[1,1]),col="lightgray",lty=4)
            abline(v=(non.parametric.inconclusive[2,1]),col="lightgray",lty=4)
            subtitle<-paste("Non-parametric inconclusive limits at",formatC(inc),"level:",formatC(non.parametric.inconclusive[1,1]),
                      "-",formatC(non.parametric.inconclusive[2,1]),".")
      }
      leg.txt<-c(leg.txt,c("Inc limits"))      
      fill.col<-c(fill.col,c("lightgray"))
      line.type<-c(line.type,4)
    }
    if(Plot.Cl==TRUE){
      if(Plot=="Both"|Plot=="Parametric"){
          lines(parametric$test.values,parametric$Se.inf.Cl,lty=5,col=2)
          lines(parametric$test.values,parametric$Se.sup.Cl,lty=5,col=2)
          lines(parametric$test.values,parametric$Sp.inf.Cl,lty=3,col=4)
          lines(parametric$test.values,parametric$Sp.sup.Cl,lty=3,col=4)      
      }
      if(Plot=="Non-parametric"){
          lines(non.parametric$test.values,non.parametric$Se.inf.Cl,lty=5,col=2)
          lines(non.parametric$test.values,non.parametric$Se.sup.Cl,lty=5,col=2)
          lines(non.parametric$test.values,non.parametric$Sp.inf.Cl,lty=3,col=4)
          lines(non.parametric$test.values,non.parametric$Sp.sup.Cl,lty=3,col=4)            
      }  
    leg.txt<-c(leg.txt,c("Se conf. band","Sp conf. band"))      
    fill.col<-c(fill.col,c(2,4))
    line.type<-c(line.type,5,3)
    }
    if(Plot.cutoff=="Se=Sp"){
      if(Plot=="Both"|Plot=="Parametric"){
          abline(v=(par.test.best.cutoff[1,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by parametric Se=Sp:",formatC(par.test.best.cutoff[1,1])))      
      }
      if(Plot=="Non-parametric"){
          abline(v=(np.test.best.cutoff[1,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Se=Sp:",formatC(np.test.best.cutoff[1,1])))      
      }

    }
    if(Plot.cutoff=="Max.Efficiency"){
       if(Plot=="Both"|Plot=="Parametric"){
          abline(v=(par.test.best.cutoff[2,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by parametric Max. Efficiency:",formatC(par.test.best.cutoff[2,1]),"."))
          #"Pop. prevalence:",formatC(pop.prevalence)))) Does not fit in the graph       
       }
       if(Plot=="Non-parametric"){
          abline(v=(np.test.best.cutoff[2,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Max. Efficiency:",formatC(np.test.best.cutoff[2,1]),"."))
          #"Pop. prevalence:",formatC(pop.prevalence)))) Does not fit in the graph       
       }
    }
    if(Plot.cutoff=="Min.MCT"){
       if(Plot=="Both"|Plot=="Parametric"){
          abline(v=(par.test.best.cutoff[3,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by minimizing paramaetric MCT:",formatC(np.test.best.cutoff[3,1],".")))
          #,"Pop. prevalence:",formatC(pop.prevalence),"Cost FN/FP:",formatC(cost))) Does not fit in the graph
       }
       if(Plot=="Non-parametric"){
          abline(v=(np.test.best.cutoff[3,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Minimizing MCT:",formatC(np.test.best.cutoff[3,1],".")))
          #,"Pop. prevalence:",formatC(pop.prevalence),"Cost FN/FP:",formatC(cost))) Does not fit in the graph
       }
    }
    legend("right",legend=leg.txt,col=fill.col,lty=line.type, bty="n")
    title(sub=subtitle,cex.sub=cex.sub)         
  }          
  rm(test.values)
  if(non.parametric.inconclusive[1,1]>non.parametric.inconclusive[2,1]){
     warning("Non-parametric lower inconclusive limit is higher than upper inconclusive limit.")
  }
  if(parametric.inconclusive[1,1]>parametric.inconclusive[2,1]){
     warning("Parametric lower inconclusive limit is higher than upper inconclusive limit.")
  }
  if(np.test.best.cutoff[1,1]>non.parametric.inconclusive[2,1]|
     np.test.best.cutoff[2,1]>non.parametric.inconclusive[2,1]|
     np.test.best.cutoff[3,1]>non.parametric.inconclusive[2,1]){
     warning("At least one fo the non-parametric best cut-off is higher then upper inconclusive limit.")
  }        
  if(np.test.best.cutoff[1,1]<non.parametric.inconclusive[1,1]|
     np.test.best.cutoff[2,1]<non.parametric.inconclusive[1,1]|
     np.test.best.cutoff[3,1]<non.parametric.inconclusive[1,1]){
     warning("At least one fo the non-parametric best cut-off is lower then lower inconclusive limit.")
  }
  if(par.test.best.cutoff[1,1]>parametric.inconclusive[2,1]|
     par.test.best.cutoff[2,1]>parametric.inconclusive[2,1]|
     par.test.best.cutoff[3,1]>parametric.inconclusive[2,1]){
     warning("At least one fo the parametric best cut-off is higher then upper inconclusive limit.")
  }        
  if(par.test.best.cutoff[1,1]<parametric.inconclusive[1,1]|
     par.test.best.cutoff[2,1]<parametric.inconclusive[1,1]|
     par.test.best.cutoff[3,1]<parametric.inconclusive[1,1]){
     warning("At least one fo the parametric best cut-off is lower then lower inconclusive limit.")
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
  invisible(reteval)
}