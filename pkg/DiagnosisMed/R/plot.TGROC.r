#' @export
plot.TGROC<-function(x,...,
                Plot="Both",
                Plot.inc.range=TRUE,
                Plot.Cl=FALSE,
                Plot.cutoff="None",
                cex=0.5,
                cex.sub=0.85){
    if(Plot!="Both" & Plot!="Non-parametric" & Plot!="Parametric" & Plot!="None"){
      stop("Plot must be set either to 'None','Both','Non-parametric' or 'Parametric'!")
    }
    if(Plot.cutoff!="Min.MCT" & Plot.cutoff!="Se=Sp" & Plot.cutoff!="Max.Efficiency" & Plot.cutoff!="None"){
      stop("Plot.cutoff must be set either to 'None','Max.Efficiency','Min.MCT' or 'Se=Sp'!")
    }            
    if(Plot!="None"){
    if(Plot=="Parametric"|Plot=="Both"){
      plot(x$parametric$test.values,x$parametric$Sensitivity,ylim=c(0,1),type="l",col=2, xlab="Test scale",
              ylab="Sensitivity & Specificity",cex=cex,lty=1)
      lines(x$parametric$test.values,x$parametric$Specificity,col=4,type="l",lty=2,cex=cex)
        if(Plot=="Both"){
             lines(x$non.parametric$test.values,x$non.parametric$Sensitivity,col=2,type="o",lty=1,cex=cex)
             lines(x$non.parametric$test.values,x$non.parametric$Specificity,col=4,type="o",lty=2,cex=cex)
        }          
      leg.txt<-c("Se", "Sp")      
      fill.col<-c(2,4)
      line.type<-c(1,2)
      subtitle<-""
    }
    if(Plot=="Non-parametric"){
       plot(x$non.parametric$test.values,x$non.parametric$Sensitivity,
         ylim=c(0,1),type="o",col=2, xlab="test scale",ylab="Sensitivity & Specificity",
         lty=1,cex=cex)
       lines(x$non.parametric$test.values,x$non.parametric$Specificity,col=4,type="o",lty=2,cex=cex)
      leg.txt<-c("Se", "Sp")      
      fill.col<-c(2,4)
      line.type<-c(1,2)
      subtitle<-""
    }  
    if(Plot.inc.range==TRUE){
      abline(h=x$inc,col="lightgray",lty=4)
      if(Plot=="Parametric"|Plot=="Both"){
            abline(v=(x$parametric.inconclusive[1,1]),col="lightgray",lty=4)
            abline(v=(x$parametric.inconclusive[2,1]),col="lightgray",lty=4)
            subtitle<-paste("Parametric inconclusive limits at",formatC(x$inc),"level:",formatC(x$parametric.inconclusive[1,1]),
                      "-",formatC(x$parametric.inconclusive[2,1]),".")
      }
      if(Plot=="Non-parametric"){
            abline(v=(x$non.parametric.inconclusive[1,1]),col="lightgray",lty=4)
            abline(v=(x$non.parametric.inconclusive[2,1]),col="lightgray",lty=4)
            subtitle<-paste("Non-parametric inconclusive limits at",formatC(x$inc),"level:",formatC(x$non.parametric.inconclusive[1,1]),
                      "-",formatC(x$non.parametric.inconclusive[2,1]),".")
      }
      leg.txt<-c(leg.txt,c("Inc limits"))      
      fill.col<-c(fill.col,c("lightgray"))
      line.type<-c(line.type,4)
    }
    if(Plot.Cl==TRUE){
      if(Plot=="Both"|Plot=="Parametric"){
          lines(x$parametric$test.values,x$parametric$Se.inf.cl,lty=5,col=2)
          lines(x$parametric$test.values,x$parametric$Se.sup.cl,lty=5,col=2)
          lines(x$parametric$test.values,x$parametric$Sp.inf.cl,lty=3,col=4)
          lines(x$parametric$test.values,x$parametric$Sp.sup.cl,lty=3,col=4)      
      }
      if(Plot=="Non-parametric"){
          lines(x$non.parametric$test.values, x$non.parametric$Se.inf.cl,lty=5,col=2)
          lines(x$non.parametric$test.values, x$non.parametric$Se.sup.cl,lty=5,col=2)
          lines(x$non.parametric$test.values, x$non.parametric$Sp.inf.cl,lty=3,col=4)
          lines(x$non.parametric$test.values, x$non.parametric$Sp.sup.cl,lty=3,col=4)            
      }  
    leg.txt<-c(leg.txt,c("Se conf. band","Sp conf. band"))      
    fill.col<-c(fill.col,c(2,4))
    line.type<-c(line.type,5,3)
    }
    if(Plot.cutoff=="Se=Sp"){
      if(Plot=="Both"|Plot=="Parametric"){
          abline(v=(x$par.test.best.cutoff[1,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by parametric Se=Sp:",formatC(x$par.test.best.cutoff[1,1])))      
      }
      if(Plot=="Non-parametric"){
          abline(v=(x$np.test.best.cutoff[1,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Se=Sp:",formatC(x$np.test.best.cutoff[1,1])))      
      }

    }
    if(Plot.cutoff=="Max.Efficiency"){
       if(Plot=="Both"|Plot=="Parametric"){
          abline(v=(x$par.test.best.cutoff[2,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by parametric Max. Efficiency:",formatC(x$par.test.best.cutoff[2,1]),"."))
          #"Pop. prevalence:",formatC(pop.prevalence)))) Does not fit in the graph       
       }
       if(Plot=="Non-parametric"){
          abline(v=(x$np.test.best.cutoff[2,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Max. Efficiency:",formatC(x$np.test.best.cutoff[2,1]),"."))
          #"Pop. prevalence:",formatC(pop.prevalence)))) Does not fit in the graph       
       }
    }
    if(Plot.cutoff=="Min.MCT"){
       if(Plot=="Both"|Plot=="Parametric"){
          abline(v=(x$par.test.best.cutoff[3,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by minimizing parametric MCT:",formatC(x$np.test.best.cutoff[3,1]),"."))
          #,"Pop. prevalence:",formatC(pop.prevalence),"Cost FN/FP:",formatC(cost))) Does not fit in the graph
       }
       if(Plot=="Non-parametric"){
          abline(v=(x$np.test.best.cutoff[3,1]),col="lightgray",lty=6)
          leg.txt<-c(leg.txt,c("Best cut-off"))      
          fill.col<-c(fill.col,c("lightgray"))
          line.type<-c(line.type,6)
          subtitle<-paste(subtitle,paste("Cut-off estimated by Minimizing MCT:",formatC(x$np.test.best.cutoff[3,1],".")))
          #,"Pop. prevalence:",formatC(pop.prevalence),"Cost FN/FP:",formatC(cost))) Does not fit in the graph
       }
    }
    legend("right",legend=leg.txt,col=fill.col,lty=line.type, bty="n")
    title(sub=subtitle,cex.sub=cex.sub)         
  }
}          