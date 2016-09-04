#' @export
# the plot commands
plot.ROC<-function(x,Plot.point="Min.ROC.Dist",cex.sub=.85,p.cex=1,...){

    if(Plot.point!="None" & Plot.point!="Min.ROC.Dist" & Plot.point!="Max.Accuracy" &
       Plot.point!="Max.DOR" & Plot.point!="Error.rate" & Plot.point!="Max.Accuracy.area" &
       Plot.point!="Max.Sens+Spec" & Plot.point!="Max.Youden" & Plot.point!="Se=Sp" & 
       Plot.point!="Min.ROC.Dist" & Plot.point!="Max.Efficiency" & Plot.point!="Min.MCT"){
     stop("The Plot.point option is not correctly set! Type '?ROC' to check possible options.")}  
 
    plot(1-x$test.diag.table$Specificity,x$test.diag.table$Sensitivity,type="l",
      col=1,xlab="1-Specificity",ylab="Sensitivity",xlim=c(0,1),ylim=c(0,1))
    grid()
    segments(0,0,1,1,col="lightgray")
    
    if(Plot.point=="None"){
      legend("bottomright",legend=(        
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       ),bty="n")}
       
    if(Plot.point=="Max.Accuracy")
      {points(1-x$test.best.cutoff[1,5],x$test.best.cutoff[1,2],col=1,pch=19,cex=p.cex)
      title(sub="Cut-off estimated by maximazing accuracy.",cex.sub=cex.sub)   
      legend("bottomright",legend=(c(paste("cut off:",formatC(x$test.best.cutoff[1,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[1,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[1,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}
    
    if(Plot.point=="Max.DOR")
      {points(1-x$test.best.cutoff[2,5],x$test.best.cutoff[2,2],col=1,pch=19,cex=p.cex)
       title(sub="Cut-off estimated by maximazing diagnostic odds ratio.",cex.sub=cex.sub)  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[2,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[2,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[2,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}
    
    if(Plot.point=="Error.rate")
      {points(1-x$test.best.cutoff[3,5],x$test.best.cutoff[3,2],col=1,pch=19,cex=p.cex)
       title(sub="Cut-off estimated by minimizing error rate.",cex.sub=cex.sub)  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[3,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[3,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[3,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}
       
    if(Plot.point=="Max.Accuracy.area")
      {points(1-x$test.best.cutoff[4,5],x$test.best.cutoff[4,2],col=1,pch=19,cex=p.cex)
       title(sub="Cut-off estimated by maximazing the area related to accuracy.",cex.sub=cex.sub)  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[4,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[4,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[4,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}

    if(Plot.point=="Max.Sens+Spec")
      {points(1-x$test.best.cutoff[5,5],x$test.best.cutoff[4,2],col=1,pch=19,cex=p.cex)
      title(sub="Cut-off value where the sum Se + Sp is maximized.",cex.sub=cex.sub)         
      legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[5,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[5,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[5,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}

    if(Plot.point=="Max.Youden")
      {points(1-x$test.best.cutoff[6,5],x$test.best.cutoff[6,2],
         col=1,pch=19,cex=p.cex)
       title(sub="Cut-off estimated by maximazing Youden Index.",cex.sub=cex.sub)         
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[6,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[6,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[6,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}

    if(Plot.point=="Se=Sp")
      {points(1-x$test.best.cutoff[7,5],x$test.best.cutoff[7,2],col=1,pch=19,cex=p.cex)
       title(sub="Cut-off value where Se is the closest to Sp.",cex.sub=cex.sub)
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[7,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[7,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[7,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
        )),bty="n")}

    if(Plot.point=="Min.ROC.Dist")
      {points(1-x$test.best.cutoff[8,5],x$test.best.cutoff[8,2],col=1,pch=19,cex=p.cex)
       title(sub="Cut-off that minimizes the distance between the curve and upper left corner.",cex.sub=cex.sub)         
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[8,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[8,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[8,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}

    if(Plot.point=="Max.Efficiency")
      {points(1-x$test.best.cutoff[9,5],x$test.best.cutoff[9,2],col=1,pch=19,cex=p.cex)
       title(sub=paste("Cut-off maximizing efficiency: population prevalence =",
             formatC(x$pop.prevalence,digits=2)),cex.sub=cex.sub)         
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[9,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[9,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[9,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}

    if(Plot.point=="Min.MCT")
      {points(1-x$test.best.cutoff[10,5],x$test.best.cutoff[10,2],col=1,pch=19,cex=p.cex)
       title(sub=paste("Cut-off minimazing MCT: population prevalence =",
             formatC(x$pop.prevalence,digits=2,format="f"),"; cost(FN)/cost(FP) =",
             formatC(x$cost,digits=2)),cex.sub=cex.sub)                  
       legend("bottomright",legend=(c(
         paste("cut off:",formatC(x$test.best.cutoff[10,1],digits=4)),
         paste("Sensitivity:",formatC(x$test.best.cutoff[10,2],digits=4,format="f")),
         paste("Specificity:",formatC(x$test.best.cutoff[10,5],digits=4,format="f")),
         paste("AUC:",formatC(x$AUC.summary[2],digits=4,format="f"))
       )),bty="n")}
}