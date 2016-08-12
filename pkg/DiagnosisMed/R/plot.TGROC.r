#' @export
plot.TGROC <- function(x,...,
                Plot.type = c("TGROC","ROC","None"),
                Plot = c("Binormal","Non-parametric","NN-parametric"),
                Plot.inc.area = TRUE,
                Plot.Cl = FALSE,
                Plot.trheshold = "None",
                threshold.arg = list(col = gray(.5), lty = 6),
                cex.sub = 0.85,
                ylab = "Sensitivity & Specificity",
                xlab = "Test scale",
                ylim = c(0,1),
                xlim = "auto",
                shade.args = list(col = gray(.8), density = 45, border = NA),
                np.Se.args = list(type = "o", col = "blue", lty = 1, cex  = .5),
                np.Se.ci.args = list(lty = 5, col = "blue"),
                np.Sp.args = list(type = "o", col = "red", lty = 2,cex = .5),
                np.Sp.ci.args = list(lty = 3, col = "red"),
                NN.Se.args = list(type = "l", col = "blue", lty = 1, cex = .5),
                NN.Se.ci.args = list(lty = 5, col = "blue"),
                NN.Sp.args = list(type = "l", col = "red", lty = 2, cex = .5),
                NN.Sp.ci.args = list(lty = 3, col = "red"),
                BN.Se.args = list(type = "l", col = "blue", lty = 1, cex = .5),
                BN.Se.ci.args = list(lty = 5, col = "blue"),
                BN.Sp.args = list(type = "l", col = "red", lty = 2, cex = .5),
                BN.Sp.ci.args = list(lty = 3, col = "red"),
                auto.legend = TRUE,
                legend.args = list(x = "top", border = NA, bty = "n", xpd = NA, inset = -.18, ncol = 2, cex = .8)){
  
  # Settign a warnings for valid values
  if(any(!(Plot.type %in% c("TGROC","ROC","None")))){
    stop("The allowed values for Plot.type argument are: 'TGROC','ROC', or 'None'.")
  }
  if(any(!(Plot %in% c("Non-parametric","NN-parametric","Binormal","None")))){
      stop("The allowed values for Plot argument are: 'None','Non-parametric', or 'NN-parametric'!")
  }
  if(any(!(Plot.trheshold %in% c("Min MCT", "Se = Sp", "Max Efficiency", "Min ROC distance", "Min Error rate", "Max DOR", "Max Accuracy area", "Max Accuracy", "Max Youden", "None")))){
      stop("Plot.trheshold argument does not have a valid input. Check documentation.")
  }
  if(!is.logical(Plot.inc.area)){
    stop("Plot.inc.area argument must be logical.")
  }
  if(!is.logical(Plot.Cl)){
    stop("Plot.Cl argument must be logical.")
  }
  if(!is.logical(auto.legend)){
    stop("legend argument must be logical.")
  }
  
  # Settings for TGROC plot   
  if(Plot.type[1] == "TGROC"){
    # setting the xlim argument
    if(any(Plot[1] == "Non-parametric")){
      if(xlim == "auto"){
          xlim <- range(c(x$SS$test.values))
      }
    }    
    if(any(Plot[1] == "NN-parametric")){
      if(xlim == "auto"){
        xlim <- range(c(x$SS$test.values,x$NN.SS$test.values))
      }
    }
    if(any(Plot[1] == "Binormal")){
      if(xlim == "auto"){
        xlim <- range(c(x$SS$test.values,x$BN.SS$test.values))
      }
    }
    
    # Opening a blank device
    plot(0,0, type="n", ..., xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim)
    # plot(0,0,type="n", xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim)
    
    # Ploting the inconclusive shade area ----------------------------------------
    if(Plot.inc.area){
    if(Plot[1] == "Non-parametric"){
      shade.args$x <- c(rep(x$np.inconclusive["Lower inconclusive","test.values"],2),rep(x$np.inconclusive["Upper inconclusive","test.values"],2))
      shade.args$y <- c(0,x$inc,x$inc,0)
      do.call(polygon,shade.args)
      subtitle <- paste0("Non-parametric inconclusive limits at ",formatC(x$inc)," level: ",formatC(x$np.inconclusive[1,1]),"-",formatC(x$np.inconclusive[2,1]),".")
    }
    if(Plot[1] == "NN-parametric"){
      shade.args$x <- c(rep(x$NN.inconclusive["Lower inconclusive","test.values"],2),rep(x$NN.inconclusive["Upper inconclusive","test.values"],2))
      shade.args$y <- c(0,x$inc,x$inc,0)
      do.call(polygon,shade.args)
      subtitle <- paste0("Parametric (neural network) inconclusive limits at",formatC(x$inc)," level: ",formatC(x$NN.inconclusive[1,1]),"-",formatC(x$NN.inconclusive[2,1]),".")
    }
    if(Plot[1] == "Binormal"){
      shade.args$x <- c(rep(x$BN.inconclusive["Lower inconclusive","test.values"],2),rep(x$BN.inconclusive["Upper inconclusive","test.values"],2))
      shade.args$y <- c(0,x$inc,x$inc,0)
      do.call(polygon,shade.args)
      subtitle <- paste0("Parametric (Binormal) inconclusive limits at",formatC(x$inc)," level: ",formatC(x$BN.inconclusive[1,1]),"-",formatC(x$BN.inconclusive[2,1]),".")
    }
  }

    # Ploting the Se and Sp lines
    if(any(Plot == "Non-parametric")){
      np.Se.args$x <- x$SS$test.values
      np.Se.args$y <- x$SS$Sensitivity
      do.call(lines,np.Se.args)
 
      np.Sp.args$x <- x$SS$test.values
      np.Sp.args$y <- x$SS$Specificity
      do.call(lines,np.Sp.args)
    }
    if(any(Plot == "NN-parametric")){
      NN.Se.args$x <- x$NN.SS$test.values
      NN.Se.args$y <- x$NN.SS$Sensitivity
      do.call(lines, NN.Se.args)
      
      NN.Sp.args$x <- x$NN.SS$test.values
      NN.Sp.args$y <- x$NN.SS$Specificity
      do.call(lines, NN.Sp.args)
    }
    if(any(Plot == "Binormal")){
      BN.Se.args$x <- x$BN.SS$test.values
      BN.Se.args$y <- x$BN.SS$Sensitivity
      do.call(lines, BN.Se.args)
      
      BN.Sp.args$x <- x$BN.SS$test.values
      BN.Sp.args$y <- x$BN.SS$Specificity
      do.call(lines, BN.Sp.args)
    }
    
    # Ploting the confidence bands -------------------------------------
    if(Plot.Cl){
      if(Plot[1] == "Non-parametric"){
        np.Se.ci.args$x <- x$SS$test.values
        np.Se.ci.args$y <- x$SS$Se.inf.cl
        do.call(lines, np.Se.ci.args)
        np.Se.ci.args$y <- x$SS$Se.sup.cl
        do.call(lines, np.Se.ci.args)
        
        np.Sp.ci.args$x <- x$SS$test.values
        np.Sp.ci.args$y <- x$SS$Sp.inf.cl
        do.call(lines, np.Sp.ci.args)
        np.Sp.ci.args$y <- x$SS$Sp.sup.cl
        do.call(lines, np.Sp.ci.args)
      }
      
      if(Plot[1] == "NN-parametric"){
        NN.Se.ci.args$x <- x$NN.SS$test.values
        NN.Se.ci.args$y <- x$NN.SS$Se.inf.cl
        do.call(lines, NN.Se.ci.args)
        NN.Se.ci.args$y <- x$NN.SS$Se.sup.cl
        do.call(lines, NN.Se.ci.args)
        
        NN.Sp.ci.args$x <- x$NN.SS$test.values
        NN.Sp.ci.args$y <- x$NN.SS$Sp.inf.cl
        do.call(lines, NN.Sp.ci.args)
        NN.Sp.ci.args$y <- x$NN.SS$Sp.sup.cl
        do.call(lines, NN.Sp.ci.args)
      }
      
      if(Plot[1] == "Binormal"){
        BN.Se.ci.args$x <- x$BN.SS$test.values
        BN.Se.ci.args$y <- x$BN.SS$Se.inf.cl
        do.call(lines, BN.Se.ci.args)
        BN.Se.ci.args$y <- x$BN.SS$Se.sup.cl
        do.call(lines, BN.Se.ci.args)
        
        BN.Sp.ci.args$x <- x$BN.SS$test.values
        BN.Sp.ci.args$y <- x$BN.SS$Sp.inf.cl
        do.call(lines, BN.Sp.ci.args)
        BN.Sp.ci.args$y <- x$BN.SS$Sp.sup.cl
        do.call(lines, BN.Sp.ci.args)
      }
    }

    # Ploting the best threshold vertical line ----------------------------
    if(Plot.trheshold != "None"){
      if(Plot[1] == "Non-parametric"){
        threshold.arg$v <- x$np.best.threshold[Plot.trheshold, 1]
        do.call(abline, threshold.arg)
        subtitle <- paste(subtitle,paste("Threshold estimated by non-parametric", Plot.threshold,":",formatC(x$np.best.threshold[1,1])))
      }
      if(Plot[1] == "Binormal"){
        threshold.arg$v <- x$BN.best.threshold[Plot.trheshold, 1]
        do.call(abline, threshold.arg)
        subtitle <- paste(subtitle,paste("Threshold estimated by parametric (Binormal)", Plot.threshold,":",formatC(x$BN.best.threshold[Plot.trheshold,1])))
      }
      if(Plot[1] == "NN-parametric"){
        threshold.arg$v <- x$NN.best.threshold[Plot.trheshold, 1]
        do.call(abline, threshold.arg)
        subtitle <- paste(subtitle,paste("Threshold estimated by parametric (neural network)", Plot.threshold,":",formatC(x$NN.best.threshold[Plot.threshold,1])))
      }
    }

    #Setting the legend arguments
    if (auto.legend) {

      legend.args$legend <- c("Se", "Sp")
      legend.args$fill <- c(NA, NA)
      legend.args$density <- c(NA, NA)
      legend.args$col <- c(BN.Se.args$col, BN.Sp.args$col)
      legend.args$lty <- c(BN.Se.args$lty, BN.Sp.args$lty)
      legend.args$density <- c(NA, NA)

      if (Plot.Cl) {
        legend.args$legend <- c(legend.args$legend, "Se conf band", "Sp conf band") 
        legend.args$fill <- c(legend.args$fill, NA, NA)
        legend.args$density <- c(legend.args$density, NA, NA)
        legend.args$col <- c(legend.args$col, BN.Se.ci.args$col, BN.Sp.ci.args$col)
        legend.args$lty <- c(legend.args$lty, BN.Se.ci.args$lty, BN.Sp.ci.args$lty)
      }

      if (Plot.inc.area) {
        legend.args$legend <- c(legend.args$legend, "Inc area") 
        legend.args$fill <- c(legend.args$fill, shade.args$col)
        legend.args$density <- c(legend.args$density, shade.args$density)
        legend.args$col <- c(legend.args$col, shade.args$col)
        legend.args$lty <- c(legend.args$lty, NA) 
      }
      
      if (Plot.trheshold != "None") {
        list(col = gray(.5), lty = 6)
        threshold.arg
        legend.args$legend <- c(legend.args$legend, "Best threshold") 
        legend.args$fill <- c(legend.args$fill, NA)
        legend.args$density <- c(legend.args$density, NA)
        legend.args$col <- c(legend.args$col, threshold.arg$col)
        legend.args$lty <- c(legend.args$lty, threshold.arg$lty)
        
      }
      
    do.call(legend,legend.args)
  }  
    title(sub = subtitle, cex.sub = cex.sub)         
  }
  if(Plot.type[1] == "ROC"){
    
  }
}