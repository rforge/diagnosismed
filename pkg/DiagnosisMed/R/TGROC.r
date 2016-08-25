#' @export
#' @import AMORE

TGROC <- function(ref,
                test,
                Cost = 1,
                CL = 0.95,
                Inconclusive = 0.95,
                Prevalence = 0,
                t.max = NULL,
                t.min = NULL,
                precision = 0.05,  
                n.neurons = c(1, 5, 1),
                learning.rate.global = 1e-2,
                momentum.global = 0.3,
                error.criterium = "LMS",
                Stao = NA,
                hidden.layer = "sigmoid",
                output.layer = "sigmoid",
                method = "ADAPTgdwm",
                report = FALSE,
                show.step = 5000,
                n.shows = 1,
                Plot.type = "TGROC",
                Plot = c("Binormal","Non-parametric"),
                Plot.inc.area = TRUE,
                Plot.Cl = FALSE,
                Plot.threshold = "None",
                cex = 0.5,
                cex.sub = 0.85,
                reverse = "auto"){

  # Preventing wrong inputs and warnings ---------------------------------------
  if (any(levels(as.factor(ref)) != c(0, 1))){
    stop("Your reference standard must be coded as 0 (absence) and 1 (presence). Check reference categories!")
  }
  if (is.null(precision) || !is.numeric(precision)) {
      stop("Precision must be set to a numeric value!")
  }
  if (Prevalence > 0){
    (pop.prevalence <- Prevalence)
  }

  # Making the non-parametric trade-off for Se and Sp (ROC analysis) -----------
  SeSp <- SS(ref, test, reverse = reverse, CL = CL)
  AUC <- np.auROCc(ref, test, reverse = reverse, CL = CL)
  
  # Setting requeired additional objects
  sample.prevalence <- SeSp$sample.prevalence
  sample.size <- SeSp$sample.size
  names(sample.prevalence) <- c("Condition prevalence in the sample")
  names(sample.size) <- c("Sample size")
  if (Prevalence == 0) { pop.prevalence <- sample.prevalence }
  names(pop.prevalence)<-c("Informed disease prevalence in the population.")
  cost <- Cost
  names(cost) <- c("Informed costs(FN)/costs(FP).")
  conf.limit <- CL
  inc <- Inconclusive
  names(inc) <- "Inconclusive tolerance level."

  # Overall test summary -------------------------------------------------------
  test.summary <- rbind(c(summary(test), sd(test)),
                        c(summary(test[which(ref == 1)]), sd(test[which(ref == 1)])),
                        c(summary(test[which(ref == 0)]), sd(test[which(ref == 0)])))
  colnames(test.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "SD")
  rownames(test.summary) <- c("Overall", "With the condition", "Without the condition")

  #Setting the decision Thresholds ---------------------------------------------
  np.best.threshold <- thresholds(SeSp, pop.prevalence = pop.prevalence, Cost = Cost)

  # Extracting the inconclusive ranges 
  np.inconclusive <- inc.limits(SeSp, Inconclusive = Inconclusive)
  
  # Fitting the Neural Network for smoothed threshold estimation
  NN.SeSp <- NN.SS(SeSp, t.max = t.max, t.min = t.min, precision = precision,  n.neurons = n.neurons,  learning.rate.global = learning.rate.global,  momentum.global = momentum.global,  error.criterium = error.criterium,  Stao = Stao, hidden.layer = hidden.layer,  output.layer = output.layer, method = method,  report = report,  show.step = show.step,  n.shows = n.shows,  CL = CL)
  
  # Extracting the NN decision thresholds --------------------------------------------------
  NN.best.threshold <- thresholds(NN.SeSp, pop.prevalence = pop.prevalence, Cost = Cost)
  
  # Extracting the inconclusive ranges -----------------------------------------------------
  NN.inconclusive <- inc.limits(NN.SeSp, Inconclusive = Inconclusive)
  
  # Fitting the Binormal SS object
  BN.SeSp <- BN.SS(ref = ref, test = test, CL = CL, t.max = t.max, t.min = t.min, precision = precision)
  
  #Setting the decision Thresholds --------------------------------------------------  
  BN.best.threshold <- thresholds(BN.SeSp, pop.prevalence = pop.prevalence, Cost = Cost)
  
  # Extracting the inconclusive ranges 
  BN.inconclusive <- inc.limits(BN.SeSp, Inconclusive = Inconclusive)
  
  # Warnings Regarding the thresholds out of the inconclusive range
  if(np.inconclusive[1,1] > np.inconclusive[2,1]){
     warning("Non-parametric lower inconclusive limit is higher than upper inconclusive limit. \n Either ROC analysis should be reversed or Inconclusive argument should have higher value.")
  }
  if(NN.inconclusive[1,1] > NN.inconclusive[2,1]){
     warning("NN-parametric lower inconclusive limit is higher than upper inconclusive limit. \n Either ROC analysis should be reversed or Inconclusive argument should have higher value.")
  }
  if(BN.inconclusive[1,1] > BN.inconclusive[2,1]){
    warning("Binormal lower inconclusive limit is higher than upper inconclusive limit. \n Inconclusive argument should have higher value.")
  }
  
  # if(any(np.best.threshold[,1] > np.inconclusive[2,1])){
  #    warning("At least one of the non-parametric threshold is higher then upper inconclusive limit.")
  # }        
  # if(any(np.best.threshold[,1] < np.inconclusive[1,1])){
  #    warning("At least one of the non-parametric threshold is lower then lower inconclusive limit.")
  # }
  # if(any(NN.best.threshold[,1] > NN.inconclusive[2,1])){
  #    warning("At least one of the NN-parametric best threshold is higher then upper inconclusive limit.")
  # }        
  # if(any(NN.best.threshold[,1] < NN.inconclusive[1,1])){
  #    warning("At least one of the NN-parametric best threshold is lower then lower inconclusive limit.")
  # }
  # if(any(BN.best.threshold[,1] > BN.inconclusive[2,1])){
  #   warning("At least one of the Binormal best threshold is higher then upper inconclusive limit.")
  # }        
  # if(any(BN.best.threshold[,1] < BN.inconclusive[1,1])){
  #   warning("At least one of the Binormal best threshold is lower then lower inconclusive limit.")
  # }
  output <- list(sample.size = sample.size,
                sample.prevalence = sample.prevalence,
                pop.prevalence = pop.prevalence,
                cost = cost,
                test.summary = test.summary,
                inc = inc,
                conf.limit = conf.limit,
                AUC = AUC,
                SS = SeSp$table,
                np.inconclusive = np.inconclusive,
                np.best.threshold = np.best.threshold,
                NN.SS = NN.SeSp$table,
                NN.inconclusive = NN.inconclusive,
                NN.best.threshold = NN.best.threshold,
                BN.SS = BN.SeSp$table,
                BN.inconclusive = BN.inconclusive,
                BN.best.threshold = BN.best.threshold
                )
  class(output) <- "TGROC"
  
  if (any(Plot.type != "None")) {
  plot(output,
       Plot = Plot,
       Plot.inc.area = Plot.inc.area,
       Plot.Cl = Plot.Cl,
       Plot.threshold = Plot.threshold,
       cex = cex,
       cex.sub = cex.sub)
   }    
  output
}