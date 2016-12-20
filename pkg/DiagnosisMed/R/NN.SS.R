#' @rdname SS
#' @import AMORE
#' @export
# Making a smoothed SS object with a Neural Network modeling-------------------
NN.SS <- function(x, t.max = NULL, t.min = NULL, precision = .01, CL = .95,
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
                  n.shows = 1){
  # Fitting the Neural Network
  net <- newff(n.neurons = n.neurons, learning.rate.global =
                 learning.rate.global, momentum.global = momentum.global,
               error.criterium = error.criterium, Stao = Stao, hidden.layer =
                 hidden.layer, output.layer = output.layer, method = method)

  # Training the Sensitivity
  net.Se <- train(net, P = x$table$test.values, T = x$table$Sensitivity,
                  error.criterium = error.criterium, report = report, show.step
                  = show.step, n.shows = n.shows)

  # Training the Specificity
  net.Sp <- train(net, P = x$table$test.values, T = x$table$Specificity,
                  error.criterium = error.criterium, report = report, show.step
                  = show.step, n.shows = n.shows)

  # Defining the values where the NN will be simulated
  if (is.null(t.max)) {
    t.max <- max(x$table$test.values)
  }
  if (is.null(t.min)) {
    t.min <- min(x$table$test.values)
  }
  test.values <- seq(t.min, t.max, precision)
  if (length(test.values) < 200) {
    warning("The number of tests values to simulate the neural network is too
            low, below 200. \n Perhaps changing t.max, t.min or precision will
            make a better a fit.")
  }

  # Simulating the NN
  NN <- data.frame(test.values = test.values)
  NN$Sensitivity <- as.numeric(sim(net.Se$net, NN$test.values))
  NN$Se.inf.cl <- NN$Sensitivity - qnorm(1 - (1 - CL) / 2) * sqrt( ( (
    NN$Sensitivity * (1 - NN$Sensitivity)) / (x$sample.size *
                                                x$sample.prevalence)))
  NN$Se.sup.cl <- NN$Sensitivity + qnorm(1 - (1 - CL) / 2) * sqrt(((
    NN$Sensitivity * (1 - NN$Sensitivity)) /
      (x$sample.size * x$sample.prevalence)))
  NN$Specificity <- as.numeric(sim(net.Sp$net, NN$test.values))
  NN$Sp.inf.cl <- NN$Specificity - qnorm(1 - (1 - CL) / 2) * sqrt(((
    NN$Specificity * (1 - NN$Specificity)) / (x$sample.size * (
      1 - x$sample.prevalence))))
  NN$Sp.sup.cl <- NN$Specificity + qnorm(1 - (1 - CL) / 2) * sqrt(((
    NN$Specificity * (1 - NN$Specificity)) / (x$sample.size * (
      1 - x$sample.prevalence))))
  NN$PLR <- NN$Sensitivity / (1 - NN$Specificity)
  NN$PLR.inf.cl <- exp(log(NN$PLR) - (qnorm(1 - ((1 - CL) / 2), mean = 0,
                   sd = 1)) * sqrt((1 - NN$Sensitivity) / ((x$sample.size *
                   x$sample.prevalence) * NN$Specificity) + (NN$Specificity) /
                   ((x$sample.size * (1 - x$sample.prevalence)) * (1 -
                   NN$Specificity))))
  NN$PLR.sup.cl <- exp(log(NN$PLR) + (qnorm(1 - ((1 - CL) / 2), mean = 0,
                   sd = 1)) * sqrt((1 - NN$Sensitivity) / ((x$sample.size *
                   x$sample.prevalence) * NN$Specificity) + (NN$Specificity) /
                   ((x$sample.size * (1 - x$sample.prevalence)) * (1 -
                   NN$Specificity))))
  #NN$NLR <- (1-NN$Specificity)/NN$Sensitivity
  #NN$NLR.inf.cl <- exp(log(NN$NLR)-(qnorm(1-((1-(1-CL))/2),mean=0,sd=1))*
  #     sqrt((NN$Sensitivity)/((x$sample.size * x$sample.prevalence)*(1-NN$Sensitivity))+(1-NN$Specificity)/
  #     ((x$sample.size * (1-x$sample.prevalence))*(NN$Specificity))))
  #NN$NLR.sup.cl <- exp(log(NN$NLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((NN$Sensitivity)/
  #     ((x$sample.size * x$sample.prevalence)*(1-NN$Sensitivity))+(1-NN$Specificity)/((x$sample.size *
  #     (1-x$sample.prevalence))*(NN$Specificity))))
  NN$TP <- NN$Sensitivity * x$sample.prevalence
  NN$FN <- (1 - NN$Sensitivity) * x$sample.prevalence
  NN$FP <- (1 - NN$Specificity) * (1 - x$sample.prevalence)
  NN$TN <- NN$Specificity * (1 - x$sample.prevalence)
  output <- list(table = NN, sample.size = x$sample.size,
                 sample.prevalence = x$sample.prevalence)
  output
  }
