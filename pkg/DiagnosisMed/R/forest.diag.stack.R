#' @rdname display
#' @import graphics
#' @export
forest.diag.stack <- function(x, type = c("SeSp","PLR"),
                              mar1 = c(5.1,1, 4.1,1),
                              mar.se = c(5.1,9, 4.1,1),
                              mar.sp = c(5.1,1, 4.1,9),
                              mar.plr = c(5.1,1, 4.1,9),
                              seg.list = list(col = "blue", lty = 1, lwd = 1, xpd = NA),
                              points.list = list(pch = 18, cex = 2, col = 1, xpd = NA, type = "p"),
                              var.list = list(x = .01, cex = 1, col = "black", font = 2, adj = 0),
                              num.list = list(x = .7, cex = .95, col = "gray30", font = 3, adj = 0, xpd = NA),
                              x.n = .65,
                              x.prev = .85,
                              adj.se = 1,
                              adj.sp = 0,
                              adj.plr = 0,
                              se.xlim = "auto",
                              sp.xlim = "auto",
                              plr.xlim = "auto",
                              se.xlab = "Sensitivity",
                              sp.xlab = "Specificity",
                              plr.xlab = "Positive likelihood ratio",
                              se.pos = .06,
                              sp.pos = .06,
                              plr.pos = .06,
                              se.col.label = "Se [95% CI]",
                              sp.col.label = "Sp [95% CI]",
                              plr.col.label = "PLR [95% CI]",
                              digits = 3,
                              grid = TRUE){
  if (any(class(x) != c("data.frame","diag.stack"))) {
    stop("'x' is not an object from diag.stack function.")
  }
  if (!all(names(x) %in% c("Tests","n","prevalence","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl"))) {
    stop("'x' do not have the expected column names. Is it from diag.subset function?")
  }
  if (type[1] != "SeSp" && type != "PLR") {
    stop("'type' must be either 'SeSp or 'PLR'.")
  }

  se.estimates = unlist(x$Sensitivity)
  se.ll = unlist(x$Se.inf.cl)
  se.ul = unlist(x$Se.sup.cl)
  sp.estimates = unlist(x$Specificity)
  sp.ll = unlist(x$Sp.inf.cl)
  sp.ul = unlist(x$Sp.sup.cl)
  plr.estimates = unlist(x$PLR)
  plr.ll = unlist(x$PLR.inf.cl)
  plr.ul = unlist(x$PLR.sup.cl)
  N <- round(unlist(x$n))
  prev <- unlist(x$prevalence)

  if (type[1] == "SeSp") {
    if (se.xlim[1] == "auto") {
      se.xlim = c(min(se.ll, na.rm = T), max(se.ul, na.rm = T))
    }
    if (sp.xlim[1] == "auto"){
      sp.xlim = c(min(sp.ll, na.rm = T), max(sp.ul, na.rm = T))
    }
  }
  if (type[1] == "PLR") {
    if (plr.xlim[1] == "auto"){
      plr.xlim = c(min(plr.ll, na.rm = T), max(plr.ul, na.rm = T))
      if (plr.xlim[2] == Inf) {
        plr.xlim[2] <- 10
        warning("plr.xlim upper limit was infinte, then it was set to 10.")
      }
    }
  }

  var.labels = x$Tests
  var.pos <- seq(1, nrow(x) * 2, 2)
  ylim = c(1, var.pos[length(var.pos)])

  opar <- par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
  on.exit(par(opar))

  if (type[1] == "SeSp") {
    par(mfrow = c(1,3))
  }
  if (type[1] == "PLR") {
    par(mfrow = c(1,2))
  }
  par(mar = mar1)
  plot(NA, NA, xlim = c(0,1), ylim = ylim , xlab = "", ylab = "",yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", xpd = NA)

  var.list$y <- var.pos
  var.list$labels <- var.labels
  var.list$las <- 1
  var.list$xpd <- NA
  do.call(text, var.list)

  num.list$labels <- N
  num.list$y <- var.pos
  num.list$x <- x.n
  do.call(text, num.list)

  num.list$labels <- sprintf(paste0("%.",digits,"f"), prev)
  num.list$x <- x.prev
  do.call(text, num.list)

  var.list$labels <- c("N","Prev")
  var.list$y <- ylim[2] + 1
  var.list$x <- c(x.n, x.prev)
  do.call(text, var.list)

  if (type[1] == "SeSp") {
    par(mar = mar.se)
    plot(NA, NA, xlim = se.xlim, ylim = ylim , xlab = se.xlab, ylab = "",yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", xpd = NA) ; axis(1)
    if (grid) grid()

    seg.list$x0 <- se.ll
    seg.list$y0 <- var.pos
    seg.list$x1 <- se.ul
    seg.list$y1 <- var.pos
    do.call(segments, seg.list)

    points.list$x <- se.estimates
    points.list$y <- var.pos
    do.call(points, points.list)

    var.list$x <- se.xlim[1] - se.pos
    var.list$y <- var.pos
    var.list$adj <- adj.se
    var.list$labels <- sprintf(paste0("%.",digits,"f [%.",digits,"f;%.",digits,"f]"),se.estimates,se.ll,se.ul)
    do.call(text, var.list)

    var.list$x <- se.xlim[1] - se.pos
    var.list$y <- ylim[2] + 1
    var.list$adj <- adj.se
    var.list$labels <- se.col.label
    do.call(text, var.list)

    par(mar = mar.sp)
    plot(NA, NA, xlim = sp.xlim, ylim = ylim , xlab = sp.xlab, ylab = "",yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", xpd = NA) ; axis(1)
    if (grid == T) grid()

    seg.list$x0 <- sp.ll
    seg.list$y0 <- var.pos
    seg.list$x1 <- sp.ul
    seg.list$y1 <- var.pos
    do.call(segments, seg.list)

    points.list$x <- sp.estimates
    points.list$y <- var.pos
    do.call(points, points.list)

    var.list$x <- sp.xlim[2] + sp.pos
    var.list$y <- var.pos
    var.list$adj <- adj.sp
    var.list$labels <- sprintf(paste0("%.",digits,"f [%.",digits,"f;%.",digits,"f]"),sp.estimates,sp.ll,sp.ul)
    do.call(text, var.list)

    var.list$y <- ylim[2] + 1
    var.list$labels <- sp.col.label
    do.call(text, var.list)
  }
}


# rm(cat.labels,cat.pos,i,ini,nLevels,tab,var.labels,var.pos,ylim)
#forest.stack(all.comb,lab=lab,cex.var=.85) #,cex.var = .8,cex.cat = .75, x.cat = .1)
