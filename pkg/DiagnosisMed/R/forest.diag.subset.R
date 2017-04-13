#' @rdname display
#' @import graphics
#' @export
forest.diag.subset <- function (x, type = c("SeSp","PLR"),
        mar1 = c(5.1,1, 4.1,1),
        mar.se = c(5.1,9, 4.1,1),
        mar.sp = c(5.1,1, 4.1,9),
        mar.plr = c(5.1,1, 4.1,9),
        seg.list = list(col = "blue", lty = 1, lwd = 1, xpd = NA),
        points.list = list(pch = 18, cex = 2, col = 1, xpd = NA, type = "p"),
        var.list = list(x = .1, cex = 1, col = "black", font = 2, adj = 0),
        cat.list = list(x = .2, cex = .95, col = "gray30", font = 3, adj = 0, xpd = NA),
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
        grid = TRUE) {
 if (any(class(x) != c("data.frame","diag.subset"))) {
   stop("'x' is not an object from diag.subset function.")
 }
 if (!all(names(x) %in% c("Variables","Levels","n","prevalence","Sensitivity","Se.inf.cl","Se.sup.cl","Specificity","Sp.inf.cl","Sp.sup.cl","PLR","PLR.inf.cl","PLR.sup.cl"))) {
   stop("'x' do not have the expected column names. Is it from diag.subset function?")
 }
 if (type[1] != "SeSp" && type != "PLR") {
   stop("'type' must be either 'SeSp or 'PLR'.")
 }

 main.pos = 1
 main.se = unlist(x[which(x$Variables == "Overall"), 5:7])
 main.sp = unlist(x[which(x$Variables == "Overall"), 8:10])
 main.plr = unlist(x[which(x$Variables == "Overall"), 11:13])
 main.n = unlist(x[which(x$Variables == "Overall"), 3])
 main.prev = unlist(x[which(x$Variables == "Overall"), 4])
 se.estimates = unlist(x$Sensitivity[-1])
 se.ll = unlist(x$Se.inf.cl[-1])
 se.ul = unlist(x$Se.sup.cl[-1])
 sp.estimates = unlist(x$Specificity[-1])
 sp.ll = unlist(x$Sp.inf.cl[-1])
 sp.ul = unlist(x$Sp.sup.cl[-1])
 plr.estimates = unlist(x$PLR[-1])
 plr.ll = unlist(x$PLR.inf.cl[-1])
 plr.ul = unlist(x$PLR.sup.cl[-1])
 N <- unlist(x$n[-1])
 prev <- unlist(x$prevalence[-1])

 if (type[1] == "SeSp") {
   if (se.xlim[1] == "auto") {
     se.xlim = c(min(unlist(c(main.se[2], se.ll)), na.rm = T), max(unlist(c(main.se[3], se.ul)), na.rm = T))
   }
   if (sp.xlim[1] == "auto"){
     sp.xlim = c(min(unlist(c(main.sp[2], sp.ll)), na.rm = T), max(unlist(c(main.sp[3], sp.ul)), na.rm = T))
   }
 }
 if (type[1] == "PLR") {
   if (plr.xlim[1] == "auto"){
     plr.xlim = c(min(unlist(c(main.plr[2], plr.ll)), na.rm = T), max(unlist(c(main.plr[3], plr.ul)), na.rm = T))
     if (plr.xlim[2] == Inf) {
       plr.xlim[2] <- 10
       warning("plr.xlim upper limit was infinte, then it was set to 10.")
     }
   }
 }


 var.labels = x$Variables[-which(is.na(x$Variables))]
 cat.labels = x$Levels[-which(is.na(x$Levels))]
 nLevels <- c(which(!is.na(x$Variables))[-1], (nrow(x) + 1))
 for (i in 1:length(nLevels)) {
   nLevels[i] <- nLevels[i + 1] - nLevels[i]
 }
 nLevels <- nLevels[which(!is.na(nLevels))]
 cat.pos <- seq(1:nLevels[1])
 for (i in 2:length(nLevels)) {
   ini <- cat.pos[length(cat.pos)]
   cat.pos <- c(cat.pos, (ini + 3):(ini + 2 + nLevels[i]))
  }
 cat.pos <- cat.pos + 2
 var.pos <- nLevels[1] + 1
 for (i in 2:length(nLevels)) {
   var.pos <- c(var.pos,var.pos[i - 1] + nLevels[i] + 2)
 }
 var.pos <- var.pos + 2
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

 var.list$y <- c(main.pos, var.pos)
 var.list$labels <- var.labels
 var.list$las <- 1
 var.list$xpd <- NA
 do.call(text, var.list)

 var.list$labels <- c(main.n, sprintf(paste0("%.",digits,"f"), main.prev))
 var.list$y <- main.pos
 var.list$x <- c(x.n, x.prev)
 do.call(text, var.list)

 var.list$labels <- c("N","Prev")
 var.list$y <- ylim[2] + 1
 do.call(text, var.list)

 cat.list$y <- cat.pos
 cat.list$x <- cat.list$x
 cat.list$las <- 1
 cat.list$labels <- cat.labels
 do.call(text, cat.list)

 cat.list$x <- x.n
 cat.list$labels <- N
 do.call(text, cat.list)

 cat.list$x <- x.prev
 cat.list$labels <- sprintf(paste0("%.",digits,"f"), prev)
 do.call(text, cat.list)

 if (type[1] == "SeSp") {
   par(mar = mar.se)
   plot(NA, NA, xlim = se.xlim, ylim = ylim , xlab = se.xlab, ylab = "",yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", xpd = NA) ; axis(1)
   if (grid) grid()

   seg.list$x0 <- main.se[2]
   seg.list$y0 <- main.pos
   seg.list$x1 <- main.se[3]
   seg.list$y1 <- main.pos
   do.call(segments, seg.list)

   points.list$x <- main.se[1]
   points.list$y <- main.pos
   do.call(points, points.list)

   var.list$x <- se.xlim[1] - se.pos
   var.list$y <- main.pos
   var.list$adj <- adj.se
   var.list$labels <- sprintf(paste0("%.",digits,"f [%.",digits,"f;%.",digits,"f]"),main.se[1], main.se[2], main.se[3])
   do.call(text, var.list)

   seg.list$x0 <- se.ll
   seg.list$y0 <- cat.pos
   seg.list$x1 <- se.ul
   seg.list$y1 <- cat.pos
   do.call(segments, seg.list)

   points.list$x <- se.estimates
   points.list$y <- cat.pos
   do.call(points, points.list)

   cat.list$x <- se.xlim[1] - se.pos
   cat.list$y <- cat.pos
   cat.list$adj <- adj.se
   cat.list$labels <- sprintf(paste0("%.",digits,"f [%.",digits,"f;%.",digits,"f]"),se.estimates,se.ll,se.ul)
   do.call(text, cat.list)

   var.list$x <- se.xlim[1] - se.pos
   var.list$y <- ylim[2] + 1
   var.list$adj <- adj.se
   var.list$labels <- se.col.label
   do.call(text, var.list)

   par(mar = mar.sp)
   plot(NA, NA, xlim = sp.xlim, ylim = ylim , xlab = sp.xlab, ylab = "",yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", xpd = NA) ; axis(1)
   if (grid == T) grid()

   seg.list$x0 <- main.sp[2]
   seg.list$y0 <- main.pos
   seg.list$x1 <- main.sp[3]
   seg.list$y1 <- main.pos
   do.call(segments, seg.list)

   points.list$x <- main.sp[1]
   points.list$y <- main.pos
   do.call(points, points.list)

   var.list$x <- sp.xlim[2] + sp.pos
   var.list$y <- main.pos
   var.list$adj <- adj.sp
   var.list$labels <- sprintf(paste0("%.",digits,"f [%.",digits,"f;%.",digits,"f]"),main.sp[1],main.sp[2],main.sp[3])
   do.call(text, var.list)

   seg.list$x0 <- sp.ll
   seg.list$y0 <- cat.pos
   seg.list$x1 <- sp.ul
   seg.list$y1 <- cat.pos
   do.call(segments, seg.list)

   points.list$x <- sp.estimates
   points.list$y <- cat.pos
   do.call(points, points.list)

   cat.list$x <- sp.xlim[2] + sp.pos
   cat.list$y <- cat.pos
   cat.list$adj <- adj.sp
   cat.list$labels <- sprintf(paste0("%.",digits,"f [%.",digits,"f;%.",digits,"f]"), sp.estimates, sp.ll, sp.ul)
   do.call(text, cat.list)

   var.list$x <- sp.xlim[2] + sp.pos
   var.list$y <- ylim[2] + 1
   var.list$adj <- adj.sp
   var.list$labels <- sp.col.label
   do.call(text, var.list)
 }
}
