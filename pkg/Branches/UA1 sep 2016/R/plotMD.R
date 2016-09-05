#' Function to plot the mixed densities of distributions of individuals with (1) and without (0) the targeted condition. 
#' @name plotMD
#' @description This plot function shows the two distributions and their overlap in a single graph.
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param n.breaks The number of breaks used to construct the histograms. This number has to be equal or lower than the discernable values in the test.
#' @param subtitle Optional subtitle
#' @details The graph shows the two distributions and their overlap. Many tests of intermediate quality have a considerable overlap. Also, the distibutions as estimated by the \code{density} function, using the gaussian kernel is shown. The intersection is indicated by a vertical line. This graph allows the visual inspection of the two distributions, as well a visual inspection of the approximation of the \code{density}, based on the gaussian kernel. When the density estimation is way off, the standard estimation of the intersection will be incorrect, and another estimation has to be supplied.
#'
#' The function \code{plotMD} can also be used for visual inspection of the Uncertain Area (see examples). Please note that the sensitivity and specificity values > .5 (including the default of .55) allows for some positive bias.
#' @return No Value returned.
#' @export
#'
#' @examples
#' # A test of intermediate quality
#' set.seed(1)
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1.2))
#' plotMD(ref, test)
#' ua = uncertain.area(ref, test) # with warning message!
#' # Add lines to indicate Uncertain Area
#' abline(v=ua[1:2])
#' select=(test <= ua[2] & test >= ua[1])
#' # plot the mixed densities for the Uncertain Area
#' plotMD(ref[select], test[select])
plotMD<-function(ref, test, n.breaks=20, subtitle=''){
  stopifnot(all(ref==0 | ref==1))
  low=min(test)
  high=max(test)
  # p1=mean(ref) # getAnywhere(hist.default); methods(hist)
  n.breaks=min(n.breaks, length(unique(test)))
  breaks=seq(low,high, (high-low)/n.breaks)
  
  hy0=hist(test[ref==0], breaks=breaks, plot=FALSE) # str(hy0)
  hy1=hist(test[ref==1], breaks=breaks, plot=FALSE)
  plot(hy0, freq=FALSE, main='Mixed Densities',
       col='skyblue',border=FALSE, xlab='test',
       xlim=c(min(test), max(test)),
       ylim = c(0, max(max(hy0$density), max(hy1$density))) )
  plot(hy1, freq=FALSE, xlim=c(min(test), max(test)),
       add=T,col=scales::alpha('red',.5),border=FALSE)
  
  d0=density(test[ref==0])
  d1=density(test[ref==1])
  # d0$y = (1-p1)*d0$y
  # d1$y = p1*d1$y
  lines(d0, col="darkblue")
  lines(d1, col="darkred")
  intersection=get.intersection(ref, test) # raw data always used
  abline(v=intersection)
  legend('topright', horiz=F, border=0, box.col='transparent' ,
         legend=c('0','overlap','1'),fill=c('skyblue', 'hotpink4','indianred1'),
         cex=.8)
  mtext(subtitle)
}
