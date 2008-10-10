plot.diag<-function(x,print=FALSE,...){
  #to do - include an error rate curve
  #consider ohter graphic parameters
  # make the scale and ticks appear only in the superior axis
  #pre-test odd p/(1-p)
  #post-test probability =pto/(1+pto)
  pre.test<-seq(0,1,by=.01)
  #pre.test.odd<-pre.test/(1-pre.test)
  #post.test.odd<-numeric(100)
  #post.test.odd<-pre.test.odd*x$PLR
  post.test<-numeric(100)
  #post.test<-post.test.odd/(1+post.test.odd)
  post.test<-((pre.test/(1-pre.test))*x$PLR)/(1+((pre.test/(1-pre.test))*x$PLR))
  #error.rate<-numeric(100)
  #ER<-((FN/(FN+TN))*p)+(((FP/(FP+TP))*(TN+FP))
  #error.rate<-((1-x$NPV)*pre.test)+(((1-x$PPV)*(x$n-pre.test))
  Result<-cbind(pre.test,post.test)
  if(print==TRUE)
    {print(Result)}
  plot(pre.test,post.test,xlab="Pre-test probability (prevalence)"
    ,ylab="Post-test proabbility (PPV)",type="l")
  #if(error.rate==T)
  #  lines(error.rate,lty=2,col=2)
  grid()
  axis(3)
  legend("bottomright",legend=paste("Positive likelihood ratio: ",formatC(x$PLR, digits=4)))
}
