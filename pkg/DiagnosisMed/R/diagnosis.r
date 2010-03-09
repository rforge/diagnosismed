diagnosis <- function(a,b=NULL,c=NULL,d=NULL,CL=0.95,print=TRUE,plot=FALSE){
  #require(epitools)
  if(is.numeric(a)){
       if(all(length(a)==1 & length(b)==1 & length(c)==1 & length(d)==1 & !is.matrix(a))){
         reference.name <- 'Not informed'
         index.name <- 'Not informed'
         tab<-as.table(cbind(rbind(d,c),rbind(b,a)))
         dimnames(tab)<-list(index.test=c("negative","positive"),reference.standard=c("negative","positive"))
         TN<-d
         FN<-b
         FP<-c
         TP<-a         
       }
       if(all(length(a) > 1 & length(b) > 1 & is.null(c) & is.null(d) & !is.matrix(a))){
          if(any(is.na(a),is.na(b))){stop('There are NAs either in index test or reference standard. Consider removing or inputing!')}
          if(nlevels(as.factor(a))!=2 | nlevels(as.factor(b))!=2){
             stop('It seems there are more levels then 0 and 1.')
          }
          if(!all(levels(as.factor(a))==c(0,1) & levels(as.factor(b))==c(0,1))){
             stop('Either the index test or the reference test is not correctly coded. 0 and 1 were expected!')
          }
          else{reference.name <- deparse(substitute(a))
               index.name <- deparse(substitute(b))
               tab<-table(b,a,dnn=c(deparse(substitute(b)),deparse(substitute(a))))
               TN<-tab[1,1]
               FN<-tab[1,2]
               FP<-tab[2,1]
               TP<-tab[2,2]
          }
       }
       if(any(is.table(a) | is.matrix(a))){
          if(!all(dim(a)==c(2,2))){
             stop('It seems the inputed table is not 2x2. Check your table output.')
          }   
          else{tab<-a
              if(is.null(dimnames(tab))){
                 reference.name <- 'Not informed'
                 index.name <- 'Not informed'
                 dimnames(tab) <- list(index.test = c('negative','positive'),reference.standard = c('negative','positive'))
              }
              else{
                 reference.name <- names(dimnames(tab)[2])
                 index.name <- names(dimnames(tab)[1])
              }
              TN<-tab[1,1]
              FN<-tab[1,2]
              FP<-tab[2,1]
              TP<-tab[2,2]
          }
       }
  }  
  if(all(any(is.factor(a) | is.character(a)) & any(is.factor(b) | is.character(b)) & !is.matrix(a))){
     if(any(is.na(a) | is.na(b))){
        stop('There seem to be NAs either in the reference standard or index test. Consider removing or inputing!')  
     }
     if(nlevels(as.factor(a))!=2 | nlevels(as.factor(b))!=2){
        stop('It seems there are more levels then negative/absence and positive/presence.')
     }          
     if(!all(levels(as.factor(a))==c("negative","positive") & levels(as.factor(b))==c("negative","positive")) &
        !all(levels(as.factor(a))==c("absence","presence") & levels(as.factor(b))==c("absence","presence"))){
          stop('It seems categories are not correctly coded in either the reference or index test.')
     }
     else{reference.name <- deparse(substitute(a))
          index.name <- deparse(substitute(b))
          tab<-table(b,a,dnn=c(deparse(substitute(b)),deparse(substitute(a))))
          TN<-tab[1,1]
          FN<-tab[1,2]
          FP<-tab[2,1]
          TP<-tab[2,2]
     }
  }       

  tabmarg<-addmargins(tab)        
  Conf.limit<-CL
  # sample size
  n<-sum(tab)
  # prevalence
  p<-(TP+FN)/n
  # sensitivity and confidence limits
  Se<-TP/(TP+FN)
  Se.cl<-as.numeric(binom.wilson(TP, TP+FN, conf.level = CL)[4:5])
  # especificity and confidence limits
  Sp<-TN/(FP+TN)
  Sp.cl<-as.numeric(binom.wilson(TN, FP+TN, conf.level = CL)[4:5])
  # positive and negative likelyhood ratios and confidence limits
  PLR<-Se/(1-Sp)
  # LR confidence limists inspired in epi.tests{epiR}
  PLR.inf.cl<-exp(log(PLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-Se)/(
    (TP+FN)*Sp)+(Sp)/((FP+TN)*(1-Sp))))
  PLR.sup.cl<-exp(log(PLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((1-Se)/(
    (TP+FN)*Sp)+(Sp)/((FP+TN)*(1-Sp))))
  NLR<-(1-Se)/Sp
  NLR.inf.cl<-exp(log(NLR)-(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((Se)/((TP+
    FN)*(1-Se))+(1-Sp)/((FP+TN)*(Sp))))
  NLR.sup.cl<-exp(log(NLR)+(qnorm(1-((1-CL)/2),mean=0,sd=1))*sqrt((Se)/((TP+
    FN)*(1-Se))+(1-Sp)/((FP+TN)*(Sp))))
  #accuracy and confidence limits
  accu<-(TP+TN)/n
  accu.cl<-as.numeric(binom.wilson(TP+TN, n, conf.level = CL)[4:5])
  # positive and negative predictive values and confidence limits
  PPV<-TP/(TP+FP)
  PPV.cl<-as.numeric(binom.wilson(TP, TP+FP, conf.level = CL)[4:5])
  NPV<-TN/(TN+FN)
  NPV.cl<-as.numeric(binom.wilson(TN, TN+FN, conf.level = CL)[4:5])
  # diagnostic odds ratio and confidence limits
  OR<-fisher.test(tab,conf.level=CL)
  DOR<-unname(OR$estimate)
  #DOR<-(TP*TN)/(FP*FN)
  DOR.inf.cl<-OR$conf.int[1]
  DOR.sup.cl<-OR$conf.int[2]
  # error rate and error trade
  #ER<-((FN/(FN+TN))*p)+(((FP/(FP+TP))*(TN+FP))
  ER<-(FN+FP)/n
  ER.cl<-as.numeric(binom.wilson(FN+FP, n, conf.level = CL)[4:5])
  ET<-(FN/FP)
  # pre-test and pos-test odds (to do)
  # area under ROC curve
  AUC<-(Se+Sp)/2
  if(plot==TRUE){
    plot(1-Sp,Se,xlim=c(0,1),ylim=c(0,1))
    segments(0,0,1-Sp,Se,col="red")
    segments(1-Sp,Se,1,1,col="red")
    grid()
  }
  #if(plot==FALSE)
  #  {ROC<-roc.from.table(tab, graph = FALSE)}
  # gives same results as AUC<-(Se+Sp)/2
  Youden<-Se+Sp-1
  Youden.inf.cl<-Youden-qnorm(CL/2)*sqrt(((Se * (1 - Se))/(TP+FN) + ((Sp * (1 - Sp))/(TN+FP))))
  Youden.sup.cl<-Youden+qnorm(CL/2)*sqrt(((Se * (1 - Se))/(TP+FN) + ((Sp * (1 - Sp))/(FP+TN))))
#  rm(ROC)
  rm(tab)
  # results evaluations
  reteval <- list(tabmarg=tabmarg,n=n,p=p,Se=Se,Se.cl=Se.cl,Sp=Sp,Sp.cl=Sp.cl,PLR=PLR,
    PLR.inf.cl=PLR.inf.cl,PLR.sup.cl=PLR.sup.cl,NLR=NLR,NLR.inf.cl=NLR.inf.cl,
    NLR.sup.cl=NLR.sup.cl,accu=accu,accu.cl=accu.cl,PPV=PPV,PPV.cl=PPV.cl,NPV=NPV,NPV.cl=NPV.cl,
    DOR=DOR,DOR.inf.cl=DOR.inf.cl,DOR.sup.cl=DOR.sup.cl,ET=ET,ER=ER,ER.cl=ER.cl,
    Youden=Youden,Youden.inf.cl=Youden.inf.cl,Youden.sup.cl=Youden.sup.cl,AUC=AUC,
    Conf.limit=Conf.limit,reference.name=reference.name,index.name=index.name)
  class(reteval) <- "diag"
  if(print==TRUE)  {print(reteval)}
  invisible(reteval)
}