check.data <- function(ref, test){
   if (length(ref) != length(test)) stop('parameters ref and test have unequal length')
   sel = stats::complete.cases(ref, test) # sum(sel)
   ref=ref[sel]
   test=test[sel]
   if (any(levels(as.factor(ref)) != c(0, 1))) {
     stop("Your reference standard must be coded as 0 (absence) and 1 (presence). Check reference categories!")
   }
   t=table(test)
   if (length(names(t)) < 20) warning('The Uncertain Interval method has been developed for continuous data. Your test has less than 20 different values.')
   return(data.frame(ref=ref, test=test))
}

# ref[1000]=NA
# test=round(test)
