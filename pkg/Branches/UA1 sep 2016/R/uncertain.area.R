#' Function for the determination of inconclusive test scores
#'
#' @name uncertain.area 
#'
#' @description This function determines an area around the intersection of the two distributions of individuals without (0) and with (1) the targeted condition. The area is restricted both by a maximum sensitivity of the test scores within the uncertain area (max.sens) and by a maximum specificity of the test scores within the uncertain area (max.spec). 
#' @param ref The reference standard. A column in a data frame or a vector indicating the classification by the reference test. The reference standard must be coded either as 0 (absence of the condition) or 1 (presence of the condition)
#' @param test The index test or test under evaluation. A column in a dataset or vector indicating the test results in a continuous scale.
#' @param max.sens (default = .55). Maximum sensitivy of the test scores within the uncertain area. A value below .5 is not allowed, while a value larger than .6 is not recommended.
#' @param max.spec (default = .55). Maximum specificity of the test scores within the uncertain area. A value below .5 is not allowed, while a value larger than .6 is not recommended.
#' @param intersection (default = NULL) When NULL, the intersection is calculated with \code{get.intersection}, which uses the kernel density method to obtain the intersection. When another value is assigned to this parameter, this value is used instead.
#' @param return.first (default = TRUE) Return only the widest possible area, given the restrictions. When FALSE all calculated areas with their sensitivity and specificity are returned. NOTE: This function does not always find a suitable area and can return a vector of NULL values.
#' @details{
#' The Uncertain Area is defined as an area below and above the intersection, with a sensitivity and specificity below a desired value (default .55).
#' 
#' In its core, the \code{uncertain.area} function is non-parametric, but it uses the gaussian kernel for estimating the intersection between the two distributions. Always check whether your results are within reason. If the results are unsatisfactoy, first check on the intersection. The \code{density} function allows for other approximations. Another estimate can be obtained by using a more suitable kernel in the \code{density} function. The parameter \code{intersection} can be used to assign the new estimate to the \code{uncertain.area} method.
#' 
#' Furthermore, only a single intersection is assumed (or an second intersection where the overlap is negligible). If another intersection exists and the overlap around this intersection is considerable, a second uncertain area may be determined by using the parameter \code{intersection}. It should be noted that in most cases, a test with more than one interection with non-neglible overlap is problematic and difficult to apply.
#' 
#' The Uncertain Area method is developed for continuous distributions, although it can be applied to tests with distinguishable categorical distributions. When a test is used with less than 20 discernible values, a warning is issued. The method may work satisfactorily, but results should always be checked carefully.
#' 
#' In general, when estimating decision thresholds, a sample of sufficient size should be used. It is recommended to use at least a sample of 100 patients with the targeted condition, and a 'healthy' sample (without the targeted condition) of the same size or larger. 
#' 
#' The Uncertain Area method is not always capable to delevir results. Clearly, when there is no overlap between the two distributions, there cannot be an uncertain area. A very small area of overlap can also li9mit the pssoibilities to find a solution. When there is no solution found, a vector of NULL values is returned.
#' 
#' Lastly, it should be noted that the Uncertain Area method has been developed recently, and future reserach may provide more satisfactory answers.
#' }
#' @return {A \code{data.frame} of
#'  \describe{
#'  \item{cp.l}{ Lower bound of the Uncertain Area.}
#'  \item{cp.h}{ Upper bound of the Uncertain Area.}
#'  \item{FN}{ Count of false negatives within the Uncertain Area.}
#'  \item{TP}{ Count of true positives within the Uncertain Area.}
#'  \item{TN}{ Count of true negatives within the Uncertain Area.}
#'  \item{FP}{ Count of false positives within the Uncertain Area.}
#'  \item{sensitivity}{ Sensitivity of the test scores within the Uncertain Area.}
#'  \item{specitivity}{ Specitivity of the test scores within the Uncertain Area.}
#' }
#' Only a single row is returned when parameter \code{return.first} = TRUE (default).}
#' @export
#'
#' @examples
#' # A test with considerable overlap, resulting in a relatively large Uncertain Area
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,2))
#' uncertain.area(ref, test)
uncertain.area <-
  function(ref,
           test,
           max.sens = .55,
           max.spec = .55,
           intersection = NULL,
           return.first = T) {
    
    df = check.data(ref, test)
    if (max.sens < .5) stop('Value < .5 invalid for max.sens')
    if (max.spec < .5) stop('Value < .5 invalid for max.spec')
    if (max.sens > .6) warning('Value > .6 not recommended for max.sens')
    if (max.spec > .6) warning('Value > .6 not recommended for max.spec')
    
    # only one relevant intersection assumed!
    # linear tests are used for determination of point of intersection
    # linear test is assumed to have a normal distribution
    if (is.null(intersection)) {
      intersection = get.intersection(df$ref, df$test)
      if (length(intersection) > 1) {
        intersection = utils::tail(intersection, n = 1)
        warning('More than one point of intersection. Highest density used.')
      } else intersection=intersection[1] # other values are ignored
    }
    
    tt <- table(df$test, df$ref) # sort test values; colSums(tt)
    
    pred = sort(unique(df$test))
    
    wi0 = rev(which(pred < intersection)) # count backwards from intersection
    wi1 = which(pred >= intersection)    # count forwards from intersection
    
    o0 = matrix(NA, ncol = 3, nrow = length(wi0))
    colnames(o0) = c('cp.l', 'TN', 'FN')
    temp = cumsum(tt[wi0, 1]) # colSums(tt[wi0[1:47],])
    o0[, 'cp.l'] = pred[wi0]  # as.numeric(names(temp))
    o0[, 'TN'] = temp
    o0[, 'FN'] = cumsum(tt[wi0, 2]) # head(o0); sprintf('%.20f', ua[1]) # include cutpoint
    # sprintf('%.20f', ua[1]); sprintf('%.20f', pred[wi0[47]])
    # ua[1] >=pred[wi0[47]]; which(pred >= ua[1]); which(abs(test-ua[1]) <= 1e-5)
    # pred[289]>=ua[1] ; test[237]>=ua[1];
    # sprintf('%.20f', pred[289]); sprintf('%.20f', test[237]); sprintf('%.20f', ua[1])
    
    o1 = matrix(NA, ncol = 3, nrow = length(wi1))
    colnames(o1) = c('cp.h', 'FP', 'TP')
    temp = cumsum(tt[wi1, 2])
    o1[, 'cp.h'] = pred[wi1] # as.numeric(names(temp))
    o1[, 'TP'] = temp
    o1[, 'FP'] = cumsum(tt[wi1, 1])  #head(o1)
    
    # find rows where TN/(TN+FP) <= max.sens
    res0 = (lapply(o0[, "TN"], function(r) {
      a = which(r <= o1[, "FP"] * max.sens / (1 - max.sens))
      ifelse(length(a) == 0, return(NA), return(a))
    }))
    
    # create matrix from list
    m = t(sapply(res0, '[', 1:max(sapply(res0, length)))) # dim(m)
    # find rows with at least one valid finding
    rcp.l = which(rowSums(is.na(m)) != ncol(m))
    cp.l = o0[rcp.l, 'cp.l'] # candidates cp.l #cp.l=o0[m, 'cp.l'] # res0[[1]]
    
    # find candidates for cp.h
    cp.h = matrix(o1[m[rcp.l, ], 'cp.h'], nrow = length(cp.l)) # length(o1[m,'cp.h']) # length(cp.h)
    #  df=data.frame(cp.l, cp.h)
    df = data.frame(cbind(cp.l, cp.h))
    df.l = reshape2::melt(df, id = 'cp.l')
    df.l = stats::na.omit(df.l)
    if (ncol(df.l) != 3)
      oo1 =
      data.frame(
        'cp.l' = NA,
        'cp.h' = NA,
        'FN' = NA,
        'TP' = NA,
        'TN' = NA,
        'FP' = NA
      )
    else {
      colnames(df.l) <- c('cp.l', 'variable', 'cp.h')
      df.l = df.l[order(df.l$cp.l), c('cp.l', 'cp.h')]
      
      m.cp.l = match(df.l$cp.l, o0[, 'cp.l'])
      m.cp.h = match(df.l$cp.h, o1[, 'cp.h'])
      oo1 = data.frame(
        'cp.l' = df.l$cp.l,
        'cp.h' = df.l$cp.h,
        'FN' = o0[m.cp.l, 'FN'],
        'TP' = o1[m.cp.h, 'TP'],
        'TN' = o0[m.cp.l, 'TN'],
        'FP' = o1[m.cp.h , 'FP']
      )
      oo1$sensitivity = oo1$TP / (oo1$TP + oo1$FN)
      oo1$specificity = oo1$TN / (oo1$FP + oo1$TN) # head(oo1)
      oo1 = oo1[oo1$specificity <= max.spec &
                  oo1$sensitivity <= max.sens & stats::complete.cases(oo1),]
      oo1 = oo1[!duplicated(oo1[, c('FN', 'TP', 'TN', 'FP')]), ] # nrow(o1)
    }
    
    # find rows where TP/(TP+FN) <= max.spec
    res1 = (lapply(o1[, "TP"], function(r) {
      a = which(r <= o0[, "FN"] * max.spec / (1 - max.spec))
      ifelse(length(a) == 0, return(NA), return(a))
    }))
    
    # create matrix from list
    m = t(sapply(res1, '[', 1:max(sapply(res1, length))))
    rcp.h = which(rowSums(is.na(m)) != ncol(m))
    cp.h = o1[rcp.h, 'cp.h'] # candidates cp.h
    
    # find candidates for cp.l
    cp.l = matrix(o0[m[rcp.h, ], 'cp.l'], nrow = length(cp.h))
    df = data.frame(cbind(cp.l, cp.h))
    df.h = reshape2::melt(df, id = 'cp.h')
    df.h = stats::na.omit(df.h)
    if (ncol(df.h) != 3)
      oo2 =
      data.frame(
        'cp.l' = NA,
        'cp.h' = NA,
        'FN' = NA,
        'TP' = NA,
        'TN' = NA,
        'FP' = NA
      )
    else {
      colnames(df.h) <-
        c('cp.h', 'variable', 'cp.l') # error! Error in `colnames<-`(`*tmp*`, value = c("cp.h", "variable", "cp.l")) :
      # 'names' attribute [3] must be the same length as the vector [1]
      #df.h=df.h[df.h$cp.l > .25,]
      df.h = df.h[order(df.h$cp.h), c('cp.l', 'cp.h')]
      # colnames(df.h) <- c('cp.l', 'cp.h')
      
      m.cp.l = match(df.h$cp.l, o0[, 'cp.l'])
      m.cp.h = match(df.h$cp.h, o1[, 'cp.h'])
      oo2 = data.frame(
        'cp.l' = df.h$cp.l,
        'cp.h' = df.h$cp.h,
        'FN' = o0[m.cp.l, 'FN'],
        'TP' = o1[m.cp.h, 'TP'],
        'TN' = o0[m.cp.l, 'TN'],
        'FP' = o1[m.cp.h , 'FP']
      )
      oo2$sensitivity = oo2$TP / (oo2$TP +
                                    oo2$FN)
      oo2$specificity = oo2$TN / (oo2$FP +
                                    oo2$TN) # head(oo2)
      oo2 = oo2[oo2$specificity <= max.spec &
                  oo2$sensitivity <= max.sens & stats::complete.cases(oo2),]
      oo2 = oo2[!duplicated(oo2[, c('FN', 'TP', 'TN', 'FP')]), ] # nrow(o2)
    }
    o3 = rbind(oo1, oo2) # nrow(o3) # oo1=data.frame()
    o3 = o3[stats::complete.cases(o3), ]
    o3 = o3[order(-o3$cp.h, o3$cp.l),]
    o3 = o3[!duplicated(o3[, c('cp.l', 'cp.h')]), ]
    o3 = o3[!duplicated(o3[, c('FN', 'TP', 'TN', 'FP')]), ] # head(o3,10) # o=t(do.call(rbind, o3))
    if (return.first |
        nrow(o3) == 0)
      return(unlist(o3[1, ]))
    else
      return(t(do.call(rbind, o3)))
  }


