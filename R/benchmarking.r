#--------------------------------
#' Computing Altman's Benchmark Scale Membership Probabilities
#' @param coeff A mandatory parameter representing the estimated value of an agreement coefficient.
#' @param se A mandatory parameter representing the agreement coefficient standard error. 
#' @param BenchDF An optional parameter that is a 3-column data frame containing the Altman's benchmark scale 
#'    information. The 3 columns are the interval lower bound, upper bound, and their interpretation. The default value 
#'    is a small file contained in the package and named \emph{altman.RData}, which describes the official Altman's 
#'    scale intervals and their interpretation. 
#' @importFrom stats na.omit pnorm
#' @return A one-column matrix containing the membership probabilities (c.f. \url{http://agreestat.com/research_papers/inter-rater\%20reliability\%20study\%20design1.pdf}) 
#' @export
altman.bf <- function(coeff,se,BenchDF=altman){
  n <- nrow(BenchDF)
  cmprob <-matrix(0,n,1)
  outprob <- matrix("",n,1)  
  truncate.fact <- pnorm((coeff+1)/se)-pnorm((coeff-1)/se)
  for (i in 1:n){
    if (i==1) {
      cmprob[i] = (pnorm((coeff-BenchDF$lb.AL[i])/se)-pnorm((coeff-BenchDF$ub.AL[i])/se))/truncate.fact
    }else{
      cmprob[i] = cmprob[i-1] + (pnorm((coeff-BenchDF$lb.AL[i])/se)-pnorm((coeff-BenchDF$ub.AL[i])/se))/truncate.fact
    }
    outprob[i] = round(cmprob[i],5)
  }
  al.out <- data.frame(altman$interp.AL,outprob)
  colnames(al.out) <- c("Altman","CumProb")
  rownames(al.out) <- paste0("(",altman$lb.AL," to ",altman$ub.AL,")")
  return(al.out)
}
#--------------------------------
#' Computing Landis-Koch Benchmark Scale Membership Probabilities 
#' @param coeff A mandatory parameter representing the estimated value of an agreement coefficient.
#' @param se A mandatory parameter representing the agreement coefficient standard error. 
#' @param BenchDF An optional parameter that is a 3-column data frame containing the Landis \& Koch's benchmark scale 
#'    information. The 3 columns are the interval lower bound, upper bound, and their interpretation. The default value 
#'    is a small file contained in the package and named \emph{landis.koch.RData}, which describes the official 
#'    Landis \& Koch's scale intervals and their interpretation.
#' @return A one-column matrix containing the membership probabilities (c.f. \url{http://agreestat.com/research_papers/inter-rater\%20reliability\%20study\%20design1.pdf}) 
#' @export
landis.koch.bf <- function(coeff,se,BenchDF=landis.koch){
  n <- nrow(BenchDF)
  cmprob <- matrix(0,n,1)
  outprob <- matrix("",n,1)
  truncate.fact <- pnorm((coeff+1)/se)-pnorm((coeff-1)/se)
  for (i in 1:n){
    if (i==1) {
      cmprob[i,1]=(pnorm((coeff-BenchDF$lb.LK[i])/se)-pnorm((coeff-BenchDF$ub.LK[i])/se))/truncate.fact
    }else{
      cmprob[i,1]=cmprob[i-1,1]+(pnorm((coeff-BenchDF$lb.LK[i])/se)-pnorm((coeff-BenchDF$ub.LK[i])/se))/truncate.fact
    }
    outprob[i]=round(cmprob[i],5)
  }
  lk.out <- data.frame(landis.koch$interp.LK,outprob)
  colnames(lk.out) <- c("Landis-Koch","CumProb")
  rownames(lk.out) <- paste0("(",landis.koch$lb.LK," to ",landis.koch$ub.LK,")")
  return(lk.out)
}
#--------------------------------
#' Computing Fleiss Benchmark Scale Membership Probabilities
#' @param coeff A mandatory parameter representing the estimated value of an agreement coefficient.
#' @param se A mandatory parameter representing the agreement coefficient standard error. 
#' @param BenchDF An optional parameter that is a 3-column data frame containing the Fleiss' benchmark scale information. 
#'    The 3 columns are the interval lower bound, upper bound, and their interpretation. The default value is a small 
#'    file contained in the package and named \emph{fleiss.RData}, which describes the fleiss' scale intervales and their
#'    interpretation. 
#' @return A one-column matrix containing the membership probabilities (c.f. \url{http://agreestat.com/research_papers/inter-rater\%20reliability\%20study\%20design1.pdf}) 
#' @export
fleiss.bf <- function(coeff,se,BenchDF=fleiss){
  n <- nrow(BenchDF)
  cmprob <-matrix(0,3,1)
  outprob <- matrix("",3,1)  
  truncate.fact <- pnorm((coeff+1)/se)-pnorm((coeff-1)/se)
  for (i in 1:3){
    if (i==1) {
      cmprob[i] = (pnorm((coeff-BenchDF$lb.FL[i])/se)-pnorm((coeff-BenchDF$ub.FL[i])/se))/truncate.fact
    }else{
      cmprob[i] = cmprob[i-1] + (pnorm((coeff-BenchDF$lb.FL[i])/se)-pnorm((coeff-BenchDF$ub.FL[i])/se))/truncate.fact
    }
    outprob[i]=round(cmprob[i],5)
  }
  fl.out <- data.frame(fleiss$interp.FL,outprob)
  colnames(fl.out) <- c("Fleiss","CumProb")
  rownames(fl.out) <- paste0("(",fleiss$lb.FL," to ",fleiss$ub.FL,")")
  return(fl.out)
}