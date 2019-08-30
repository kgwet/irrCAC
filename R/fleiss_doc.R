#' Dataset describing Fleiss' Benchmarking Scale
#' 
#' This dataset contains information describing Fleiss' scale for benchmarking chance-corrected agreement 
#' coefficients such as Gwet AC1/AC2, Kappa and many others.
#' 
#' @format Each row of this dataset describes an interval and the interpretation of the magnitude it represents.
#'  \describe{
#'     \item{lb.FL}{The interval lower bound}
#'     \item{ub.FL}{The interval upper bound}
#'     \item{interp.FL}{The interpretation}
#'  }   
#'  
#' @source Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John Wiley & Sons.
"fleiss"