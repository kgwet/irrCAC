#' Dataset describing the Altman's Benchmarking Scale
#' 
#' This dataset contains information describing the Altman scale for benchmarking chance-corrected agreement 
#' coefficients such as Gwet AC1/AC2, Kappa and many others.
#' 
#' @format Each row of this dataset describes an interval and the interpretation of the magnitude it represents.
#'  \describe{
#'     \item{lb.AL}{The interval lower bound}
#'     \item{ub.AL}{The interval upper bound}
#'     \item{interp.AL}{The interpretation}
#'  }   
#'  
#' @source Altman, D.G. (1991). \emph{Practical Statistics for Medical Research}. Chapman and Hall.
#' 
"altman"