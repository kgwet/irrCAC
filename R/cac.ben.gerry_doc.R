#' Ratings of 12 units from 2 raters named Ben and Gerry
#' 
#' This dataset contains ratings that 2 raters named Ben and Gerry assigned to 12 units distributed in 2 groups 
#' "G1" and "G2".
#' @usage cac.ben.gerry[,c(3,4)] or cac.ben.gerry[,c("Ben","Gerry")]
#' @format Each row of this dataset describes an interval and the interpretation of the magnitude it represents.
#'  \describe{
#'     \item{Group}{Group Name}
#'     \item{Units}{Unit number}
#'     \item{Ben}{Ben's Ratings}
#'     \item{Gerry}{Gerry's Ratings}
#'  }   
#' The first 2 columns "Group" and "Units" play a descriptive role here and are not used by any fucntion included in this
#' package. One will typically use \code{cac.ben.gerry[,c(3,4)]} or \code{cac.ben.gerry[,c("Ben","Gerry")]} as input 
#' dataset.
"cac.ben.gerry"