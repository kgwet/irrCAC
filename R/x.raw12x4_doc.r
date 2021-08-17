#' This dataset contains raw categorical ratings that 4 raters assigned to
#' 12 subjects.
#' 
#' Dataset of raw ratings assigned to 12 units by 4 raters. Each row is 
#' associated with a unit identifier and each column to a rater. The columns
#' are named Rater1, Rater2, Rater3 and Rater4 and contain the categories into
#' which the unit was assigned. Category values are 1, 2, 3, 4, 5 and some 
#' ratings are missing.
#' 
#' @format A data frame of 12 rows and 5 columns containing integer values 
#' 1, 2, 3, 4 and 5. Missing data points are represented with a dot (".").
#' integers.
#'  \describe{
#'     \item{Units}{Patient's identifier}
#'     \item{Rater1}{Category into which Rater1 classified the unit}
#'     \item{Rater2}{Category into which Rater2 classified the unit}
#'     \item{Rater3}{Category into which Rater3 classified the unit}
#'     \item{Rater4}{Category into which Rater4 classified the unit}
#'  }   
#'  
#' @source K. Gwet, PhD.
#' 
"x.raw12x4"