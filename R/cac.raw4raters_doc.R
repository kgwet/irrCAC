#' Rating Data from 4 Raters and 12 Subjects.
#' 
#' This dataset contains data from a reliability experiment where 5 observers scored 15 units on a 4-point numeric scale 
#' based on the values 0, 1, 2 and 3.  
#' 
#' @format This dataset contains ratings obtained from an experiment where 4 raters classified 12 subjects into 5 
#' possible categories labeled as 1, 2, 3, 4, and 5. None of the 4 raters scored all 12 units. Therefore, 
#' some missing ratings in the form of "NA" appear in each of the columns associated with the 4 raters.
#' 
#' Note that only the the 4 last columns are to be used with the functions included in this package.  The first 
#' column only plays a descriptive role and is not used in any calculation.
#'  \describe{
#'     \item{Units}{This variable repsents the unit number.}
#'     \item{Rater1}{All ratings from rater 1}
#'     \item{Rater2}{All ratings from rater 2}
#'     \item{Rater3}{All ratings from rater 3}
#'     \item{Rater4}{All ratings from rater 4}
#'  }   
#'  
#' @source Gwet, K.L. (2014) \emph{Handbook of Inter-Rater Reliability}, 4th Edition, page #120. Advanced Analytics, LLC. 
"cac.raw4raters"