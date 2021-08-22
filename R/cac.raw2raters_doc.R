#' Dataset of raw ratings by 2 raters and 12 subjects.
#' 
#' This dataset contains 12 rows and 2 columns. Each row represents a subject
#' and the 2 columns labeled as "rater1" and "rater2" contain the ratings 
#' produced by the 2 raters.  
#' 
#' @format This dataset contains ratings obtained from an experiment where 2 
#' raters classified 12 subjects into 3 categories A, B, and C. Two of the 12
#' subjects were rated by a single rater. Consequently, this dataset contains
#' 2 missing ratings that are identified with the symbol <NA>.
#' 
#'  \describe{
#'     \item{rater1}{All ratings from rater 1}
#'     \item{rater2}{All ratings from rater 2}
#'  }   
#'  
#' @source Gwet, K.L. (2021) \emph{Handbook of Inter-Rater Reliability: 
#' Volume 1}, 5th Edition. AgreeStat Analytics. 
"cac.raw2raters"