#' Distribution of 223 Psychiatric Patients by Type of of Psychiatric Disorder and Diagnosis Method.
#' 
#' This dataset shows the distribution of 223 psychiatric patients by diagnosis category and by the method used to 
#' obtain the diagnosis. The first method named ``Clinical Diagnosis" (also known as ``Facility Diagnosis") is used 
#' in a service facility (e.g. public hospital, or a community unit) and does not rely on a rigorous application of 
#' research criteria. The second method known as ``Research Diagnosis" is based on a strict application of research 
#' criteria. Column 1 contains the diagnosis categories into which patients are classified with Method 1.  The first
#' row on the other hand, shows categories into which patients are classified with Method 2. 
#' 
#' @format This dataset contains a 4x4 squared table. The first column is never used in the calculations and only 
#' contains row names. Only the last 4 columns are used for computing agreement coefficients. 
#'  \describe{
#'     \item{Diagnosis}{Pregnancy Type. This variable is shown here for information only and is never used by any function 
#'                 in the irrCAC package.}
#'     \item{Schizophrenia}{Ectopic Pregnancy}
#'     \item{Bipolar.Disorder}{Abnormal Intrauterine Pregnancy}
#'     \item{Depression}{Normal Intrauterine Pregnancy}
#'     \item{Other}{Normal Intrauterine Pregnancy}
#'  }   
#'  
#' @source Gwet, K.L. (2014). \emph{Handbook of Inter-Rater Reliability}, 4th Edition. Advanced Analytics, LLC.
"cont4x4diagnosis"