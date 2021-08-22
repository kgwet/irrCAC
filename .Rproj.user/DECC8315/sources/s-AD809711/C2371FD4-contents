#								  AGREE.COEFF3.RAW.R
#						 		    (August 31, 2021)
#Description: This script file contains a series of R functions for computing various agreement coefficients
#	      for multiple raters (2 or more) when the input data file is in the form of nxr matrix or data frame showing 
#             the actual ratings each rater (column) assigned to each subject (in row). That is n = number of subjects, and r = number of raters.
#             A typical table entry (i,g) represents the rating associated with subject i and rater g. 
#Author: Kilem L. Gwet, Ph.D. (gwet@agreestat.com)
#-----------------------------------------------------------------
#     EXAMPLES OF SIMPLE CALLS OF THE MAIN FUNCTIONS:
# > pa.coeff.raw(YourRatings)       # to obtain the percent agreement coefficient
# > gwet.ac1.raw(YourRatings)       # to obtain gwet's AC1 coefficient
# > fleiss.kappa.raw(YourRatings)   # to obtain fleiss' unweighted generalized kappa coefficient
# > krippen.alpha.raw(YourRatings)  # to obtain krippendorff's unweighted alpha coefficient
# > conger.kappa.raw(YourRatings)   # to obtain conger's unweighted generalized kappa coefficient
# > bp.coeff.raw(YourRatings)       # to obtain Brennan-Prediger unweighted coefficient
#
#===========================================================================================
#pa.coeff.raw: Percent agreement coefficient and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#======================================================================================
#' Percent agreement among multiple raters (2, 3, +) when the input data represent the raw ratings reported for each subject and each rater.  
#' @param ratings An nxr matrix / data frame of ratings where each column represents one rater and each row one subject.
#' @param weights is a mandatory parameter that is either a string variable or a matrix. 
#' The string describes one of the predefined weights and must take one of the values 
#' ("quadratic", "ordinal", "linear", "radical", "ratio", "circular", "bipolar"). 
#' If this parameter is a matrix then it must be a square matri qxq where q is the number 
#' of posssible categories where a subject can be classified. If some of the q possible 
#' categories are not used, then it is strobgly advised to specify the complete list of 
#' possible categories as a vector in parametr categ.labels. Otherwise, the program may not work. 
#' @param categ.labels An optional vector parameter containing the list of all possible ratings. It may be useful 
#' in case some of the possibe ratings are not used by any rater, they will still be used when calculating agreement 
#' coefficients. The default value is NULL. In this case, only categories reported by the raters are used in the
#' calculations.
#' @param conflev An optional parameter representing the confidence level associated with the confidence interval. Its default value is 0.95.
#' @param N An optional parameter representing the population size (if any). It may be use to perform the final population correction to the variance.  Its default value is infinity. 
#' @importFrom stats pt qt
#' @return A data list containing 3 objects: (1) a one-row data frame containing the estimates, (2) the weight matrix 
#' used in the calculations, and (3) the categories used in the analysis.  The data frame of estimates contains the 
#' following variables "coeff.name" (coefficient name), "pa" (the percent agreement), "pe" (percent chance-agreement-always equals 0),
#' "coeff.val" (agreement coefficient = pa), coeff.se (the percent agreement standard error), "conf.int" (the percent agreement confidence interval),
#' "p.value"(the percent agreement p-value), "w.name"(the weights' identification). 
#' @examples 
#' #The dataset "cac.raw4raters" comes with this package. Analyze it as follows:
#' cac.raw4raters
#' pa.coeff.raw(cac.raw4raters) #Percent agreement, precision measures, weights & categories
#' pa.coeff.raw(cac.raw4raters)$est #Yields percent agreement with precision measures
#' pa <- pa.coeff.raw(cac.raw4raters)$est$coeff.val #Yields percent agreement alone.
#' pa
#' pa.coeff.raw(cac.raw4raters, weights = "quadratic") #weighted percent agreement/quadratic weights
#' @export
pa.coeff.raw <- function(ratings,weights="unweighted",categ.labels=NULL,conflev=0.95,N=Inf){ 
  pa<-NA;pe<-NA;coeff.val<-NA;coeff.se<-NA;conf.int<-NA;p.value<-NA;w.name<-NA;weights.mat<-NA
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){
    ratings.mat <- trim(toupper(ratings.mat))
    ratings.mat[ratings.mat==''] <- NA_character_
  }
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # finite population correction 
  
  # creating a vector containing all categories used by the raters
  
  if (is.null(categ.labels)){
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  }else categ <- toupper(categ.labels)
  q <- length(categ)

  # creating the weights matrix
  if (is.character(weights)){
    w.name <- weights
    if (weights=="quadratic") weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal") weights.mat<-ordinal.weights(categ)
    else if (weights=="linear") weights.mat<-linear.weights(categ)
    else if (weights=="radical") weights.mat<-radical.weights(categ)
    else if (weights=="ratio") weights.mat<-ratio.weights(categ)
    else if (weights=="circular") weights.mat<-circular.weights(categ)
    else if (weights=="bipolar") weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else{
    w.name <- "Custom weights"
    weights.mat= as.matrix(weights)
  } 
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    categ.is.k <- (ratings.mat==categ[k])
    agree.mat[,k] <- (replace(categ.is.k,is.na(categ.is.k),FALSE)) %*% rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # calculating percent agreement coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  pe <- 0
  coeff.val <- pa
  
  if (q>=2){    
    # calculating variance, stderr & p-value of gwet's ac1 coefficient
    
    den.ivec <- ri.vec*(ri.vec-1)
    den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
    pa.ivec <- (n/n2more)*(sum.q/den.ivec)
    var.pa<-NA;stderr<-NA;stderr.est<-NA;p.value<-NA;lcb<-NA;ucb<-NA
    if (n>=2){
      var.pa <- ((1-f)/(n*(n-1))) * sum((pa.ivec - pa)^2)
      stderr <- sqrt(var.pa)# pa's standard error
      stderr.est <- round(stderr,5)
      p.value <- 2*(1-pt(abs(pa/stderr),n-1))
      lcb <- pa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
      ucb <- min(1,pa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound 
    }
    conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
    coeff.se <- stderr.est
  }
  obs.count <- dplyr::summarise(as.data.frame(1-is.na(ratings)),across(1:r,sum))
  tot.obs <- sum((1-is.na(ratings)))
  coeff.name <- "Percent Agreement"
  df.out <- data.frame(coeff.name,pa,pe,coeff.val,coeff.se,
                       conf.int,p.value,tot.obs,w.name)
  return(list("est" = df.out,"weights" = weights.mat,
              "categories"=categ,"obs"=obs.count))
}
#===========================================================================================
#gwet.ac1.raw: Gwet's AC1/Ac2 coefficient (Gwet(2008)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. L. (2008). ``Computing inter-rater reliability and its variance in the presence of high
#		agreement." British Journal of Mathematical and Statistical Psychology, 61, 29-48.
#============================================================================================
#' Gwet's AC1/AC2 agreement coefficient among multiple raters (2, 3, +) when the input data represent the raw ratings reported for each subject and each rater.
#' @param ratings An nxr matrix / data frame of ratings where each column represents one rater and each row one subject.
#' @param weights is a mandatory parameter that is either a string variable or a matrix. 
#' The string describes one of the predefined weights and must take one of the values 
#' ("quadratic", "ordinal", "linear", "radical", "ratio", "circular", "bipolar"). 
#' If this parameter is a matrix then it must be a square matri qxq where q is the number 
#' of posssible categories where a subject can be classified. If some of the q possible 
#' categories are not used, then it is strobgly advised to specify the complete list of 
#' possible categories as a vector in parametr categ.labels. Otherwise, the program may not work.
#' @param categ.labels An optional vector parameter containing the list of all possible ratings. It may be useful 
#' in case some of the possibe ratings are not used by any rater, they will still be used when calculating agreement 
#' coefficients. The default value is NULL. In this case, only categories reported by the raters are used in the
#' calculations.
#' @param conflev An optional parameter representing the confidence level associated with the confidence interval. Its default value is 0.95.
#' @param N An optional parameter representing the population size (if any). It may be use to perform the final population correction to the variance.  Its default value is infinity. 
#' @return A data list containing 3 objects: (1) a one-row data frame containing various statistics including the 
#' requested agreement coefficient, (2) the weight matrix used in the calculations if any, and (3) the categories 
#' used in the analysis. These could be categories reported by the raters, or those that were available to the raters
#' whether they used them or not.  The output data frame contains the following variables: "coeff.name" 
#' (coefficient name), "pa" (the percent agreement), "pe" (the percent chance agreement), coeff.val (the agreement coefficient 
#' estimate-AC1 or AC2), "coeff.se" (the standard error), "conf.int" (AC1/AC2 confidence interval), "p.value"
#' (Gwet AC1/AC2 p-value), "w.name"(the weights' identification). 
#' @examples 
#' #The dataset "cac.raw4raters" comes with this package. Analyze it as follows:
#' cac.raw4raters
#' gwet.ac1.raw(cac.raw4raters) #AC1 coefficient, precision measures, weights & categories
#' gwet.ac1.raw(cac.raw4raters)$est #Yields AC1 coefficient with precision measures
#' ac1 <- gwet.ac1.raw(cac.raw4raters)$est$coeff.val #Yields AC1 coefficient alone.
#' ac1
#' gwet.ac1.raw(cac.raw4raters, weights = "quadratic") #AC2 coefficient with quadratic wts
#' @references Gwet, K. L. (2008). ``Computing inter-rater reliability and its variance in the presence of high
#' agreement." \emph{British Journal of Mathematical and Statistical Psychology}, 61, 29-48.
#' @export
gwet.ac1.raw <- function(ratings,weights="unweighted",categ.labels=NULL,conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings)
  if (is.character(ratings.mat)){
    ratings.mat <- trim(toupper(ratings.mat))
    ratings.mat[ratings.mat==''] <- NA_character_
  }
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # finite population correction 

  # creating a vector containing all categories used by the raters
 
  if (is.null(categ.labels)){
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  }else categ <- toupper(categ.labels)
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
      w.name <- weights
      if (weights=="quadratic") weights.mat<-quadratic.weights(categ)
      else if (weights=="ordinal") weights.mat<-ordinal.weights(categ)
      else if (weights=="linear") weights.mat<-linear.weights(categ)
      else if (weights=="radical") weights.mat<-radical.weights(categ)
      else if (weights=="ratio") weights.mat<-ratio.weights(categ)
      else if (weights=="circular") weights.mat<-circular.weights(categ)
      else if (weights=="bipolar") weights.mat<-bipolar.weights(categ)
      else weights.mat<-identity.weights(categ)
  }else{
    w.name <- "Custom Weights"
    weights.mat= as.matrix(weights)
  }
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    categ.is.k <- (ratings.mat==categ[k])
    agree.mat[,k] <- (replace(categ.is.k,is.na(categ.is.k),FALSE)) %*% rep(1,r)    
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  # calculating gwet's ac1 coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  if (q>=2){pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))}else pe=1e-15
  gwet.ac1 <- (pa-pe)/(1-pe)
  gwet.ac1.est <- round(gwet.ac1,5)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  ac1.ivec.x <- ac1.ivec - 2*(1-gwet.ac1) * (pe.ivec-pe)/(1-pe)
  var.ac1<-NA;stderr<-NA;stderr.est<-NA;p.value<-NA;lcb<-NA;ucb<-NA
  if (n>=2){
    var.ac1 <- ((1-f)/(n*(n-1))) * sum((ac1.ivec.x - gwet.ac1)^2)
    stderr <- sqrt(var.ac1)# ac1's standard error
    stderr.est <- round(stderr,5)
    p.value <- 2*(1-pt(abs(gwet.ac1/stderr),n-1))
    lcb <- gwet.ac1 - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,gwet.ac1 + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  if (q==1) coeff.name <- "AC1"
  else{
    if (sum(weights.mat)==q) coeff.name <- "AC1"
    else coeff.name <- "AC2"
  }
  coeff.val <- gwet.ac1.est
  coeff.se <- stderr.est
  obs.count <- dplyr::summarise(as.data.frame(1-is.na(ratings)),across(1:r,sum))
  tot.obs <- sum((1-is.na(ratings)))
  df.out <- data.frame(coeff.name,pa,pe,coeff.val,coeff.se,conf.int,
                       p.value,tot.obs,w.name)
  return(list("est" = df.out,"weights" = weights.mat,
              "categories"=categ,"obs"=obs.count))
}

#=====================================================================================
#fleiss.kappa.raw: This function computes Fleiss' generalized kappa coefficient (see Fleiss(1971)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John Wiley & Sons.
#======================================================================================
#' Fleiss' generalized kappa among multiple raters (2, 3, +) when the input data represent the raw ratings reported for each subject and each rater.
#' @param ratings An nxr matrix / data frame of ratings where each column represents one rater and each row one subject.
#' @param weights is a mandatory parameter that is either a string variable or a matrix. 
#' The string describes one of the predefined weights and must take one of the values 
#' ("quadratic", "ordinal", "linear", "radical", "ratio", "circular", "bipolar"). 
#' If this parameter is a matrix then it must be a square matri qxq where q is the number 
#' of posssible categories where a subject can be classified. If some of the q possible 
#' categories are not used, then it is strobgly advised to specify the complete list of 
#' possible categories as a vector in parametr categ.labels. Otherwise, the program may not work.
#' @param categ.labels An optional vector parameter containing the list of all possible ratings. It may be useful in 
#' case some of the possibe ratings are not used by any rater, they will still be used when calculating agreement 
#' coefficients. The default value is NULL. In this case, only categories reported by the raters are used in the
#' calculations.
#' @param conflev An optional parameter representing the confidence level associated with the confidence interval. Its default value is 0.95.
#' @param N An optional parameter representing the population size (if any). It may be use to perform the final population correction to the variance.  Its default value is infinity. 
#' @return A data list containing 3 objects: (1) a one-row data frame containing various statistics including the 
#' requested agreement coefficient, (2) the weight matrix used in the calculations if any, and (3) the categories 
#' used in the analysis. These could be categories reported by the raters, or those that were available to the raters
#' whether they used them or not.  The output data frame contains the following variables: "coeff.name" 
#' (coefficient name-here it will be "Fleiss' Kappa"), "pa" (the percent agreement), "pe" (the percent chance agreement), coeff.val 
#' (the agreement coefficient estimate-Fleiss' Kappa), "coeff.se" (the standard error), "conf.int" (Fleiss Kappa's confidence 
#' interval), "p.value"(Fleiss Kappa's p-value), "w.name"(the weights' identification).
#' @examples 
#' #The dataset "cac.raw4raters" comes with this package. Analyze it as follows:
#' cac.raw4raters
#' fleiss.kappa.raw(cac.raw4raters) #Fleiss' kappa, precision measures, weights & categories
#' fleiss.kappa.raw(cac.raw4raters)$est #Yields Fleiss' kappa with precision measures
#' fleiss <- fleiss.kappa.raw(cac.raw4raters)$est$coeff.val #Yields Fleiss' kappa alone.
#' fleiss
#' fleiss.kappa.raw(cac.raw4raters, weights = "quadratic") #weighted Fleiss' kappa/quadratic wts
#' @references Fleiss, J. L. (1981). \emph{Statistical Methods for Rates and Proportions}. John Wiley \& Sons.
#' @export
fleiss.kappa.raw <- function(ratings,weights="unweighted",categ.labels=NULL,conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){
    ratings.mat <- trim(toupper(ratings.mat))
    ratings.mat[ratings.mat==''] <- NA_character_
  }
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # finite population correction 

  # creating a vector containing all categories used by the raters
 
  if (is.null(categ.labels)){
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  }else categ <- toupper(categ.labels)
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
    w.name <- weights
    if (weights=="quadratic") weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal") weights.mat<-ordinal.weights(categ)
    else if (weights=="linear") weights.mat<-linear.weights(categ)
    else if (weights=="radical") weights.mat<-radical.weights(categ)
    else if (weights=="ratio")weights.mat<-ratio.weights(categ)
    else if (weights=="circular") weights.mat<-circular.weights(categ)
    else if (weights=="bipolar") weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else{
    w.name <- "Custom Weights"
    weights.mat= as.matrix(weights)
  }
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    categ.is.k <- (ratings.mat==categ[k])
    agree.mat[,k] <- (replace(categ.is.k,is.na(categ.is.k),FALSE)) %*% rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating fleiss's generalized kappa coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  if (q>=2){pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))}else pe=1e-15
  fleiss.kappa <- (pa-pe)/(1-pe)
  fleiss.kappa.est <- round(fleiss.kappa,5)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  kappa.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2
  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.vec
  kappa.ivec.x <- kappa.ivec - 2*(1-fleiss.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.fleiss<-NA;stderr<-NA;stderr.est<-NA;p.value<-NA;lcb<-NA;ucb<-NA
  if (n>=2){
    var.fleiss <- ((1-f)/(n*(n-1))) * sum((kappa.ivec.x - fleiss.kappa)^2)
    stderr <- sqrt(var.fleiss)# kappa's standard error
    stderr.est <- round(stderr,5)
    p.value <- 2*(1-pt(abs(fleiss.kappa/stderr),n-1))
    lcb <- fleiss.kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,fleiss.kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Fleiss' Kappa"
  coeff.val <- fleiss.kappa.est
  coeff.se <- stderr.est
  obs.count <- dplyr::summarise(as.data.frame(1-is.na(ratings)),across(1:r,sum))
  tot.obs <- sum((1-is.na(ratings)))
  df.out <- data.frame(coeff.name,pa,pe,coeff.val,coeff.se,conf.int,
                       p.value,tot.obs,w.name)
  return(list("est" = df.out,"weights" = weights.mat,
              "categories"=categ,"obs"=obs.count))
}
#=====================================================================================
#krippen.alpha.raw: This function computes Krippendorff's alpha coefficient (see Krippendorff(1970, 1980)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The algorithm used to compute krippendorff's alpha is very different from anything that was published on this topic. Instead,
#it follows the equations presented by K. Gwet (2012)
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. (2014). Handbook of Inter-Rater Reliability: The Definitive Guide to Measuring the Extent of Agreement Among 
#	Multiple Raters, 4th Edition.  Advanced Analytics, LLC;
#Krippendorff (1970). "Bivariate agreement coefficients for reliability of data." Sociological Methodology,2,139-150
#Krippendorff (1980). Content analysis: An introduction to its methodology (2nd ed.), New-bury Park, CA: Sage.
#======================================================================================
#' Krippendorff's alpha coefficient for an arbitrary number of raters (2, 3, +) when the input data represent the raw ratings reported for each subject and each rater.
#' @param ratings An nxr matrix / data frame of ratings where each column represents one rater and each row one subject.
#' @param weights is a mandatory parameter that is either a string variable or a matrix. 
#' The string describes one of the predefined weights and must take one of the values 
#' ("quadratic", "ordinal", "linear", "radical", "ratio", "circular", "bipolar"). 
#' If this parameter is a matrix then it must be a square matri qxq where q is the number 
#' of posssible categories where a subject can be classified. If some of the q possible 
#' categories are not used, then it is strobgly advised to specify the complete list of 
#' possible categories as a vector in parametr categ.labels. Otherwise, the program may not work.
#' @param categ.labels An optional vector parameter containing the list of all possible ratings. It may be useful in 
#' case some of the possibe ratings are not used by any rater, they will still be used when calculating agreement 
#' coefficients. The default value is NULL. In this case, only categories reported by the raters are used in the
#' calculations.
#' @param conflev An optional parameter representing the confidence level associated with the confidence interval. Its default value is 0.95.
#' @param N An optional parameter representing the population size (if any). It may be use to perform the final population correction to the variance.  Its default value is infinity. 
#' @return A data list containing 3 objects: (1) a one-row data frame containing various statistics including the 
#' requested agreement coefficient-in this case, Krippendorff's alpha, (2) the weight matrix used in the calculations if any, and 
#' (3) the vector of categories used in the analysis. These could be categories reported by the raters, or those that were available 
#' to the raters whether they used them or not.  The output data frame contains the following variables: "coeff.name" (coefficient 
#' name), "pa" (the percent agreement), "pe" (the percent chance agreement), coeff.val (Krippendorff's alpha estimate), "coeff.se 
#' (standard error), conf.int" (Krippendorff alpha's confidence interval),"p.value" (Krippendorff alpha's p-value), "w.name" 
#' (the weights' identification).
#' @examples 
#' #The dataset "cac.raw4raters" comes with this package. Analyze it as follows:
#' cac.raw4raters
#' krippen.alpha.raw(cac.raw4raters) #Alpha coeff. , precision measures, weights & categories
#' krippen.alpha.raw(cac.raw4raters)$est #Krippendorff's alpha with precision measures
#' alpha <- krippen.alpha.raw(cac.raw4raters)$est$coeff.val #Krippendorff's alpha alone.
#' alpha
#' krippen.alpha.raw(cac.raw4raters, weights = "quadratic") #weighted alpha/ quadratic wts
#' @references Gwet, K. (2014). \emph{Handbook of Inter-Rater Reliability: The Definitive Guide to Measuring the Extent of Agreement Among Multiple Raters}, 4th Edition.  Advanced Analytics, LLC.\cr\cr
#' Krippendorff (1970). ``Bivariate agreement coefficients for reliability of data." \emph{Sociological Methodology},2,139-150.\cr\cr
#' Krippendorff (1980). \emph{Content analysis: An introduction to its methodology} (2nd ed.), New-bury Park, CA: Sage.
#' @export
krippen.alpha.raw <- function(ratings,weights="unweighted",categ.labels=NULL,conflev=0.95,N=Inf){ 
  pa<-NA;pe<-NA;coeff.val<-NA;coeff.se<-NA;conf.int<-NA;p.value<-NA;w.name<-NA;weights.mat<-NA
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){
    ratings.mat <- trim(toupper(ratings.mat))
    ratings.mat[ratings.mat==''] <- NA_character_
  }
  n0 <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n0/N # finite population correction 

  # creating a vector containing all categories used by the raters
 
  if (is.null(categ.labels)){
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  }else categ <- toupper(categ.labels)
  q <- length(categ)

  # creating the weights matrix
  
  if (is.character(weights)){
    w.name <- weights
    if (weights=="quadratic") weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal") weights.mat<-ordinal.weights(categ)
    else if (weights=="linear") weights.mat<-linear.weights(categ)
    else if (weights=="radical") weights.mat<-radical.weights(categ)
    else if (weights=="ratio") weights.mat<-ratio.weights(categ)
    else if (weights=="circular") weights.mat<-circular.weights(categ)
    else if (weights=="bipolar") weights.mat<-bipolar.weights(categ)
    else weights.mat <- identity.weights(categ)
  }else{
    w.name <- "Cutsom Weights"
    weights.mat= as.matrix(weights)
  } 
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n0,ncol=q)
  for(k in 1:q){
    categ.is.k <- (ratings.mat==categ[k])
    agree.mat[,k] <- (replace(categ.is.k,is.na(categ.is.k),FALSE)) %*% rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # calculating krippendorff's alpha coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  agree.mat <- agree.mat[(ri.vec>=2),]
  agree.mat.w <- agree.mat.w[(ri.vec>=2),]
  ri.vec <- ri.vec[(ri.vec>=2)]
  ri.mean <- mean(ri.vec)
  n <- nrow(as.matrix(agree.mat))
  epsi <- 1/sum(ri.vec)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%as.matrix(rep(1,q))
  paprime <-sum(sum.q/(ri.mean*(ri.vec-1)))/n
  pa <- (1-epsi)*paprime + epsi
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/ri.mean))
  if (q>=2){pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))}else  pe=1e-15
  krippen.alpha <- (pa-pe)/(1-pe)
  krippen.alpha.est <- round(krippen.alpha,5)
  krippen.alpha.prime <- (paprime-pe)/(1-pe)
    
  if (q>=2){
    # calculating variance, stderr & p-value of krippendorff's alpha coefficient
    pa.ivec <- sum.q/(ri.mean*(ri.vec-1)) - paprime*(ri.vec-ri.mean)/ri.mean
    krippen.ivec <- (pa.ivec-pe)/(1-pe)
    pi.vec.wk. <- weights.mat%*%pi.vec
    pi.vec.w.k <- t(weights.mat)%*%pi.vec
    pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2
    pe.ivec <- (agree.mat%*%pi.vec.w)/ri.mean - pe * (ri.vec-ri.mean)/ri.mean
    krippen.ivec.x <- krippen.ivec - 2*(1-krippen.alpha.prime) * (pe.ivec-pe)/(1-pe)
    var.krippen<-NA;stderr<-NA;stderr.est<-NA;p.value<-NA;lcb<-NA;ucb<-NA
    if (n>=2){
      var.krippen <- ((1-f)/(n*(n-1))) * sum((krippen.ivec.x - krippen.alpha.prime)^2)
      stderr <- sqrt(var.krippen)# alpha's standard error
      stderr.est <- round(stderr,5)
      p.value <- 2*(1-pt(abs(krippen.alpha/stderr),n0-1))
      lcb <- krippen.alpha - stderr*qt(1-(1-conflev)/2,n0-1) # lower confidence bound
      ucb <- min(1,krippen.alpha + stderr*qt(1-(1-conflev)/2,n0-1)) # upper confidence bound
    }
    conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
    coeff.se <- stderr.est
  }
  coeff.val <- krippen.alpha.est
  obs.count <- dplyr::summarise(as.data.frame(1-is.na(ratings)),
                                across(1:r,sum))
  tot.obs <- sum((1-is.na(ratings)))
  coeff.name <- "Krippendorff's Alpha"
  df.out <- data.frame(coeff.name,pa,pe,coeff.val,coeff.se,
                       conf.int,p.value,tot.obs,w.name)
  return(list("est" = df.out,"weights" = weights.mat,
              "categories"=categ,"obs"=obs.count))
}
#===========================================================================================
#conger.kappa.raw: Conger's kappa coefficient (see Conger(1980)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Conger, A. J. (1980), ``Integration and Generalization of Kappas for Multiple Raters,"
#		Psychological Bulletin, 88, 322-328.
#======================================================================================
#' Conger's generalized kappa coefficient for an arbitrary number of raters (2, 3, +) when the input data represent the raw ratings reported for each subject and each rater.
#' @param ratings An nxr matrix / data frame of ratings where each column represents one rater and each row one subject.
#' @param weights is a mandatory parameter that is either a string variable or a matrix. 
#' The string describes one of the predefined weights and must take one of the values 
#' ("quadratic", "ordinal", "linear", "radical", "ratio", "circular", "bipolar"). 
#' If this parameter is a matrix then it must be a square matri qxq where q is the number 
#' of posssible categories where a subject can be classified. If some of the q possible 
#' categories are not used, then it is strobgly advised to specify the complete list of 
#' possible categories as a vector in parametr categ.labels. Otherwise, the program may not work.
#' @param categ.labels An optional vector parameter containing the list of all possible ratings. It may be useful in 
#' case some of the possibe ratings are not used by any rater, they will still be used when calculating agreement 
#' coefficients. The default value is NULL. In this case, only categories reported by the raters are used in the
#' calculations.
#' @param conflev An optional parameter representing the confidence level associated with the confidence interval. Its default value is 0.95.
#' @param N An optional parameter representing the population size (if any). It may be use to perform the final population correction to the variance.  Its default value is infinity. 
#' @return A data list containing 3 objects: (1) a one-row data frame containing various statistics including the requested agreement
#' coefficient, (2) the weight matrix used in the calculations if any, and (3) A vector of categories used in the analysis. These 
#' could be categories reported by the raters, or those available to the raters whether they used them or not.  The output data frame
#' contains the following variables: "coeff.name" (coefficient name), "pa" (the percent agreement), "pe" (the percent chance 
#' agreement), coeff.val (Conger's Kappa estimate), "coeff.se" (standard error), "conf.int" (Conger Kappa's confidence 
#' interval), "p.value"(agreement coefficient's p-value), "w.name"(the weights' identification).
#' @examples 
#' #The dataset "cac.raw4raters" comes with this package. Analyze it as follows:
#' cac.raw4raters
#' conger.kappa.raw(cac.raw4raters) #Conger's kappa, precision stats, weights & categories
#' conger.kappa.raw(cac.raw4raters)$est #Conger's kappa with precision measures
#' conger <- conger.kappa.raw(cac.raw4raters)$est$coeff.val #Yields Conger's kappa alone.
#' conger
#' conger.kappa.raw(cac.raw4raters, weights = "quadratic") #weighted Conger's kappa
#' @references Conger, A. J. (1980), ``Integration and Generalization of Kappas for Multiple Raters," \emph{Psychological Bulletin}, 88, 322-328.
#' @export
conger.kappa.raw <- function(ratings,weights="unweighted",categ.labels=NULL,conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){
    ratings.mat <- trim(toupper(ratings.mat))
    ratings.mat[ratings.mat==''] <- NA_character_
  }
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # finite population correction 

  # creating a vector containing all categories used by the raters
 
  if (is.null(categ.labels)){
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  }else categ <- toupper(categ.labels)
  q <- length(categ)
  # creating the weights matrix

  if (is.character(weights)){
    w.name <- weights
    if (weights=="quadratic") weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal") weights.mat<-ordinal.weights(categ)
    else if (weights=="linear") weights.mat<-linear.weights(categ)
    else if (weights=="radical") weights.mat<-radical.weights(categ)
    else if (weights=="ratio") weights.mat<-ratio.weights(categ)
    else if (weights=="circular") weights.mat<-circular.weights(categ)
    else if (weights=="bipolar") weights.mat<-bipolar.weights(categ)
    else weights.mat <- identity.weights(categ)
  }else{
    w.name <- "Custom Weights"
    weights.mat <- weights
  } 
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    categ.is.k <- (ratings.mat==categ[k])
    agree.mat[,k] <- (replace(categ.is.k,is.na(categ.is.k),FALSE)) %*% rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  # creating the rxq rater-category matrix representing the distribution of subjects by rater and category

  classif.mat <- matrix(0,nrow=r,ncol=q)
  for(k in 1:q){
  	with.mis <-(t(ratings.mat)==categ[k])
  	without.mis <- replace(with.mis,is.na(with.mis),FALSE)
  	classif.mat[,k] <- without.mis%*%rep(1,n)
  }
  # calculating conger's kappa coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  ng.vec <- classif.mat%*%rep(1,q)
  pgk.mat <- classif.mat/(ng.vec%*%rep(1,q))
  p.mean.k <- (t(pgk.mat)%*%rep(1,r))/r 
  s2kl.mat <- (t(pgk.mat)%*%pgk.mat - r * p.mean.k%*%t(p.mean.k))/(r-1)
  if (q>=2){pe <- sum(weights.mat * (p.mean.k%*%t(p.mean.k) -  s2kl.mat/r))}else  pe=1e-15
  conger.kappa <- (pa-pe)/(1-pe)
  conger.kappa.est <- round(conger.kappa,5)

  # calculating variance, stderr & p-value of conger's kappa coefficient
  
  bkl.mat <- (weights.mat+t(weights.mat))/2
  pe.ivec1 <- r*(agree.mat%*%t(t(p.mean.k)%*%bkl.mat))
  pe.ivec2 = rep(0,n)
  
  lamda.ig.mat=matrix(0,n,r)
  if (is.numeric(ratings.mat)){
  	epsi.ig.mat <-1-is.na(ratings.mat)
  	epsi.ig.mat <- replace(epsi.ig.mat,is.na(epsi.ig.mat),FALSE)
  }else{
  	epsi.ig.mat <- 1-(ratings.mat=="")
  	epsi.ig.mat <- replace(epsi.ig.mat,is.na(epsi.ig.mat),FALSE)
  }
  for(k in 1:q){
  	lamda.ig.kmat=matrix(0,n,r)
  	for(l in 1:q){
  		delta.ig.mat <- (ratings.mat==categ[l])
  		delta.ig.mat <- replace(delta.ig.mat,is.na(delta.ig.mat),FALSE)
  		lamda.ig.kmat <- lamda.ig.kmat + weights.mat[k,l] * (delta.ig.mat - (epsi.ig.mat - rep(1,n)%*%t(ng.vec/n)) * (rep(1,n)%*%t(pgk.mat[,l])))
  	}
  	lamda.ig.kmat = lamda.ig.kmat*(rep(1,n)%*%t(n/ng.vec))
  	lamda.ig.mat = lamda.ig.mat+ lamda.ig.kmat*(r*mean(pgk.mat[,k]) - rep(1,n)%*%t(pgk.mat[,k]))
  }
  pe.ivec <- (lamda.ig.mat%*%rep(1,r)) / (r*(r-1))
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  pe.r2 <- pe*(ri.vec>=2)
  conger.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe) 
  conger.ivec.x <- conger.ivec - 2*(1-conger.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.conger<-NA;stderr<-NA;stderr.est<-NA;p.value<-NA;lcb<-NA;ucb<-NA
  if (n>=2){
    var.conger <- ((1-f)/(n*(n-1))) * sum((conger.ivec.x - conger.kappa)^2)
    stderr <- sqrt(var.conger)# conger's kappa standard error
    stderr.est <- round(stderr,5)
    p.value <- 2*(1-pt(abs(conger.kappa/stderr),n-1))
    lcb <- conger.kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,conger.kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Conger's Kappa"
  coeff.val <- conger.kappa.est
  coeff.se <- stderr.est
  obs.count <- dplyr::summarise(as.data.frame(1-is.na(ratings)),across(1:r,sum))
  tot.obs <- sum((1-is.na(ratings)))
  df.out <- data.frame(coeff.name,pa,pe,coeff.val,coeff.se,conf.int,
                       p.value,tot.obs,w.name)
  return(list("est" = df.out,"weights" = weights.mat,
              "categories"=categ,"obs"=obs.count))
}
#===========================================================================================
#bp.coeff.raw: Brennan-Prediger coefficient (see Brennan & Prediger(1981)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, misuses, and alternatives."
#           Educational and Psychological Measurement, 41, 687-699.
#======================================================================================
#' Brennan \& Prediger's (BP) agreement coefficient for an arbitrary number of raters (2, 3, +) when the input data represent the raw ratings reported for each subject and each rater.
#' @param ratings An nxr matrix / data frame of ratings where each column represents one rater and each row one subject.
#' @param weights is a mandatory parameter that is either a string variable or a matrix. 
#' The string describes one of the predefined weights and must take one of the values 
#' ("quadratic", "ordinal", "linear", "radical", "ratio", "circular", "bipolar"). 
#' If this parameter is a matrix then it must be a square matri qxq where q is the number 
#' of posssible categories where a subject can be classified. If some of the q possible 
#' categories are not used, then it is strobgly advised to specify the complete list of 
#' possible categories as a vector in parametr categ.labels. Otherwise, the program may not work.
#' @param categ.labels An optional vector parameter containing the list of all possible ratings. It may be useful in case some of the
#' possibe ratings are not used by any rater, they will still be used when calculating agreement coefficients. The default value is 
#' NULL. In this case, only categories reported by the raters are used in the
#' calculations.
#' @param conflev An optional parameter representing the confidence level associated with the confidence interval. Its default value 
#' is 0.95.
#' @param N An optional parameter representing the population size (if any). It may be use to perform the final population correction
#' to the variance.  Its default value is infinity. 
#' @return A data list containing 3 objects: (1) a one-row data frame containing various statistics including the requested agreement
#' coefficient, (2) the weight matrix used in the calculations if any, and (3) A vector of categories used in the analysis. These 
#' could be categories reported by the raters, or those available to the raters whether they used them or not.  The output data frame
#' contains the following variables: "coeff.name" (coefficient name), "pa" (the percent agreement), "pe" (the percent chance 
#' agreement), coeff.val (Brennan-Prediger coefficient estimate), "coeff.se" (standard error), "conf.int" (the confidence interval), 
#' "p.value"(Brennan-Prediger coefficient's p-value), "w.name"(the weights' identification).
#' @examples 
#' #The dataset "cac.raw4raters" comes with this package. Analyze it as follows:
#' cac.raw4raters
#' bp.coeff.raw(cac.raw4raters) #BP coefficient, precision measures, weights & categories
#' bp.coeff.raw(cac.raw4raters)$est #Brennan-Prediger coefficient with precision measures
#' bp <- bp.coeff.raw(cac.raw4raters)$est$coeff.val #Yields Brennan-Prediger coefficient alone.
#' bp
#' bp.coeff.raw(cac.raw4raters, weights = "quadratic") #weighted Brennan-Prediger coefficient
#' @references Brennan, R.L., \& Prediger, D. J. (1981). ``Coefficient Kappa: some uses, misuses, and alternatives." \emph{Educational and Psychological Measurement}, 41, 687-699.
#' @export
bp.coeff.raw <- function(ratings,weights="unweighted",categ.labels=NULL,conflev=0.95,N=Inf){ 
  ratings.mat <- as.matrix(ratings) 
  if (is.character(ratings.mat)){
    ratings.mat <- trim(toupper(ratings.mat))
    ratings.mat[ratings.mat==""] <- NA_character_
  }
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # finite population correction 

  # creating a vector containing all categories used by the raters
 
  if (is.null(categ.labels)){
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  }else categ <- toupper(categ.labels)
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
    w.name <- weights
    if (weights=="quadratic") weights.mat<-quadratic.weights(categ)
    else if (weights=="ordinal") weights.mat<-ordinal.weights(categ)
    else if (weights=="linear") weights.mat<-linear.weights(categ)
    else if (weights=="radical") weights.mat<-radical.weights(categ)
    else if (weights=="ratio") weights.mat<-ratio.weights(categ)
    else if (weights=="circular") weights.mat<-circular.weights(categ)
    else if (weights=="bipolar") weights.mat<-bipolar.weights(categ)
    else weights.mat<-identity.weights(categ)
  }else{
    w.name <- "Custom Weights"
    weights.mat= as.matrix(weights)
  } 
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    categ.is.k <- (ratings.mat==categ[k])
    agree.mat[,k] <- (replace(categ.is.k,is.na(categ.is.k),FALSE)) %*% rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating gwet's ac1 coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  if (q>=2){pe <- sum(weights.mat) / (q^2)} else pe=1e-15
  bp.coeff <- (pa-pe)/(1-pe)
  bp.coeff.est <- round(bp.coeff,5)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  pe.r2 <- pe*(ri.vec>=2)
  bp.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  
  var.bp<-NA;stderr<-NA;stderr.est<-NA;p.value<-NA;lcb<-NA;ucb<-NA
  if (n>=2){
    var.bp <- ((1-f)/(n*(n-1))) * sum((bp.ivec - bp.coeff)^2)
    stderr <- sqrt(var.bp)# BP's standard error
    stderr.est <- round(stderr,5)
    p.value <- 2*(1-pt(abs(bp.coeff/stderr),n-1))
    lcb <- bp.coeff - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,bp.coeff + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Brennan-Prediger"
  coeff.val <- bp.coeff.est
  coeff.se <- stderr.est
  obs.count <- dplyr::summarise(as.data.frame(1-is.na(ratings)),across(1:r,sum))
  tot.obs <- sum((1-is.na(ratings)))
  df.out <- data.frame(coeff.name,pa,pe,coeff.val,coeff.se,conf.int,
                       p.value,tot.obs,w.name)
  return(list("est" = df.out,"weights" = weights.mat,
              "categories"=categ,"obs"=obs.count))
}
#
#-----  Additional functions needed to run the main functions. If the main functions must be included in another R script, then
# the user will need to add these additional functions to the new script file.
#
# ==============================================================
# trim(x): This is an r function for trimming leading and trealing blanks
# ==============================================================
#' An r function for trimming leading and trealing blanks
#' @param  x is a string variable.
#' @return A string variable where leading and trealing blanks are trimmed.
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) 
}
# # ==============================================================
# # The following functions generate various weight matrices used 
# # in the weighted or unweighted analyses.
# # ==============================================================
# identity.weights<-function(categ){
#   weights<-diag(length(categ))
#   return (weights)
# }
# #-----------------------------
# quadratic.weights<-function(categ){
#   q<-length(categ)
#   weights <- diag(q)
#   if (is.numeric(categ)) { 
#     categ.vec <- sort(categ)
#   }
#   else {
#     categ.vec<-1:length(categ)
#   }
#   xmin<-min(categ.vec)
#   xmax<-max(categ.vec)
#   for(k in 1:q){
#     for(l in 1:q){
#       weights[k,l] <- 1-(categ.vec[k]-categ.vec[l])^2/(xmax-xmin)^2 
#     }
#   }
#   return (weights)
# }
# linear.weights<-function(categ){
#   q<-length(categ)
#   weights <- diag(q)
#   if (is.numeric(categ)) { 
#     categ.vec <- sort(categ)
#   }
#   else {
#     categ.vec<-1:length(categ)
#   }
#   xmin<-min(categ.vec)
#   xmax<-max(categ.vec)
#   for(k in 1:q){
#     for(l in 1:q){
#       weights[k,l] <- 1-abs(categ.vec[k]-categ.vec[l])/abs(xmax-xmin)
#     }
#   }
#   return (weights)
# }
# #--------------------------------
# radical.weights<-function(categ){
#   q<-length(categ)
#   weights <- diag(q)
#   if (is.numeric(categ)) { 
#     categ.vec <- sort(categ)
#   }
#   else {
#     categ.vec<-1:length(categ)
#   }
#   xmin<-min(categ.vec)
#   xmax<-max(categ.vec)
#   for(k in 1:q){
#     for(l in 1:q){
#       weights[k,l] <- 1-sqrt(abs(categ.vec[k]-categ.vec[l]))/sqrt(abs(xmax-xmin))
#     }
#   }
#   return (weights)
# }
# #--------------------------------
# ratio.weights<-function(categ){
#   q<-length(categ)
#   weights <- diag(q)
#   if (is.numeric(categ)) { 
#     categ.vec <- sort(categ)
#   }
#   else {
#     categ.vec<-1:length(categ)
#   }
#   xmin<-min(categ.vec)
#   xmax<-max(categ.vec)
#   for(k in 1:q){
#     for(l in 1:q){
#       weights[k,l] <- 1-((categ.vec[k]-categ.vec[l])/(categ.vec[k]+categ.vec[l]))^2 / ((xmax-xmin)/(xmax+xmin))^2
#     }
#   }
#   return (weights)
# }
# 
# #--------------------------------
# circular.weights<-function(categ){
#   q<-length(categ)
#   weights <- diag(q)
#   if (is.numeric(categ)) { 
#     categ.vec <- sort(categ)
#   }
#   else {
#     categ.vec<-1:length(categ)
#   }
#   xmin<-min(categ.vec)
#   xmax<-max(categ.vec)
#   U = xmax-xmin+1
#   for(k in 1:q){
#     for(l in 1:q){
#       weights[k,l] <- (sin(pi*(categ.vec[k]-categ.vec[l])/U))^2
#     }
#   }
#   weights <- 1-weights/max(weights)
#   return (weights)
# }
# 
# #--------------------------------
# bipolar.weights<-function(categ){
#   q<-length(categ)
#   weights <- diag(q)
#   if (is.numeric(categ)) { 
#     categ.vec <- sort(categ)
#   }
#   else {
#     categ.vec<-1:length(categ)
#   }
#   xmin<-min(categ.vec)
#   xmax<-max(categ.vec)
#   for(k in 1:q){
#     for(l in 1:q){
#       if (k!=l)
#         weights[k,l] <- (categ.vec[k]-categ.vec[l])^2 / (((categ.vec[k]+categ.vec[l])-2*xmin)*(2*xmax-(categ.vec[k]+categ.vec[l])))
#       else weights[k,l] <- 0
#     }
#   }
#   weights <- 1-weights/max(weights)
#   return (weights)
# }
# #--------------------------------
# ordinal.weights<-function(categ){
#   q<-length(categ)
#   weights <- diag(q)
#   categ.vec<-1:length(categ)
#   for(k in 1:q){
#     for(l in 1:q){
#       nkl <- max(k,l)-min(k,l)+1
#       weights[k,l] <- nkl * (nkl-1)/2
#     }
#   }
#   weights <- 1-weights/max(weights)
#   return (weights)
# }
