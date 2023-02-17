#						  AGREE.COEFF2.R
#					     (August 25, 2019)
#Description: This script file contains a series of R functions for computing various agreement coefficients
#		  for 2 raters when the input data file is in the form of 2x2 contingency table showing the distribution
#             of subjects by rater, and by category.
#Author: Kilem L. Gwet, Ph.D.
#

#====================================================================================
#kappa2.table: Cohen's kappa (Cohen(1960)) coefficient and its standard error for 2 raters when input dataset is a contingency table
#-------------
#The input data "ratings" is a qxq contingency table showing the distribution of
#subjects by rater, when q is the number of categories.
#======================================================================================
#' Kappa coefficient for 2 raters
#' @param ratings A square or contingency table of ratings (assume no missing ratings). See the 2 datasets 
#' "cont3x3abstractors" and "cont4x4diagnosis" that come with this package as examples.
#' @param weights An optional matrix that contains the weights used in the weighted analysis.
#' @param conflev An optional confidence level for confidence intervals. The default value is the traditional 0.95.
#' @param N An optional population size.  The default value is infinity.
#' @importFrom stats pt qt
#' @return A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
#' @examples 
#' #The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
#' kappa2.table(cont3x3abstractors) #Yields Cohen's kappa along with precision measures
#' kappa <- kappa2.table(cont3x3abstractors)$coeff.val #Yields Cohen's kappa alone.
#' kappa
#' q <- nrow(cont3x3abstractors) #Number of categories
#' kappa2.table(cont3x3abstractors,weights = quadratic.weights(1:q))#weighted kappa/quadratic wts
#' @export
kappa2.table <- function(ratings,weights = identity.weights(1:ncol(ratings)),conflev=0.95,N=Inf){
  ratings <- as.matrix(ratings)
  if(dim(ratings)[1] != dim(ratings)[2]){
  	stop('The contingency table should have the same number of rows and columns!')
  }
  if (ncol(ratings) != ncol(weights)){
    stop('The weight matrix has fewer columns than the contingency table. Please revise your input weights!')
  }
  n <- sum(ratings) # number of subjects
  f <- n/N # finite population correction
  q <- ncol(ratings) # number of categories
  pa <- sum(weights * ratings/n) # percent agreement
  pk. <- (ratings%*%rep(1,q))/n
  p.l <- t((t(rep(1,q))%*%ratings)/n)
  pe <- sum(weights*(pk.%*%t(p.l)))
  kappa <- (pa - pe)/(1 - pe) # weighted kappa
  
  # 2 raters special case variance
  
  pkl <- ratings/n
  pb.k <- weights %*% p.l
  pbl. <- t(weights) %*% pk.
  sum1 <- 0
  for(k in 1:q){
    for(l in 1:q){
      sum1 <- sum1 + pkl[k,l]* (weights[k,l]-(1-kappa)*(pb.k[k] + pbl.[l]))^2
    }
  }
  var.kappa <- ((1-f)/(n*(1-pe)^2)) * (sum1 - (pa-2*(1-kappa)*pe)^2)
  var.kappa <- max(var.kappa,1e-100)
  stderr <- sqrt(var.kappa)# kappa standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(kappa/stderr,n-1)
    lcb <- kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Cohen's Kappa"
  coeff.val <- kappa
  coeff.se <- stderr
  coeff.ci <- conf.int
  coeff.pval <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,coeff.val,coeff.se,coeff.ci,coeff.pval))
}

#scott2.table: Scott's pi coefficient (Scott(1955)) and its standard error for 2 raters when input dataset is a contingency table
#-------------
#The input data "ratings" is a qxq contingency table showing the distribution of
#subjects by rater, when q is the number of categories.
#======================================================================================
#' Scott's coefficient for 2 raters
#' @param ratings A square table of ratings (assume no missing ratings).
#' @param weights An optional matrix that contains the weights used in the weighted analysis. By default, this parameter contaings the identity weight matrix, which leads to the unweighted analysis.
#' @param conflev An optional parameter that specifies the confidence level used for constructing confidence intervals. By default the function assumes the standard value of 95\%.
#' @param N An optional parameter representing the finite population size if any. It is used to perform the finite population correction to the standard error. It's default value is infinity. 
#' @return A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
#' @examples 
#' #The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
#' scott2.table(cont3x3abstractors) #Yields Scott's Pi coefficient along with precision measures
#' scott <- scott2.table(cont3x3abstractors)$coeff.val #Yields Scott's coefficient alone.
#' scott
#' q <- nrow(cont3x3abstractors) #Number of categories
#' scott2.table(cont3x3abstractors,weights = quadratic.weights(1:q)) #weighted Scott's coefficient
#' @export
scott2.table <- function(ratings,weights=identity.weights(1:ncol(ratings)),conflev=0.95,N=Inf){
  ratings <- as.matrix(ratings)
  if(dim(ratings)[1] != dim(ratings)[2]){
  	stop('The contingency table should have the same number of rows and columns!')
  }
  if (ncol(ratings) != ncol(weights)){
    stop('The weight matrix has fewer columns than the contingency table. Please revise your input weights!')
  }
  n <- sum(ratings) # number of subjects
  f <- n/N # finite population correction
  q <- ncol(ratings) # number of categories
  pa <- sum(weights * ratings/n) # percent agreement

  pk. <- (ratings%*%rep(1,q))/n
  p.l <- t((t(rep(1,q))%*%ratings)/n)
  pi.k <- (pk.+p.l)/2
  pe <- sum(weights*(pi.k%*%t(pi.k)))
  scott <- (pa - pe)/(1 - pe) # weighted scott's pi coefficint

	# 2 raters special case variance

  pkl <- ratings/n	     #p_{kl}
  pb.k <- weights %*% p.l    #\ov{p}_{+k}
  pbl. <- t(weights) %*% pk. #\ov{p}_{l+}
  pbk  <- (pb.k + pbl.)/2    #\ov{p}_{k}
  sum1 <- 0
  for(k in 1:q){
  	for(l in 1:q){
  		sum1 <- sum1 + pkl[k,l] * (weights[k,l]-(1-scott)*(pbk[k] + pbk[l]))^2
  	}
  }
  var.scott <- ((1-f)/(n*(1-pe)^2)) * (sum1 - (pa-2*(1-scott)*pe)^2)
  var.scott <- max(var.scott,1e-100)
  stderr <- sqrt(var.scott)# Scott's standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(scott/stderr,n-1)
    lcb <- scott - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,scott + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  } 
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Scott's Pi"
  coeff.val <- scott
  coeff.se <- stderr
  coeff.ci <- conf.int
  coeff.pval <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,coeff.val,coeff.se,coeff.ci,coeff.pval))
}

#gwet.ac1.table: Gwet's AC1/Ac2 coefficient (Gwet(2008)) and its standard error for 2 raters when input dataset is a contingency table
#-------------
#The input data "ratings" is a qxq contingency table showing the distribution of
#subjects by rater, when q is the number of categories.
#======================================================================================
#' Gwet's AC1/AC2 coefficient for 2 raters
#' @param ratings A square table of ratings (assume no missing ratings).
#' @param weights An optional matrix that contains the weights used in the weighted analysis. By default, this 
#' parameter contaings the identity weight matrix, which leads to the unweighted analysis.
#' @param conflev An optional parameter that specifies the confidence level used for constructing confidence 
#' intervals. By default the function assumes the standard value of 95\%.
#' @param N An optional parameter representing the finite population size if any. It is used to perform the finite
#' population correction to the standard error. It's default value is infinity.
#' @return A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
#' @examples 
#' #The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
#' gwet.ac1.table(cont3x3abstractors) #Yields AC1 along with precision measures
#' ac1 <- gwet.ac1.table(cont3x3abstractors)$coeff.val #Yields AC1 coefficient alone.
#' ac1
#' q <- nrow(cont3x3abstractors) #Number of categories
#' gwet.ac1.table(cont3x3abstractors,weights = quadratic.weights(1:q)) #AC2 with quadratic weights
#' @export
gwet.ac1.table <- function(ratings,weights=identity.weights(1:ncol(ratings)),conflev=0.95,N=Inf){
  ratings <- as.matrix(ratings)
  if(dim(ratings)[1] != dim(ratings)[2]){
  	stop('The contingency table should have the same number of rows and columns!')
  }
  if (ncol(ratings) != ncol(weights)){
    stop('The weight matrix has fewer columns than the contingency table. Please revise your input weights!')
  }
  n <- sum(ratings) # number of subjects
  f <- n/N # finite population correction
  q <- ncol(ratings) # number of categories
  pa <- sum(weights * ratings/n) # percent agreement

  pk. <- (ratings%*%rep(1,q))/n
  p.l <- t((t(rep(1,q))%*%ratings)/n)
  pi.k <- (pk.+p.l)/2
  tw <- sum(weights)
  pe <- tw * sum(pi.k *(1-pi.k))/(q*(q-1))
  gwet.ac1 <- (pa - pe)/(1 - pe) # gwet's ac1/ac2 coefficint

	# calculation of variance - standard error - confidence interval - p-value

  pkl <- ratings/n	     #p_{kl}
  sum1 <- 0
  for(k in 1:q){
  	for(l in 1:q){
  		sum1 <- sum1 + pkl[k,l] * (weights[k,l]-2*(1-gwet.ac1)*tw*(1-(pi.k[k] + pi.k[l])/2)/(q*(q-1)))^2
  	}
  }
  var.gwet <- ((1-f)/(n*(1-pe)^2)) * (sum1 - (pa-2*(1-gwet.ac1)*pe)^2)
  var.gwet <- max(var.gwet,1e-100)
  stderr <- sqrt(var.gwet)# ac1's standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(gwet.ac1/stderr,n-1)
    lcb <- gwet.ac1 - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,gwet.ac1 + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  if (sum(weights)==q) coeff.name <- "Gwet's AC1"
  else coeff.name <- "Gwet's AC2"
  coeff.val <- gwet.ac1
  coeff.se <- stderr
  coeff.ci <- conf.int
  coeff.pval <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,coeff.val,coeff.se,coeff.ci,coeff.pval))
}
#bp2.table: Brennan-Prediger coefficient (Brennan & Prediger (1981)) and its standard error for 2 raters when input dataset is a contingency table
#-------------
#The input data "ratings" is a qxq contingency table showing the distribution of
#subjects by rater, when q is the number of categories.
#======================================================================================
#' Brenann-Prediger coefficient for 2 raters
#' @param ratings A square table of ratings (assume no missing ratings).
#' @param weights An optional matrix that contains the weights used in the weighted analysis. By default, this 
#' parameter contaings the identity weight matrix, which leads to the unweighted analysis.
#' @param conflev An optional parameter that specifies the confidence level used for constructing confidence 
#' intervals. By default the function assumes the standard value of 95\%.
#' @param N An optional parameter representing the finite population size if any. It is used to perform the finite
#' population correction to the standard error. It's default value is infinity.
#' @examples 
#' #The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
#' bp2.table(cont3x3abstractors) #Yields Brennan-Prediger's coefficient along with precision measures
#' bp <- bp2.table(cont3x3abstractors)$coeff.val #Yields Brennan-Prediger coefficient alone.
#' bp
#' q <- nrow(cont3x3abstractors) #Number of categories
#' bp2.table(cont3x3abstractors,weights = quadratic.weights(1:q)) #Weighted BP coefficient
#' @return A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
#' @export
bp2.table <- function(ratings,weights=identity.weights(1:ncol(ratings)),conflev=0.95,N=Inf){
  ratings <- as.matrix(ratings)
  if(dim(ratings)[1] != dim(ratings)[2]){
  	stop('The contingency table should have the same number of rows and columns!')
  }
  if (ncol(ratings) != ncol(weights)){
    stop('The weight matrix has fewer columns than the contingency table. Please revise your input weights!')
  }
  n <- sum(ratings) # number of subjects
  f <- n/N # finite population correction
  q <- ncol(ratings) # number of categories
  pa <- sum(weights * ratings/n) # percent agreement

  tw <- sum(weights)
  pe <- tw/(q^2)
  bp.coeff <- (pa - pe)/(1 - pe) # Brennan-Prediger coefficint

	# calculation of variance - standard error - confidence interval - p-value

  pkl <- ratings/n	     #p_{kl}
  sum1 <- 0
  for(k in 1:q){
  	for(l in 1:q){
  		sum1 <- sum1 + pkl[k,l] * weights[k,l]^2
  	}
  }
  var.bp <- ((1-f)/(n*(1-pe)^2)) * (sum1 - pa^2)
  var.bp <- max(var.bp,1e-100)
  stderr <- sqrt(var.bp)# bp's standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(bp.coeff/stderr,n-1)
    lcb <- bp.coeff - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,bp.coeff + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Brennan-Prediger"
  coeff.val <- bp.coeff
  coeff.se <- stderr
  coeff.ci <- conf.int
  coeff.pval <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,coeff.val,coeff.se,coeff.ci,coeff.pval))
}

#krippen2.table: Scott's pi coefficient (Scott(1955)) and its standard error for 2 raters when input dataset is a contingency table
#-------------
#The input data "ratings" is a qxq contingency table showing the distribution of
#subjects by rater, when q is the number of categories.
#======================================================================================
#' Krippendorff's Alpha coefficient for 2 raters
#' @param ratings A square table of ratings (assume no missing ratings).
#' @param weights An optional matrix that contains the weights used in the weighted analysis. By default, this 
#' parameter contaings the identity weight matrix, which leads to the unweighted analysis.
#' @param conflev An optional parameter that specifies the confidence level used for constructing confidence 
#' intervals. By default the function assumes the standard value of 95\%.
#' @param N An optional parameter representing the finite population size if any. It is used to perform the finite
#' population correction to the standard error. It's default value is infinity.
#' @return A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
#' @examples 
#' #The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
#' krippen2.table(cont3x3abstractors) #Krippendorff's alpha along with precision measures
#' alpha <- krippen2.table(cont3x3abstractors)$coeff.val #Krippendorff's alpha alone.
#' alpha
#' q <- nrow(cont3x3abstractors) #Number of categories
#' krippen2.table(cont3x3abstractors,weights = quadratic.weights(1:q)) #Weighted alpha coefficient
#' @export
krippen2.table <- function(ratings,weights=identity.weights(1:ncol(ratings)),conflev=0.95,N=Inf){
  ratings <- as.matrix(ratings)
  if(dim(ratings)[1] != dim(ratings)[2]){
  	stop('The contingency table should have the same number of rows and columns!')
  }
  if (ncol(ratings) != ncol(weights)){
    stop('The weight matrix has fewer columns than the contingency table. Please revise your input weights!')
  }
  n <- sum(ratings) # number of subjects
  f <- n/N # finite population correction
  q <- ncol(ratings) # number of categories
  epsi = 1/(2*n)
  pa0 <- sum(weights * ratings/n)
  pa <- (1-epsi)*pa0 + epsi # percent agreement

  pk. <- (ratings%*%rep(1,q))/n
  p.l <- t((t(rep(1,q))%*%ratings)/n)
  pi.k <- (pk.+p.l)/2
  pe <- sum(weights*(pi.k%*%t(pi.k)))
  kripp.coeff <- (pa - pe)/(1 - pe) # weighted Krippen's alpha coefficint

	# calculating variance

  pkl <- ratings/n	     #p_{kl}
  pb.k <- weights %*% p.l    #\ov{p}_{+k}
  pbl. <- t(weights) %*% pk. #\ov{p}_{l+}
  pbk  <- (pb.k + pbl.)/2    #\ov{p}_{k}
  kcoeff <- (pa0 - pe)/(1 - pe)
  sum1 <- 0
  for(k in 1:q){
	for(l in 1:q){
		sum1 <- sum1 + pkl[k,l] * (weights[k,l]-(1-kcoeff)*(pbk[k] + pbk[l]))^2
	}
  }
  var.kripp <- ((1-f)/(n*(1-pe)^2)) * (sum1 - (pa0-2*(1-kcoeff)*pe)^2)
  var.kripp <- max(var.kripp,1e-100)
  stderr <- sqrt(var.kripp)# Kripp. alpha's standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(kripp.coeff/stderr,n-1)
    lcb <- kripp.coeff - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,kripp.coeff + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Krippendorff's Alpha"
  coeff.val <- kripp.coeff
  coeff.se <- stderr
  coeff.ci <- conf.int
  coeff.pval <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,coeff.val,coeff.se,coeff.ci,coeff.pval))
}



#bangdiwala.fn: Bangdiwala's B coefficient (see Shankar & Bangdiwala, 2014) and its standard error for 2 raters when input dataset is a contingency table
#-------------
#The input data "ratings" is a qxq contingency table showing the distribution of
#subjects by rater, when q is the number of categories.
#======================================================================================
#' Bangdiwala B coefficient for 2 raters
#' @param ratings A square table of ratings (assume no missing ratings).
#' @param weights An optional matrix that contains the weights used in the weighted analysis. By default, this parameter contaings the identity weight matrix, which leads to the unweighted analysis.
#' @param conflev An optional parameter that specifies the confidence level used for constructing confidence intervals. By default the function assumes the standard value of 95\%.
#' @param N An optional parameter representing the finite population size if any. It is used to perform the finite population correction to the standard error. It's default value is infinity. 
#' @return A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
#' @examples 
#' #The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
#' bangdiwala.fn(cont3x3abstractors) #Yields Scott's Pi coefficient along with precision measures
#' Bcoeff <- bangdiwala.fn(cont3x3abstractors)$coeff.val #Yields Scott's coefficient alone.
#' Bcoeff
#' q <- nrow(cont3x3abstractors) #Number of categories
#' @export
bangdiwala.table <- function(ratings,conflev=0.95,N=Inf){
  ratings <- as.matrix(ratings)
  if(dim(ratings)[1] != dim(ratings)[2]){
    stop('The contingency table should have the same number of rows and columns!')
  }
  n <- sum(ratings) # number of subjects
  f <- n/N # finite population correction
  q <- ncol(ratings) # number of categories
  pk. <- rowSums(ratings)/n
  p.k <- colSums(ratings)/n
  b1 <- sum((diag(ratings)/n)**2)
  b2 <- c(pk.%*%p.k)
  bcoeff <- as.vector(b1/b2) # Bangdiwala B coefficient

  # Variance of Bangdiwala's B coefficient
  
  pi.k <- (pk.+p.k)/2
  pkl = ratings/n
  pkk <- diag(ratings)/n
  var1 = 2*sum((pkk^2)*(pkk-2*bcoeff*pi.k))
  var2 = (bcoeff^2)*sum(pi.k*pk.*p.k + (pkl%*%pk.)*p.k)
  var = 2*(var1+var2)/(n*(b2^2))
  varB <- (1-f)*var
  varB <- max(varB,1e-100)
  stderr <- sqrt(varB)  #standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(bcoeff/stderr,n-1)
    lcb <- bcoeff - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,bcoeff + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  } 
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Bangdiwala's B"
  coeff.val <- bcoeff
  coeff.se <- stderr
  coeff.ci <- conf.int
  coeff.pval <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,coeff.val,coeff.se,coeff.ci,coeff.pval))
}


#bangdiwala2RR.fn: Bangdiwala's B coefficient (see Shankar & Bangdiwala, 2014) and 
#                  its standard error for 2 raters when input dataset is a made up of
#                  2 columns of raw data.
#======================================================================================
#' Bangdiwala B coefficient for 2 raters when input dataset is made up of 2 columns of raw data.
#' @param fra.ratings.raw A dataframe with 2 columns of raw ratings.
#' @param conflev An optional parameter that specifies the confidence level used for constructing confidence intervals. By default the function assumes the standard value of 95\%.
#' @param N An optional parameter representing the finite population size if any. It is used to perform the finite population correction to the standard error. It's default value is infinity. 
#' @return A data frame containing the following 9 variables: coeff.name, b1, b2, 
#' coeff.val, coeff.se, conf.int, p.value, n and name of the weight used.
#' @examples 
#' The dataset cac.ben.gerry comes with this package. Analyze it as follows:
#' bangdiwala2RR.fn(cac.ben.gerry[,c(3,4)]), using only the last 2 columns.
#' The result will be following:
#' coeff.name        pa        pe coeff.val  coeff.se  conf.int     p.val tot.obs
#' 1 Bangdiwala''s B 0.1322314 0.2066116      0.64 0.2158518 (0.159,1) 7.083e-03      11
#' w.name
#' 1 Identity
#' @export
bangdiwala2RR.fn <- function(fra.ratings.raw,conflev=0.95,N=Inf){
  tib.ratings <- fra.ratings.raw
  tib.ratings[tib.ratings==''] <- NA_character_
  fr.data <- tib.ratings %>%
    count(tib.ratings[1],tib.ratings[2]) %>%
    drop_na()
  ratings <- long2wide.fn(fr.data)
  n <- nrow(tib.ratings) # number of subjects
  if (n<=N) f <- n/N # finite population correction
  else f=0
  q <- ncol(ratings) # number of categories
  pk. <- rowSums(ratings)/n
  p.k <- colSums(ratings)/n
  b1 <- sum((diag(ratings)/n)**2)
  b2 <- c(pk.%*%p.k)
  bcoeff <- as.vector(b1/b2) # Bangdiwala B coefficient
  
  # Variance of Bangdiwala's B coefficient
  
  pi.k <- (pk.+p.k)/2
  pkl = ratings/n
  pkk <- diag(ratings)/n
  var1 = 2*sum((pkk^2)*(pkk-2*bcoeff*pi.k))
  var2 = (bcoeff^2)*sum(pi.k*pk.*p.k + (pkl%*%pk.)*p.k)
  var = 2*(var1+var2)/(n*(b2^2))
  varB <- (1-f)*var
  varB <- max(varB,1e-100)
  stderr <- sqrt(varB)  #standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(bcoeff/stderr,n-1)
    lcb <- bcoeff - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,bcoeff + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  } 
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Bangdiwala's B"
  pa <- b1
  pe <- b2
  coeff.val <- bcoeff
  coeff.se <- stderr
  tot.obs <- n
  w.name <- "Identity"
  #p.value <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,pa,pe,coeff.val,coeff.se,conf.int,p.value,tot.obs,w.name))
}


#pa2.table: Percent agreement coefficient and its standard error for 2 raters when input dataset is a contingency table
#-------------
#The input data "ratings" is a qxq contingency table showing the distribution of
#subjects by rater, when q is the number of categories.
#======================================================================================
#' Percent Agreement coefficient for 2 raters
#' @param ratings A square table of ratings (assume no missing ratings).
#' @param weights An optional matrix that contains the weights used in the weighted analysis. By default, this 
#' parameter contains the identity weight matrix, which leads to the unweighted analysis.
#' @param conflev An optional parameter that specifies the confidence level used for constructing confidence 
#' intervals. By default the function assumes the standard value of 95\%.
#' @param N An optional parameter representing the finite population size if any. It is used to perform the finite
#' population correction to the standard error. It's default value is infinity.
#' @return A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
#' @examples 
#' #The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
#' pa2.table(cont3x3abstractors) #Yields percent agreement along with precision measures
#' pa <- pa2.table(cont3x3abstractors)$coeff.val #Yields percent agreement alone.
#' pa
#' q <- nrow(cont3x3abstractors) #Number of categories
#' pa2.table(cont3x3abstractors,weights = quadratic.weights(1:q)) #Weighted percent agreement
#' @export
pa2.table <- function(ratings,weights=identity.weights(1:ncol(ratings)),
                      conflev=0.95,N=Inf){
  ratings <- as.matrix(ratings)
  if(dim(ratings)[1] != dim(ratings)[2]){
    stop('The contingency table should have the same number of rows and columns!')
  }
  if (ncol(ratings) != ncol(weights)){
    stop('The weight matrix has fewer columns than the contingency table. Please revise your input weights!')
  }
  n <- sum(ratings) # number of subjects
  f <- n/N # finite population correction
  q <- ncol(ratings) # number of categories
  pa <- sum(weights * ratings/n) # percent agreement

  # calculation of variance - standard error - confidence interval - p-value

  pkl <- ratings/n	     #p_{kl}
  sum1 <- 0
  for(k in 1:q){
    for(l in 1:q){
      sum1 <- sum1 + pkl[k,l] * weights[k,l]^2
    }
  }
  var.pa <- ((1-f)/n) * (sum1 - pa^2)
  var.pa <- max(var.pa,1e-100)
  stderr <- sqrt(var.pa) # pa standard error
  p.value <- NA;lcb <- NA;ucb <- NA
  if (n>=2){
    p.value <- 1-pt(pa/stderr,n-1)
    lcb <- pa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
    ucb <- min(1,pa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  }
  conf.int <- paste0("(",round(lcb,3),",",round(ucb,3),")")
  coeff.name <- "Percent Agreement"
  coeff.val <- pa
  coeff.se <- stderr
  coeff.ci <- conf.int
  coeff.pval <- format(p.value,digits=4,nsmall=3,scientific = TRUE)
  return(data.frame(coeff.name,coeff.val,coeff.se,coeff.ci,coeff.pval))
}
#======================================================================================
#' long2wide.fn: This function transforms a 3-column dataset of frequencies to a square
#' matrix or a contingency table. This function uses the freq.supp.fn() function.
#' @param freqs.long A 3-column data frame, where the first 2 variables represent the categories
#' that both raters have actually used when classifying the subjects. The third and last 
#' variable is generally named "n" and represents the count of subjects that classified into
#' the 2 associated categories by both raters.
#' @return A matrix that represents a contingency showing the distribution of subjects by
#' rater and category.
#' @examples 
#' #The dataset "freqs.data" comes with this package. Analyze it as follows:
#' long2wide.fn(freqs.data) #Yields a 5x5 matrix
#' This will produce the following 5x5 matrix:
#' > long2wide.fn(freqs.data)
#' a b c d e
#' a 0 1 0 1 0
#' b 0 2 0 0 0
#' c 0 0 3 0 0
#' d 0 1 0 1 0
#' e 0 0 0 0 1
#' @export
long2wide.fn <- function(freqs.long){
  if (is.character(as.matrix(freqs.long))) freqs.long[freqs.long==""] <- NA
  freq.tab <- freqs.long %>%
    drop_na()
  cat.vec <- as.vector(unique(na.omit(c(freqs.long[[1]],freqs.long[[2]]))))
  freq.tab1 <-freq.supp.fn(freq.tab,cat.vec)
  freq.tab2 <- as_tibble(freq.tab1[order(freq.tab1[[1]],freq.tab1[[2]]),])
  freq.tab2 <- mutate(freq.tab2,n=as.integer(freq.tab2[[3]]))
  dfra.mat <- pivot_wider(data=freq.tab2,names_from = colnames(freq.tab2)[2],
                          values_from = n,values_fill = 0)
  sfra.mat <- as.matrix(select(dfra.mat,all_of(cat.vec)))
  rownames(sfra.mat) <- colnames(sfra.mat)
  return(sfra.mat)
}
#=========================================================================================
#' freq.supp.fn: This function reads a 3-variable input data file containing unique pairs
#' of categories along with their frequency of occurrences, and outputs a similar file
#' where all possible pairs of categories are represented, some with a frequency of 
#' occurrence of 0. 
#' @param freq.data The input data file containing all unique combinations of reported 
#' categories.
#' @param categories.vec A vector containing the complete set of categories available to 
#' raters (e.g. "a", "b", "c", "d", "e"). The raters will not necessarily use all of 
#' these categories.
#' @return This function returns a complete data frame containing all possible combinations of
#' of categories in the categories.vec vector. Newly-added combinations of categories will
#' have a frequency occurrence of 0.
#' @examples 
#' #The dataset "freqs.data" comes with this package. Analyze it as follows:
#' freq.supp.fn(freqs.data) 
#' Executing this command will yield the following data frame:
#'   Ben   Gerry n
#' <chr> <chr> <chr>
#'   a     b     1    
#'   a     d     1
#'   b     b     2    
#'   c     c     3    
#'   d     b     1    
#'   d     d     1    
#'   e     e     1    
#'   a     a     0
#' @export
freq.supp.fn <- function(freq.data,categories.vec){
  cnames <- colnames(freq.data)
  v1.supp <- NULL
  v2.supp <- NULL
  freq.data.all <- freq.data
  r1categ <- sort(as.character(unique(na.omit(freq.data[[1]]))))
  r2categ <- sort(as.character(unique(na.omit(freq.data[[2]]))))
  categories.vec <- sort(as.character(categories.vec))
  if (!identical(categories.vec,r1categ)){
    v1.supp <- setdiff(categories.vec,r1categ)
    n.v2supp<-length(v1.supp)
    freq.v1 <- cbind(v1.supp,r2categ[1:n.v2supp],rep(0,n.v2supp))
    freq.v1 <- setNames(as_tibble(freq.v1),cnames)
    freq.v1[[3]] <- as.integer(as.character(freq.v1[[3]]))
  }
  if (!identical(categories.vec,r2categ)){
    v2.supp <- setdiff(categories.vec,r2categ)
    n.v3supp <- length(v2.supp)
    freq.v2 <- cbind(r1categ[1:n.v3supp],v2.supp,rep(0,n.v3supp))
    freq.v2 <- setNames(as_tibble(freq.v2),cnames)
    freq.v2[[3]] <- as.integer(as.character(freq.v2[[3]]))
  }
  if (length(v1.supp)>0) freq.data.all <- rbind(freq.data.all,freq.v1)
  if (length(v2.supp)>0) freq.data.all <- rbind(freq.data.all,freq.v2)
  return(freq.data.all)
}