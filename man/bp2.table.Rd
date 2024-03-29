% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agree.coeff2.r
\name{bp2.table}
\alias{bp2.table}
\title{Brenann-Prediger coefficient for 2 raters}
\usage{
bp2.table(
  ratings,
  weights = identity.weights(1:ncol(ratings)),
  conflev = 0.95,
  N = Inf
)
}
\arguments{
\item{ratings}{A square table of ratings (assume no missing ratings).}

\item{weights}{An optional matrix that contains the weights used in the weighted analysis. By default, this 
parameter contaings the identity weight matrix, which leads to the unweighted analysis.}

\item{conflev}{An optional parameter that specifies the confidence level used for constructing confidence 
intervals. By default the function assumes the standard value of 95\%.}

\item{N}{An optional parameter representing the finite population size if any. It is used to perform the finite
population correction to the standard error. It's default value is infinity.}
}
\value{
A data frame containing the following 5 variables: coeff.name coeff.val coeff.se coeff.ci coeff.pval.
}
\description{
Brenann-Prediger coefficient for 2 raters
}
\examples{
#The dataset "cont3x3abstractors" comes with this package. Analyze it as follows:
bp2.table(cont3x3abstractors) #Yields Brennan-Prediger's coefficient along with precision measures
bp <- bp2.table(cont3x3abstractors)$coeff.val #Yields Brennan-Prediger coefficient alone.
bp
q <- nrow(cont3x3abstractors) #Number of categories
bp2.table(cont3x3abstractors,weights = quadratic.weights(1:q)) #Weighted BP coefficient
}
