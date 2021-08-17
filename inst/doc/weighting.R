## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(irrCAC)

## -----------------------------------------------------------------------------
  identity.weights(1:3)
  radical.weights(1:3)
  linear.weights(1:3)
  ordinal.weights(1:3)
  quadratic.weights(1:3)
  circular.weights(1:3)
  bipolar.weights(1:3)

## -----------------------------------------------------------------------------
  cont3x3abstractors
  q <- nrow(cont3x3abstractors)
  kappa2.table(cont3x3abstractors,weights = quadratic.weights(1:q))
  scott2.table(cont3x3abstractors,weights = quadratic.weights(1:q))
  gwet.ac1.table(cont3x3abstractors,weights = quadratic.weights(1:q))
  bp2.table(cont3x3abstractors,weights = quadratic.weights(1:q))
  krippen2.table(cont3x3abstractors,weights = quadratic.weights(1:q))
  pa2.table(cont3x3abstractors,weights = quadratic.weights(1:q))

## -----------------------------------------------------------------------------
  pa.coeff.raw(cac.raw4raters,weights = "quadratic")$est
  gwet.ac1.raw(cac.raw4raters,weights = "quadratic")$est
  fleiss.kappa.raw(cac.raw4raters,weights = "quadratic")$est
  krippen.alpha.raw(cac.raw4raters,weights = "quadratic")$est
  conger.kappa.raw(cac.raw4raters,weights = "quadratic")$est
  bp.coeff.raw(cac.raw4raters,weights = "quadratic")$est

## -----------------------------------------------------------------------------
  q <- ncol(distrib.6raters)
  gwet.ac1.dist(distrib.6raters,weights = quadratic.weights(1:q))
  fleiss.kappa.dist(distrib.6raters,weights = quadratic.weights(1:q))
  krippen.alpha.dist(distrib.6raters,weights = quadratic.weights(1:q))
  bp.coeff.dist(distrib.6raters,weights = quadratic.weights(1:q))

