## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(irrCAC)

## -----------------------------------------------------------------------------
  data(package="irrCAC")

## -----------------------------------------------------------------------------
  cont3x3abstractors
  kappa2.table(cont3x3abstractors)
  scott2.table(cont3x3abstractors)
  gwet.ac1.table(cont3x3abstractors)
  bp2.table(cont3x3abstractors)
  krippen2.table(cont3x3abstractors)
  pa2.table(cont3x3abstractors)

## -----------------------------------------------------------------------------
  ac1 <- gwet.ac1.table(cont3x3abstractors)$coeff.val

## -----------------------------------------------------------------------------
distrib.6raters
gwet.ac1.dist(distrib.6raters)
fleiss.kappa.dist(distrib.6raters)
krippen.alpha.dist(distrib.6raters)
bp.coeff.dist(distrib.6raters)

## -----------------------------------------------------------------------------
  alpha <- krippen.alpha.dist(distrib.6raters)$coeff

## -----------------------------------------------------------------------------
  ac1 <- gwet.ac1.dist(cac.dist4cat[,2:4])$coeff

## -----------------------------------------------------------------------------
  cac.raw4raters

## -----------------------------------------------------------------------------
pa.coeff.raw(cac.raw4raters)
gwet.ac1.raw(cac.raw4raters)
fleiss.kappa.raw(cac.raw4raters)
krippen.alpha.raw(cac.raw4raters)
conger.kappa.raw(cac.raw4raters)
bp.coeff.raw(cac.raw4raters)

## -----------------------------------------------------------------------------
ac1 <- gwet.ac1.raw(cac.raw4raters)$est
ac1

## -----------------------------------------------------------------------------
ac1 <- gwet.ac1.raw(cac.raw4raters)$est
ac1$coeff.val

