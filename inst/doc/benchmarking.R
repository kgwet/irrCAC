## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(irrCAC)

## -----------------------------------------------------------------------------
landis.koch
altman
fleiss

## -----------------------------------------------------------------------------
  ac1 <- gwet.ac1.raw(cac.raw4raters)$est
  data.frame(ac1$coeff.val, ac1$coeff.se)
  landis.koch.bf(ac1$coeff.val, ac1$coeff.se)
  altman.bf(ac1$coeff.val, ac1$coeff.se)
  fleiss.bf(ac1$coeff.val, ac1$coeff.se)

