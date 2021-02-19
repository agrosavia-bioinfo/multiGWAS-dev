## ---- include = FALSE---------------------------------------------------------
eval <- requireNamespace("updog", quietly = TRUE)
if (eval) {
  eval <- packageVersion("updog") >= "2.0.2"
}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = eval
)

## ----setup--------------------------------------------------------------------
library(updog)
library(ldsep)
library(corrplot)
data("uit", package = "ldsep")
class(uit)

## -----------------------------------------------------------------------------
ploidy <- 4
gp <- format_multidog(x = uit, varname = paste0("Pr_", 0:ploidy))
class(gp)
dim(gp)

## -----------------------------------------------------------------------------
gl <- format_multidog(x = uit, varname = paste0("logL_", 0:ploidy))
gp2 <- gl_to_gp(gl = gl)

## -----------------------------------------------------------------------------
ldout <- ldfast(gp = gp, type = "r2")

## -----------------------------------------------------------------------------
corrplot(corr = ldout$ldmat, 
         method = "color", 
         type = "upper",
         diag = FALSE,
         tl.pos = "n")

## -----------------------------------------------------------------------------
corrplot(corr = ldout$semat, 
         method = "color", 
         type = "upper",
         diag = FALSE,
         tl.pos = "n", 
         is.corr = FALSE)

## -----------------------------------------------------------------------------
graphics::hist(x = ldout$rr, 
               main = "Reliability Ratios", 
               xlab = "Reliability Ratios")

