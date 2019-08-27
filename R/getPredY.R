# getPredY.R
#
# helper functions for cobraPhaseII.R and updateSaveCobra.R
#
getPredY <- function(xNew,cobra) {
  predy <- interpRBF(xNew, cobra$fitnessSurrogate)      
  if (cobra$PLOG[length(cobra$PLOG)]) {
    predy <- plogReverse(predy,tail(cobra$pShift,1))
  }
  return (predy)
}
getPredY0 <- function(xNew,fitnessSurrogate,cobra) {
  predy <- interpRBF(xNew, fitnessSurrogate)      
  if (cobra$PLOG[length(cobra$PLOG)]) {
    predy <- plogReverse(predy,tail(cobra$pShift,1))
  }
  return (predy)
}
getPredY1 <- function(xNew,fitnessSurrogate,cobra) {
  predy <- interpRBF(xNew, fitnessSurrogate)      
  #if (cobra$PLOG[length(cobra$PLOG)]) {
  #  predy <- plogReverse(predy,tail(cobra$pShift,1))
  #}
  return (predy)
}

