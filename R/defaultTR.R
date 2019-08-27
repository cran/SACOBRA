#defaultTR.r
#'
#' Default settings for the trust-region part of COBRA.
#'
#' Sets default values for the trust-region part \code{cobra$TRlist} of SACOBRA.  \cr
#' With the call \code{\link{setOpts}(myTR,defaultTR())} it is possible to extend a partial list 
#' \code{myTR} to a list containing all \code{TR}-elements (the missing ones are taken from 
#' \code{defaultTR()}).
#'
#' @return a list with the following elements 
#'      \item{shape}{["cube"] Shape of the trust region can be chosen between cube or a sphere \code{[cube|sphere]}}
#'      \item{radiMin}{[0.01] A value betwwen 0 and 1, minimum fraction of the width of the search space to be used as radius of the trust region}
#'      \item{radiMax}{[0.8] A value between 0 and 1, maximum fraction of the width of the search space to be used as radius of the trust region }
#'      \item{radiInit}{[0.1] Initial radius of trust region}
#'      \item{center}{[cobra$xbest] Center of the trust region can be the current bwst solution or the new solution\code{[xbest|xnew]}}
#' @seealso   \code{\link{setOpts}}, \code{\link{trustRegion}}
#' @export
#'
defaultTR<-function(){
  tr<-list(shape="cube",
           radiMin=0.01,
           radiMax=0.8,
           radiInit=0.1,
           center="xbest")
  return(tr)  
}

