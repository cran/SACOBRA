#
#Samineh Bagheri, Wolfgang Konen
#Cologne University of Applied Sciences
#
#Oct,2015 - Aug,2019
#modifyEquCons.R

#defaultEquMu
#'
#'Default settings for equality handling mechanism
#'
#'Sets suitable defaults for the equality handling part of SACOBRA. \cr
#' \cr The EH technique transforms each equality constraint \eqn{h(\vec{x})=0} into two inequality  
#'        constraints \eqn{h(\vec{x})-\mu <0} and \eqn{-h(\vec{x})-\mu<0 } with an adaptively  
#'        decaying margin \eqn{\mu}.
#'\cr \cr If \code{refine} parameter is set to TRUE, then a refine mechanism is applied to shift the best found solution within the equality margin \eqn{\mu} toward the feasible subspace by minimizing the sum of squared constraint surrogates with a conjugate gradient method.        
#'\deqn{\mbox{Minimize}\quad \sum_i( max(0,g_i(x))^2 )  + \sum_j ( h_j(x)^2 )}
#'
#'
#Details
#
#' With the call \code{\link{setOpts}(equHandle,defaultEquMu())} it is possible to extend a partial list 
#' \code{equHandle} list which is set by user to a list containing all \code{equHandle}-elements (the missing ones are taken from 
#' \code{defaultEquMu()}). 
#'These settings are used by \code{\link{cobraInit}} for initializing the equality margin \eqn{\mu} and by the internal functions \code{\link{updateCobraEqu}} and \code{\link{modifyMu}}. 
#'The minimization step of refine mechanism is done by \code{L-BFGS-B} method in \code{optim} function from \code{stats} package.
#'
#' @return equHandle, a list with the following elements: 
#'    \item{active}{  [TRUE] if set to TRUE the equality-handling (EH) technique is activated. The EH 
#'        technique transforms each equality constraint \eqn{h(\vec{x})=0} into two inequality  
#'        constraints \eqn{h(\vec{x})-\mu <0} and \eqn{-h(\vec{x})-\mu<0 } with an adaptively  
#'        decaying margin \eqn{\mu}.  }
#'    \item{equEpsFinal}{[1e-07] lower bound for margin \eqn{\mu}. \code{equEpsFinal} should be set
#'        to a small but non-zero value (larger than machine accuracy).  }
#'    \item{initType}{["TAV"] the equality margin \eqn{\mu} can be initialized with one of these approaches:\cr ["TAV"|"TMV"|"EMV"|"useGrange"]\cr
#'    \strong{TAV}:  (Total Absolute Violation) takes the median of the sum of violations of the initial population.
#'    \cr \strong{TMV}: (Total Maximum Violation) takes the median of the maximum violation of the initial population
#'    \cr \strong{EMV}: takes the median of the maximum violation of equality constraints of the initial population
#'    \cr \strong{useGrange}: takes the average of the ranges of the equality constraint functions}
#'    \item{epsType}{["SAexpFunc"] type of the function used to modify margin \eqn{\mu} during the optimization process can be one of ["SAexpFunc"|"expFunc"|"Zhang"|"CONS"]. see \code{\link{modifyMu}}.}
#'    \item{dec}{[1.5] decay factor for margin \eqn{\mu}. see \code{\link{modifyMu}}}
#'    \item{refine}{[TRUE] enables the refine mechanism f the equality handling mechanism.}
#'    \item{refineMaxit}{maximum number of iterations used in the refine step. Note that the refine step runs on the surrogate models and does not impose any extra real function evaluation.}
#' @seealso \code{\link{updateCobraEqu}}, \code{\link{modifyMu}}
#' @export
#'
defaultEquMu<-function(){
  equHandle<-list(  active=TRUE,
                    equEpsFinal=1e-07,
                    initType="TAV", #"useGrange", "TAV", "TMV", "EMV"
                    epsType="expFunc", 
                    dec=1.5,
                    refine=TRUE,
                    refineMaxit=100
                  )
  return(equHandle)
}


#updateCobraEqu
#'
#' Calculate cobra$xbest,fbest,ibest given the currently evaluated points.
#' 
#' First, we calculate the set S of points being feasible in the inequalities and feasible
#' with margin cobra$currentEps in the equalities. If S is non-empty, take the point from S 
#' with best objective as cobra$xbest. If S is empty (there are no feasible points): If the new
#' point has lower maxViol than the maxViol of cobra$xbest, take it as the new cobra$xbest. 
#' 
#' Called from updateSaveCobra.
#' 
#' @param cobra     an object of class COBRA (a long list)
#' @param xNew      the last evaluated point
#' 
#' @return \code{cobra}, an object of class COBRA, with potentially modified cobra$xbest, 
#'          cobra$fbest, cobra$ibest, cobra$noProgressCount
#' 
#' @keywords internal   
#' 
updateCobraEqu<-function(cobra,xNew){

  currentEps<-cobra$currentEps[(length(cobra$currentEps))]
  equMargin<-currentEps
  xNewIndex<-length(cobra$numViol)
  temp<-cobra$Gres
  #temp[,cobra$equIndex]<-temp[,cobra$equIndex]-currentEps

  #   # /WK/ this currentFeas-calculation has a bug for mixed equality-inequality constraints
  #   #      (the inequality constraints are compared to equMargin, but should be compared to 0)
  #   NEW_CURRENTFEAS=T
  #   if (!NEW_CURRENTFEAS) {
  #     temp<-cbind(temp,-temp[,cobra$equIndex])
  #     currentMaxViols<-sapply(1:nrow(temp) , function(i){
  #       y<-max(0,temp[i,])
  #       return(y)
  #     })    
  #     currentFeas<-which(currentMaxViols < equMargin )    
  #   } else {
  #
  # /WK/the correct new version: we check whether 
  #
  #          g_i(x) <= 0,  h_j(x) - equMargin <= 0,    -h_j(x) - equMargin <= 0
  #
  # for each row of cobra$Gres is valid and set only rows fulfilling this to currentFeas==TRUE
  temp<-cbind(temp,-temp[,cobra$equIndex])
  equ2Index <- c(cobra$equIndex,cobra$nConstraints+(1: length(cobra$equIndex)))
  temp[,equ2Index] <- temp[,equ2Index] - equMargin
  currentMaxViols<-sapply(1:nrow(temp) , function(i){
    y<-max(0,temp[i,])
    return(y)
  })  
  currentFeas<-which(currentMaxViols <= 0 )

  #   # /WK/2016/04/: The while-loop below was buggy for the following reasons: 
  #   #   - It could lead to a non-terminating loop for mixed equality-inequality constraints, 
  #   #     if one or several *inequalities* were the reason for infeasibility
  #   #   - It would not be clear, after increasing equMargin enough to have at least one equality-
  #   #     feasible solution, how to balance this with the infeasibility of the inequalities.
  #   while(length(currentFeas)==0){
  #     equMargin<-equMargin*1.1
  #     temp<-cobra$Gres
  #     temp<-cbind(temp,-temp[,cobra$equIndex])
  #     equ2Index <- c(cobra$equIndex,cobra$nConstraints+(1: length(cobra$equIndex)))
  #     temp[,equ2Index] <- temp[,equ2Index] - equMargin
  #     currentMaxViols<-sapply(1:nrow(temp) , function(i){
  #       y<-max(0,temp[i,])
  #       return(y)
  #     })  
  #     currentFeas<-which(currentMaxViols <= 0 )
  #   }
  #
  #   #   **Solution**: instead of the while-loop we follow now a simpler solution:
  #   #   If length(currentFeas)==0, we do the same thing as in updateSaveCobra:
  #   #   If the new point has a smaller maxViol then we take it as xbest, otherwise
  #   #   we leave the triple (xbest,fbest,ibest) as before (e.g. as set by cobraInit.R, 
  #   #   line 400: From all points with minimum number of violated constraints, take 
  #   #   the one with smallest Fres.)
  if (length(currentFeas)==0) {
    
    cobra$ibest<-which(currentMaxViols==min(currentMaxViols))[1]
    cobra$fbest<-cobra$Fres[cobra$ibest]
    cobra$xbest<-cobra$A[cobra$ibest,]
    return(cobra)
  }
  
  # Otherwise, if length(currentFeas)>0, we take among the feasible solutions the one with 
  # minimum objective Fres:
  fminInd<-which(cobra$Fres[currentFeas]==min(cobra$Fres[currentFeas]))
    # The new cobra$ibest might refer to a previous solution which differs from cobra$ibest 
    # so far and is NOT the ibest of the new point!
    # Why? - The so-far ibest might be a solution for an older equMargin band which is no longer 
    # valid in the current iteration. Then fminInd searches among the now valid solutions
    # the new best Fres. An older solution might come into play.
    # This is the reason why there can be in cobra$df a line where df$Best changes to a new value, 
    # but this value is NOT the df$y of the current iteration (as it used to be the case 
    # for inequality constraint handling).
  cobra$ibest<-currentFeas[fminInd[1]]
  cobra$xbest<-cobra$A[cobra$ibest,]
  cobra$fbest<-cobra$Fres[cobra$ibest]
  cobra$currentFeas<-currentFeas
  
  return(cobra)
}




#modifyMu
#'
#' Modify equality margin \eqn{\mu}
#' 
#' @param   Cfeas       counter feasible iterates
#' @param   Cinfeas     counter infeasible iterates
#' @param   Tfeas       threshold counts
#' @param   currentEps  current value for \eqn{\mu}
#' @param   cobra       list of class COBRA
#' 
#' @return  \code{currentEps}, the modified value for \eqn{\mu}
#' 
#' @keywords internal   
#' 
modifyMu<-function(Cfeas,Cinfeas,Tfeas,currentEps,cobra){
  alg<-cobra$equHandle$epsType
  switch(alg,
         expFunc={ # old funcCon function which is refered in equHandle-test file
           currentEps<-max(currentEps[length(currentEps)]/(cobra$equHandle$dec), cobra$equHandle$equEpsFinal)
         },
        SAexpFunc={ # Self-Adjsuting expFunc
          currentEps<-max(mean(c(currentEps[length(currentEps)]/(cobra$equHandle$dec),cobra$trueMaxViol[cobra$ibest]*cobra$finMarginCoef)), cobra$equHandle$equEpsFinal)
        },
         funcDim={
           currentEps<-(cobra$currentEps[1]*(1/cobra$equHandle$dec)^((nrow(cobra$A)-3*ncol(cobra$A))/((Tfeas^2)/2-1)))+cobra$equHandle$equEpsFinal
         },
         funcSDim={ 
           currentEps<-(cobra$currentEps[1]*(1/cobra$equHandle$dec)^((nrow(cobra$A)-3*ncol(cobra$A))/Tfeas))+cobra$equHandle$equEpsFinal
         },
         Zhang={
           currentEps<-max(currentEps*(1-(length(cobra$currentFeas)/nrow(cobra$A))),cobra$equHandle$equEpsFinal)
         },
         CONS={
           currentEps<-cobra$equHandle$equEpsFinal
         },
         
         "invalid expression"
  )
  
  return(currentEps)
}