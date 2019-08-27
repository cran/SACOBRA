#
#
#Samineh Bagheri, Patrick Koch, Wolfgang Konen
#Cologne University of Applied Science
#
#27.April.2014
#cobraPhaseI.R
#phase I
#Find a feasible point
#
#
#

#' Find a feasible solution.
#'
#' Find a feasible solution using the COBRA optimizer phase I by searching new infill points.
#' Please note that this phase can be skipped by setting the \code{cobra$skipPhaseI} parameter to \code{TRUE} in the initialization phase \code{\link{cobraInit}()}
#'
#' @param cobra an object of class COBRA, this is a (long) list containing all settings
#'      from \code{\link{cobraInit}}
#'
#' @return \code{cobra}, an object of class COBRA
#'  
#' @seealso   \code{\link{cobraPhaseII}}, \code{\link{cobraInit}}
#' @author Wolfgang Konen, Samineh Bagheri, Patrick Koch, Cologne University of Applied Sciences
#' @export
#' 
cobraPhaseI <- function(cobra){
  
  # -----------------------------------------------------------------------------------------------
  # ---------------- helper functions cobraPhaseI -------------------------------------------------
  # -----------------------------------------------------------------------------------------------
  
  ### these are now in innerFuncs.R
  
  # -----------------------------------------------------------------------------------------------
  # ---------------- end helper functions cobraPhaseI ---------------------------------------------
  # -----------------------------------------------------------------------------------------------
  
  
  print("no Feasible point available")
  print("PHASE I Started")
  phase<-"PHASE I"
  

  ##################################Initializing cobra parameters####################################
  #xbestIndex<-which(cobra$Fres==cobra$fbest)    # /WK/ we have cobra$ibest
  fn=cobra$fn
  dimension=cobra$dimension
  nConstraints <- cobra$nConstraints
  cobra$EPS <- rep(cobra$epsilonInit,cobra$nConstraints) # Initializing Margin for all constraints
  n<-nrow(cobra$A)
  #iteration<-cobra$iteration     # /WK/04/2016: obsolete
  feasibleSolutionExists<-(0 %in% cobra$numViol)
  #predY initialization should move to init code for phase I and phase II
  predY = rep(NA,cobra$initDesPoints) # structure to store surrogate optimization results
  cobra$predC = matrix(nrow=cobra$initDesPoints,ncol=cobra$nConstraints) # matrix to store predicted constraint values
  constraintPrediction = NULL # actual constraint value prediction
  optimizerConvergence = rep(1,cobra$initDesPoints) # vector to store optimizer convergence
  cobra$optimizationTime <- rep(0,cobra$initDesPoints)
  feas <- sapply(1:nrow(cobra$Gres), FUN = function(i) !any(cobra$Gres[i,]>0)) # feasibility of initial design
  #feval initialization should move to init code for phase I and phase II
  feval <- rep(NA,cobra$initDesPoints) # structure to store function evaluations on surrogate
  fbestArray <- rep(cobra$fbest,cobra$initDesPoints)
  penaF <- cobra$penaF    # /WK/
  sigmaD <- cobra$sigmaD;  # /WK/
  cobra$Ffeas <- cobra$fbest      # /WK/ needed?
  cobra$Fall <- min(cobra$Fres)   # /WK/ needed?
  cobra$nCobyla <- 0;               # only for ISRESCOBY
  cobra$nIsres <- 0;                # only for ISRESCOBY
  
  ########################################################################################################
  # STEP4:                                                                                               #
  #  While a feasible Point is not found do the loop                                                     #
  ########################################################################################################
  
  while(!feasibleSolutionExists && n<cobra$feval){   
    gc()
    #iteration<-iteration+1       # /WK/04/2016: obsolete, we have 'n'
    
    ################################################
    # STEP4.1: Update RBF for m constraints       #
    ################################################
    #constraintSurrogates <- list()  
    cobra$A <- as.matrix(cobra$A)
    cobra$Gres <- as.matrix(cobra$Gres)
    cobra$constraintSurrogates <- trainCubicRBF(cobra$A,cobra$Gres,ptail=cobra$ptail,squares=cobra$squares)
    cobra$fitnessSurrogate <- trainCubicRBF(cobra$A,cobra$Fres,ptail=cobra$ptail,squares=cobra$squares)
    

    ################################################
    # STEP4.2:Determine distance requirement(XI) #
    ################################################
    
    gama<-cobra$XI[(nrow(cobra$A) %% length(cobra$XI))+1]     
    cobra$ro<-gama*cobra$l
    
    ################################################
    # STEP4.3: select next point                   #
    ################################################
    #optimization sub problem
    
    ptm <- proc.time()
    subMin <- list()
    cat(cobra$seqOptimizer, " optimization on surrogate...\n")
    cobra <- checkIfCobraOptimizable(cobra);
    
    sw=switch(cobra$seqOptimizer,
           COBYLA={ subMin<- nloptr::cobyla(cobra$xbest,fn=subProb2PhaseI,lower=cobra$lower,upper=cobra$upper,hin=gCOBRAPhaseI,control=list(maxeval=cobra$seqFeval),cobra=cobra); subMin$feval=subMin$iter },
           ISRES ={ subMin<- isres2(cobra$xbest,fn=subProb2PhaseI,lower=cobra$lower,upper=cobra$upper,hin=gCOBRAPhaseI, maxeval=cobra$seqFeval, cobra=cobra); subMin$feval=subMin$iter},
           HJKB = { subMin<- dfoptim::hjkb(cobra$xbest,fn=subProbPhaseI,lower=cobra$lower,upper=cobra$upper,control=list(maxfeval=cobra$seqFeval),cobra=cobra) },
           NMKB = { subMin<- nmkb2(cobra$xbest,fn=subProbPhaseI,lower=cobra$lower,upper=cobra$upper, control=list(maxfeval=cobra$seqFeval,tol=cobra$seqTol),cobra=cobra) },
           ISRESCOBY = { subMin<-isresCobyla(cobra$xbest,fn=subProb2PhaseI,hin=gCOBRAPhaseI, cobra=cobra); subMin$feval=subMin$iter; cobra$nCobyla=subMin$nCobyla; cobra$nIsres=subMin$nIsres; "ok" },
           #ACTIVECMA  = { subMin <- ActiveOnePlusOneCMAES(subProbPhaseI, length(cobra$xbest), cobra$xbest, opts=list(esname="ActiveCMAES", lb=cobra$lower, ub=cobra$upper, maxFunEvals=cobra$seqFeval)); 
           #               subMin$convergence <- 1;
           #             }
           "InvalidOptimizer"
    )
    testit::assert(sprintf("cobraPhaseI: Wrong value %s for seqOptimizer",cobra$seqOptimizer),
                   sw!="InvalidOptimizer")
    cobra$optimizationTime <- c(cobra$optimizationTime, (proc.time()-ptm)[3])
    cat("Optimization time for", subMin$feval, " iterations:", cobra$optimizationTime[length(cobra$optimizationTime)], "seconds\n")  
    cat("Predicted infill value:",subMin$value,"\n")
    if (penaF[1]*penaF[2] < penaF[3])
      penaF[1] <- penaF[1]*penaF[2]

    CHECKDIST=T               
    if (CHECKDIST) {
      # check whether the solution returned fulfills the distance requirement
      # 
      ed = cobra$ro - distLine(subMin$par,cobra$constraintSurrogates$xp)
      violatedDist = which(ed>0)
      if (length(violatedDist)>0) {
        #
        # If distance requirement is not fulfilled, increase sigmaD[1] (the empirical penalty factor 
        # in fitFuncPenalRBF or subProbPhaseI). This influences currently only NMKB (not COBYLA),  
        # because sigmaD is only used in fitFuncPenalRBF or subProbPhaseI.
        #
        if(sigmaD[1]*sigmaD[2]<sigmaD[3]) {
          sigmaD[1] <- sigmaD[1]*sigmaD[2]
          cat("***   Increasing sigmaD to: ",sigmaD[1],"at iteration",n,"  ***\n")          
        }
      }
    }

    
    ################################################
    # STEP4.4: Evaluate real functions             #
    ################################################
    cobra$fe<-cobra$fe+1
    predY <- c(predY,subMin$value)
    feval <- c(feval,subMin$feval)
    optimizerConvergence <- c(optimizerConvergence,subMin$convergence)
    cobra$predC <- rbind(cobra$predC,constraintPrediction)
    xNew<-subMin$par
    #cat("X:",xNew,"\n")
    xNew <- unlist(sapply(1:cobra$dimension, FUN=function(i)max(cobra$lower[i],xNew[i])))
    xNew <- unlist(sapply(1:cobra$dimension, FUN=function(i)min(cobra$upper[i],xNew[i])))
    
    ### -- /WK/2016/04: rescaling is implicit now, I think
    ##rescaling
    #if(cobra$rescale){      # /WK/ bug fix: added missing cobra$rescale clause
    #  xNewTemp<-sapply(1:dimension , function(i){  
    #    y<-scales::rescale(xNew[i],to=c(cobra$originalLower[i],cobra$originalUpper[i]),from=c(cobra$lower[i],cobra$upper[i]))
    #    return(y)
    #  })
    #} else {
    #  xNewTemp <- xNew
    #}
    #xNewEval<-fn(xNewTemp)
    
    xNewEval<-fn(xNew)
    #cat("F/Max(C):",xNewEval[1],"/",max(xNewEval[2:length(xNewEval)]),"\n")
    
    ## /WK/ OLD version, did not reflect the case of mixed equality-inequality-constraints
    ##      and did not reflect cobra$conTol
    #newNumViol<-length(which((unlist(xNewEval)[-1])>0))  # number of constraint violations for new point
    #newMaxViol<-max(0,max((unlist(xNewEval)[-1])) )      # maximum violation
    
    ## /WK/ NEW version for calculation of newMaxViol and newNumViol
    ##
    if (cobra$equHandle$active){
      currentEps<-cobra$currentEps[1]
      temp<-xNewEval[-1]
      # /WK/the new version: we check whether 
      #
      #          g_i(x) <= 0,  h_j(x) - currentEps <= 0,    -h_j(x) - currentEps <= 0
      #
      # for the approximation newPredC with cobra$constraintSurrogates and set ev1$newNumViol to the
      # number of violated constraints.
      # NOTE that temp is also used for ev1$newMaxViol below.
      temp<-c(temp,-temp[cobra$equIndex])
      equ2Index <- c(cobra$equIndex,cobra$nConstraints+(1: length(cobra$equIndex)))
      temp[equ2Index] <- temp[equ2Index] - currentEps
      newNumViol<-length(which(temp > cobra$conTol)) # Calculating number of constraint Violations for new point #0 change to conTol

      M <-  max(0,max(temp))  # Calculating maximum violation
      if(M <= cobra$conTol) M=0 
      newMaxViol <- M  
    }
    else # i.e. if (!cobra$equHandle$active)
    {  
      newNumViol<-length(which(xNewEval[-1] > cobra$conTol)) # number of constraint violations for new point #0 change to conTol
      if((max(0,max((xNewEval[-1])) )) > cobra$conTol){      # maximum violation
        newMaxViol<-max(0,max((xNewEval[-1])) )  
      }else{
        newMaxViol<-0
      }
    } # (cobra$equHandle$active)
    feas = c(feas, newNumViol < 1 )
    
    
    ################################################
    # STEP4.5: Update Information                  #
    ################################################
    cobra$A<-rbind(cobra$A,xNew)
    cobra$Fres <- c(cobra$Fres, xNewEval[1])
    newConstraintLine <- xNewEval[2:(cobra$nConstraints+1)]#sapply(2:(m+1), FUN=function(i)(c(Gres[i],xNewEval[i])))
    cobra$Gres = rbind(cobra$Gres,newConstraintLine)
    cobra$numViol<-c(cobra$numViol,newNumViol)
    cobra$maxViol<-c(cobra$maxViol,newMaxViol)
    cobra$phase<-c(cobra$phase,phase)
    xNewIndex<-length(cobra$numViol)
    if(nrow(cobra$A) %% cobra$verboseIter == 0) {
      s<-sprintf("%s.[%d]: %f %f | %f | %f" , 
                 phase ,nrow(cobra$A), cobra$A[xNewIndex,1] ,cobra$A[xNewIndex,2] , xNewEval[1] , newMaxViol)
      print(s)
    }
    
    
    ################################################
    # STEP4.6: Update best Point                   #
    ################################################
    if( (cobra$numViol[xNewIndex] < cobra$numViol[cobra$ibest]) || 
          ((cobra$numViol[xNewIndex]==cobra$numViol[cobra$ibest]) 
           && (cobra$maxViol[xNewIndex]< cobra$maxViol[cobra$ibest])) ){
      
      cobra$xbest<-xNew
      cobra$fbest<-cobra$Fres[xNewIndex]
      cobra$ibest<-xNewIndex
    }
    
    
    
    cobra$fbestArray<-c(cobra$fbestArray,cobra$fbest)
    cobra$xbestArray<-rbind(cobra$xbestArray,cobra$xbest)
    feasibleIndices <- which(sapply(1:nrow(cobra$Gres),FUN=function(i)(all(cobra$Gres[i,]<0))))
    # xbestIndex<-which.min(cobra$Fres[feasibleIndices])                      # finding index of the best point so far
    
    n<-nrow(cobra$A)
    
    
    # result data frame
    df <- data.frame(cobra$Fres, 
                     predY, 
                     feas, 
                     feasPred=rep(NA,length(feas)),   # /WK/ bug fix
                     cobra$numViol,
                     cobra$maxViol,
                     feval, 
                     cobra$fbestArray,
                     #cobra$A,  
                     #cobra$Gres, 
                     #predC, 
                     cobra$seqOptimizer,
                     cobra$optimizationTime,
                     optimizerConvergence, 
                     row.names=NULL)
    colnames(df) = c("y", 
                     "predY", 
                     "feasible",
                     "feasPred",                # /WK/ bug fix
                     "nViolations",
                     "maxViolation",
                     "FEval",
                     "Best",
                     #sprintf("x%i",1:cobra$dimension),
                     #sprintf("c%i",1:cobra$nConstraints), 
                     #sprintf("predC%i",1:cobra$nConstraints),
                     "optimizer",
                     "optimizationTime",
                     "conv")
    df <- cbind(iter=1:nrow(df),df)
    df <- cbind(df,seed=cobra$cobraSeed)
    cobra$df <- df
    cobra$phase1DesignPoints <- nrow(df)
    
    # reset some 'inner' variables of cobra so that the list returned is simpler:
    if (!cobra$saveSurrogates) {
      cobra$constraintSurrogates <- NULL;
      cobra$fitnessSurrogate <- NULL
    }
    
    if (cobra$saveIntermediate) {
      # save intermediate results
      # cobraResult = list(cobra=cobra, df=df, constraintSurrogates=cobra$constraintSurrogates, fn=fn) 
      cobraResult = cobra
      if (is.na(file.info("results")$isdir)) dir.create("results")    # if directory "results" does not exist, create it
      save(cobraResult, file=sprintf("results/cobra-%s-%s-%i.RData",cobra$fName,cobra$seqOptimizer,cobra$cobraSeed))
    }
    

    
    feasibleSolutionExists<-(0 %in% cobra$numViol)
  }
  
  print("END of PHASE I")
  
  return(cobra)
  
  
}