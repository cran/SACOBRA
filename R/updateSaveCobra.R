# updateSaveCobra.R
#
#' Update and save cobra
#'
#' Helper for \code{\link{cobraPhaseII}}: make some assertion tests, 
#' update elements in object \code{cobra}, including data frames \code{df} and \code{df2}, 
#' and - if \code{cobra$saveIntermediate==TRUE} - save cobra in subdir \code{results/}.
#' Most importantly \code{cobra$xbest, cobra$fbest, cobra$ibest} are updated. They characterize the 
#' best feasible point (least violating point if no feasible point was found so far) and
#' influence the next starting point.
#'
#' Note: the elements \code{A, Fres, Gres} of \code{cobra} are set in \code{updateInfoAndCounters},
#' an internal function of \code{\link{cobraPhaseII}}.
#'
#' @param cobra an object of class COBRA, this is a (long) list containing all settings
#'        from \code{\link{cobraPhaseII}}
#' @param ev1   a list filled by calls to \code{\link{evalReal}}. We need here the elements xNew,
#'        feas, feasPred, feval, optimizerConvergence, optimizationTime, predY, predVal  
#' @param subMin   see \code{\link{cobraPhaseII}} 
#' @param sigmaD   see \code{\link{cobraPhaseII}} 
#' @param penaF   see \code{\link{cobraPhaseII}}
#' @param gama   see \code{\link{cobraPhaseII}} 
#' @param EPS   see \code{\link{cobraPhaseII}}
#' @param fitFuncPenalRBF   helper function from \code{\link{cobraPhaseII}}
#' @param distRequirement   helper function from \code{\link{cobraPhaseII}}
#' @param fitnessSurrogate  the model used for \code{predSoluFunc}
#'
#' @return \code{cobra}, an object of class COBRA, 
#'    enhanced here by the following elements (among others):
#'      \item{\code{df}}{  data frame with summary of the optimization run (see \code{\link{cobraPhaseII}})}
#'      \item{\code{df2}}{  data frame with additional summary information (see \code{\link{cobraPhaseII}})}
#'      \item{\code{dftr}}{  data frame with additional summary information for TR(see \code{\link{cobraPhaseII}})}
#'      \item{\code{fbest}}{ the best feasible objective value found }
#'      \item{\code{xbest}}{ the point in input space yielding the best feasible objective value }
#'      \item{\code{ibest}}{ the corresponding iteration number (row of \code{cobra$df}, of \code{cobra$A}) }
#'      \item{\code{fbestArray}}{ vector of all \code{fbest} }
#'      \item{\code{xbestArray}}{ vector of all \code{xbest} }
#'      
#' @keywords internal   
#' @seealso   \code{\link{cobraPhaseII}}, \code{\link{cobraInit}}
#'      
updateSaveCobra <- function(cobra,ev1,subMin,sigmaD,penaF,gama,EPS,
                            fitFuncPenalRBF,distRequirement,
                            fitnessSurrogate=cobra$fitnessSurrogate)
{
  xNew=ev1$xNew;
  feas=ev1$feas;
  feasPred=ev1$feasPred;
  feval=ev1$feval;
  optimizerConvergence=ev1$optimizerConvergence;
  predY=ev1$predY;
  predVal=ev1$predVal;
  df_RS  <-  cobra$df$RS
  #browser()
  diff<-nrow(cobra$A)-length(predY)
  #deprecated
  # if(cobra$LocalExpensiveSearch && nrow(cobra$A)>length(predY)){
  # 
  #   ev1$feas<-feas<-c(feas,rep(NA,diff))
  #   ev1$feasPred<-feasPred<-c(feasPred,rep(NA,diff))
  #   ev1$feval<-feval<-c(feval,rep(NA,diff))
  #   ev1$optimizerConvergence<-optimizerConvergence<-c(optimizerConvergence,rep(NA,diff))
  #   ev1$predY<-predY<-c(predY,rep(NA,diff))
  #   ev1$predVal<-predVal<-c(predVal,rep(NA,diff)) 
  #   ev1$optimizationTime<-c(ev1$optimizationTime,rep(NA,diff))
  #   df_RS  <-  c(cobra$df$RS,rep(NA,diff))
  # }
  
  ro<-gama*cobra$l
  df_RS  <-  c(df_RS,!all(cobra$xbest==cobra$xStart))
  if (cobra$DEBUG_XI) {
    df_fxStart <- c(cobra$df$fxStart,cobra$fn(cobra$xStart)[1])
    df_fxbest <-  c(cobra$df$fxbest,cobra$fn(cobra$xbest)[1])
    df_RS2 <-  c(cobra$df$RS2,cobra$DEBUG_RS)    
  }
  
  if (cobra$WRITE_XI) {
    sumViol <- distRequirement(xNew,cobra$fitnessSurrogate,ro)$sumViol
    if (is.null(cobra$df)) {
      df_XI <- c(rep(NA,cobra$initDesPoints),gama)
      df_XIsumViol <- c(rep(NA,cobra$initDesPoints),sumViol)
    } else {
      df_XI <- c(cobra$df$XI,gama)
      df_XIsumViol <- c(cobra$df$XIsumViol,sumViol)
    }
    # cat("gama:",gama,"\n")
    
  }
  xNewIndex<-length(cobra$numViol)
  
  #testit::assert(cobra$fbest<=min(cobra$fbestArray[which(cobra$numViol==0)]))
  
  ### /WK/ Bug fix: the 3 commented lines below were mind-buggingly complex and not 
  ### correct. In some cases, cobra$fbestArray (a.k.a df$Best) would *increase* in value!!
  ### Changed it to a simpler logic with the help of cobra$ibest, which is - if set - 
  ### the index to the so-far-best feasible solution. If there is no feasible solution yet 
  ### it is NULL.
  ###
  #    index<-which(cobra$Fres==cobra$fbest)
  #    if(cobra$numViol[tail(index,1)]==0){       #If the so-far-best is a feasible point
  #      if((cobra$numViol[xNewIndex]==0 )&&(cobra$Fres[xNewIndex] < cobra$fbest )){  #If xNew is feasible and even better
  #SB: this condition is always true therefore replaced with Samineh's suggestion   
  #if (!is.null(cobra$ibest)){   # if we have an everBestFeasible at index ibest and fitness value fbest
  #testit::assert(cobra$fn(cobra$A[cobra$ibest,])[1]==cobra$fbest)
  #Samineh's suggestion:
  
  
  if(cobra$equHandle$active && (length(cobra$equIndex)!=0)){
    # ... if we handle equality constraints by 'equMargin & two inequality constraints' 
    # and have equality constraints in the actual problem: 
    # calculate cobra$xbest,fbest,ibest via updateCobraEqu:
    cobra<-updateCobraEqu(cobra,xNew)      # in modifyEquCons.R
  }else if(cobra$numViol[cobra$ibest]==0){    # if the so-far best is feasible...
    #testit::assert(cobra$fn(cobra$A[cobra$ibest,])[1]==cobra$fbest)
    testit::assert(cobra$Fres[cobra$ibest]==cobra$fbest)
    if((cobra$numViol[xNewIndex]==0 )&&(cobra$Fres[xNewIndex] < cobra$fbest )){  #If xNew is feasible and even better
      cobra$xbest<-xNew
      cobra$fbest<-cobra$Fres[xNewIndex]
      cobra$ibest<-xNewIndex    
      #cobra$noProgressCount<-0
    }
  }else{                                # if we have not yet an everBestFeasible ...
    if(cobra$numViol[xNewIndex]==0){# the new point is feasible then we select it
      cobra$xbest<-xNew                 # ... take xNew, if it is feasible
      cobra$fbest<-cobra$Fres[xNewIndex]
      cobra$ibest<-xNewIndex
     # cobra$noProgressCount<-0
    }else{        #SB12.10.2015: the new solution is infeasible: look for the best infeasible solu
      if(cobra$maxViol[xNewIndex] < cobra$maxViol[cobra$ibest]){
        cobra$xbest<-xNew                 # ... take xNew, if it has smaller maxViol
        cobra$fbest<-cobra$Fres[xNewIndex]
        cobra$ibest<-xNewIndex
      #  cobra$noProgressCount<-0
      }
    }
    
  }
  # If we do not have a feasible point AND xNew has a larger maxMivol than ibest, then leave the 
  # triple (xbest,fbest,ibest) at the setting of cobraInit.R, line 400: From all points 
  # with minimum number of violated constraints, take the one with smallest Fres.
  
  #--- only debug ---
  #tn = fn(cobra$xbest)[-1]
  #cat(sprintf("max.viol(xbest)=%9.4g\n",max(0,max(tn))))
 
  cobra$fbestArray<-c(cobra$fbestArray,cobra$fbest)
  cobra$xbestArray<-rbind(cobra$xbestArray,cobra$xbest)
  feasibleIndices <- which(sapply(1:nrow(cobra$Gres),FUN=function(i)(all(cobra$Gres[i,]<0))))
  xbestIndex<-which.min(cobra$Fres[feasibleIndices])                      # finding index of the best point so far
  
  # only diagnostics, needed for cobra$df & cobra$df2 /WK/
  solu <- cobra$solu; 
  if (is.null(solu)) {
    solu=subMin$par;
    soluOrig <- inverseRescale(subMin$par,cobra);
  } else {
    if (cobra$rescale) 
      if (is.matrix(solu)) {
        solu <- t(sapply(1:nrow(solu),function(i){ forwardRescale(solu[i,],cobra)}))
      } else {
        solu <- forwardRescale(solu,cobra);
      }
    soluOrig <- cobra$solu;
  }
  # now solu is always in *rescaled* input space
  
  predSoluFunc <- function(x)getPredY0(x,fitnessSurrogate,cobra);
  if (is.matrix(solu)) {      # in case of multiple global optima in solu:
    predSolu <- sapply(1:nrow(solu),function(i){ predSoluFunc(solu[i,])}) ;
    predSoluPenal <- sapply(1:nrow(solu),function(i){ fitFuncPenalRBF(solu[i,])}) ;
  } else {
    predSolu <- predSoluFunc(solu);      
    predSoluPenal <- fitFuncPenalRBF(solu);      
  }
  predSolu <- min(predSolu)   # Why min? - In case of multiple global optima: predSolu is the 
  # value of fitFuncPenalRBF at the best solution solu
  predSoluPenal <- min(predSoluPenal) 
  if (is.null(cobra$df)) {
    df_predSolu <- c(rep(NA,cobra$initDesPoints),predSolu)
  } else {
    df_predSolu <- c(cobra$df$predSolu,predSolu)
    
  }
  #deprecated
  # if(cobra$LocalExpensiveSearch && nrow(cobra$A)>length(df_predSolu)){
  #   df_predSolu<-c(df_predSolu,rep(NA,diff))
  #   df_XI <- c(df_XI,rep(NA,diff))
  #   df_XIsumViol <- c(df_XIsumViol,rep(NA,diff))
  # }
  
  origA = t(sapply(1:nrow(cobra$A),function(i){ inverseRescale(cobra$A[i,],cobra) }))
  if (is.matrix(solu)) {          # this is for the case with multple solutions (like in G11)
    da=sapply(1:nrow(solu),function(i){ distLine(solu[i,],cobra$A) })
    # da has nrow(solu) columns, each column has the distance of the ith solu to all points cobra$A.
    # Select this column which has the element with minimum distance.
    ind=which.min(apply(da,2,min))
    distA = da[,ind]
    do=sapply(1:nrow(soluOrig),function(i){ distLine(soluOrig[i,],origA) })
    distOrig = do[,ind]
  } else {
    distA = distLine(solu,cobra$A)       # distance in rescaled space, distLine: see RbfInter.R
    distOrig = distLine(soluOrig,origA)  # distance in original space
  }
  
  # result data frame
  testit::assert("predY",length(cobra$Fres)==length(predY))
  if(cobra$CONSTRAINED){
  testit::assert("feas",length(cobra$Fres)==length(feas))
  testit::assert("feasPred",length(cobra$Fres)==length(feasPred))
  testit::assert("cobra$numViol",length(cobra$Fres)==length(cobra$numViol))
  testit::assert("cobra$maxViol",length(cobra$Fres)==length(cobra$maxViol))
  }
  testit::assert("cobra$fbestArray",length(cobra$Fres)==length(cobra$fbestArray))
  testit::assert("ev1$optimizationTime",length(cobra$Fres)==length(ev1$optimizationTime))
  testit::assert("optimizerConvergence",length(cobra$Fres)==length(optimizerConvergence))
  testit::assert("optimizerConvergence",length(df_predSolu)==length(optimizerConvergence))
  if(cobra$CONSTRAINED){
     df <- data.frame(y=cobra$Fres, 
                   predY=predY,           # surrogate fitness
                   #predSolu=interpRBF(solu,cobra$fitnessSurrogate),  # OLD and WRONG
                   predSolu=df_predSolu,
                   feasible=feas, 
                   feasPred=feasPred,
                   nViolations=cobra$numViol,
                   maxViolation=cobra$maxViol,
                   FEval=feval, 
                   Best=cobra$fbestArray,
                   optimizer=rep(cobra$seqOptimizer,length(cobra$Fres)),
                   optimizationTime=ev1$optimizationTime,
                   conv=optimizerConvergence,
                   dist=distA,
                   distOrig=distOrig,
                   RS=df_RS,           # TRUE, if it is an iteration with random start point
                   row.names=NULL
      )
  }else{
    if(is.null(cobra$df$realfeval))
      realfeval<-c()
    else
      realfeval<-cobra$df$realfeval
    
    
    df <- data.frame(y=cobra$Fres, 
                     predY=predY,           # surrogate fitness
                     predSolu=df_predSolu,
                     feasible=T,
                     FEval=feval, 
                     realfeval=c(realfeval,nrow(get("ARCHIVE",envir=intern.archive.env))),
                     Best=cobra$fbestArray,
                     optimizer=rep(cobra$seqOptimizer,length(cobra$Fres)),
                     optimizationTime=ev1$optimizationTime,
                     conv=optimizerConvergence,
                     dist=distA,
                     distOrig=distOrig,
                     RS=df_RS,           # TRUE, if it is an iteration with random start point
                     row.names=NULL
    )
    
  }
  
  if (cobra$WRITE_XI) {
    df$XI=df_XI
    df$XIsumViol=df_XIsumViol
  }
  
  if (cobra$DEBUG_XI) {
    firstSolu <- solu;
    if (is.matrix(solu)) firstSolu <- solu[1,]; 
    optimum <- cobra$fn(firstSolu)[1];      
    
    df$fxStart=df_fxStart # objective function at xStart
    df$fxbest=df_fxbest   # objective function at xbest
    df$exbest=df_fxbest - optimum   # error (objective function - optimum) at xbest
    df$RS2=df_RS2         # the same
    df$iter2=1:nrow(df)
    df$errFy=df$y - optimum # the error of the optimizer result in every iteration
    #if(tail(df_RS,1)==TRUE) browser()
    #browser()
    testit::assert(df$RS==df_RS2)
    if (any(df$fxbest[!df$RS]!=df$fxStart[!df$RS])) {
      browser()
      df$fxbest[!df$RS]=df$fxStart[!df$RS]  # symptomatic fix for the next assert
    }
    testit::assert(df$fxbest[!df$RS]==df$fxStart[!df$RS])
  }
  df <- cbind(iter=1:nrow(df),df)
  df <- cbind(df,seed=cobra$cobraSeed)
  cobra$df <- df
  cobra$df2 <- rbind(cobra$df2,data.frame(
    iter=tail(df$iter,1),
    predY=tail(predY,1),           # surrogate fitness at current point xNew  
    predVal=tail(predVal,1),       # surrogate fitness + penalty at xNew
    predSolu=predSolu,             # surrogate fitness at solu (only diagnostics). 
    predSoluPenal=predSoluPenal,   # surrogate fitness + penalty at solu (only diagnostics). 
    sigmaD=sigmaD[1],
    penaF=penaF[1],
    XI=gama,
    fBest=tail(df$Best,1),
    EPS=EPS[1]
    #,fBest2=fitFuncPenalRBF(xNew)    # not the same as predVal, since penaF or sigmaD might have changed (!)
    #,fSolu=fitFuncPenalRBF(min(solu))# not the same as predSolu for the same reason  
    #,feas=feas,
    #,data.frame(xNew,row.names=NULL)
  ))
  
  cobra$dftr<-rbind(cobra$dftr,data.frame(
    TRiter=cobra$TRiter,
    TRapprox=cobra$TRapprox,
    TFR=cobra$TFR,
    TSR=cobra$TSR,
    Fratio=cobra$Fratio,
    TRpop=cobra$TRpop,
    TRdelta=cobra$TRdelta
  ))
  
  # consistency check for data frames df and df2:
  msg <- "updateSaveCobra: wrong nrow for df and df2";
  if (is.null(cobra$phase1DesignPoints)) {
 #   testit::assert(msg,nrow(cobra$df)==nrow(cobra$df2)+cobra$initDesPoints)
  } else {
  #  testit::assert(msg,nrow(cobra$df)==nrow(cobra$df2)+cobra$phase1DesignPoints)      
  }
  
  if (cobra$saveIntermediate) {
    # save intermediate results
    # cobraResult = list(cobra=cobra, df=df, constraintSurrogates=cobra$constraintSurrogates, fn=fn) 
    cobraResult = cobra
    if (is.na(file.info("results")$isdir)) dir.create("results")    # if directory "results" does not exist, create it
    save(cobraResult, file=sprintf("results/cobra-%s-%s-%i.RData",cobra$fName,cobra$seqOptimizer,cobra$cobraSeed))
  }
  cobra$ev1<-ev1
  return(cobra)
} # updateSaveCobra()
