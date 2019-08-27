#
#Samineh Bagheri, Patrick Koch, Wolfgang Konen
#Cologne University of Applied Sciences
#
#April,2014 - Nov,2015
#cobraPhaseII.R
#
#' Improve the feasible solution by searching new infill points
#'
#' Improve the feasible solution using the SACOBRA optimizer phase II
#' by searching new infill points with the help of RBF surrogate models. 
#' May be even called if no feasible solution is found yet, then phase II will try to find
#' feasible solutions. \cr
#' The problem to solve iteratively is: \cr
#' \deqn{ \mbox{Minimize}\quad  f(\vec{x}) , \vec{x} \in [\vec{a},\vec{b}] \subset \mathbf{R}^d }
#' \deqn{ \mbox{subject to}\quad g_i(\vec{x}) \le 0, i=1,\ldots,m    }
#' \deqn{ \mbox{~~~~~~~~~~}\quad\quad h_j(\vec{x}) = 0, j=1,\ldots,r.    } \cr
#' In this phase the main optimization steps are repeated in a loop as long as the budget is not exhausted.
#' In every iteration the surrogate models are updated and an optimization on the surrogates is done in order 
#' to find a better feasible solution.
#'
#' @param cobra an object of class COBRA, this is a (long) list containing all settings
#'        from \code{\link{cobraInit}}
#'
#' @return \code{cobra}, an object of class COBRA from \code{\link{cobraInit}}, 
#'    enhanced here by the following elements (among others):
#'      \item{\code{fn}}{ function accepting a \code{d}-dimensional vector \eqn{\vec{x}} and
#'            returning an \code{(1+m+r)}-vector \code{c(}\eqn{f,g_1,\ldots,g_m,h_1,\ldots,h_r}\code{)}. This
#'            function may be a rescaled and plog-transformed version of the original \code{fn} 
#'            passed into \code{\link{cobraInit}}. The original \code{fn} is stored in 
#'            \code{cobra$originalFn}. }
#'      \item{\code{df}}{  data frame with summary of the optimization run (see below)}
#'      \item{\code{df2}}{  data frame with additional summary information (see below)}
#'      \item{\code{dftr}}{  data frame with additional summary information for TR (see below)}
#'      \item{\code{A}}{ \code{(feval x d)}-matrix containing all evaluated points 
#'            in input space. If rescale==TRUE, all points are in \strong{rescaled} input space. }
#'      \item{\code{Fres}}{ a vector of the objective values of all evaluated points }
#'      \item{\code{Gres}}{ a \code{(feval x m)}-matrix of the constraint values of all evaluated points }
#'      \item{\code{predC}}{ a \code{(feval x m)}-matrix with the prediction of  
#'            \code{cobra$constraintSurrogates} at all evaluated points  }
#'      \item{\code{fbest}}{ the best feasible objective value found }
#'      \item{\code{xbest}}{ the point in input space yielding the best feasible objective value }
#'      \item{\code{ibest}}{ the corresponding iteration number (row of cobra$df, of cobra$A)}
#'      \item{\code{PLOG}}{ If TRUE, then the objective surrogate model is trained on the 
#'            \code{\link{plog}}-transformed objective function. }
#Note   
#'   Note that \code{cobra$Fres}, \code{cobra$fbest}, \code{cobra$fbestArray} and similar contain 
#'   always the objective values of the orignial function \code{cobra$fn[1]}. (The surrogate models 
#'   may be trained on a \code{\link{plog}}-transformed version of this function.)
#'   
#'   \code{feval} = \code{cobra$feval} is the maximum number of function evaluations.\cr 
#'   
#'   The data frame \code{cobra$df} contains one row per iteration with columns 
#'   \describe{
#'      \item{iter}{ iteration index }
#'      \item{y}{   true objective value \code{Fres} }
#'      \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predY} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{predSolu}{  surrogate objective value at best-known solution \code{cobra$solu}, if given. 
#'            If \code{cobra$solu} is NULL, take the current point instead. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predSolu} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{feasible}{ boolean indicating the feasibiltiy of infill point }
#'      \item{feasPred}{ boolean indicating if each infill point is feasible for \code{cobra$constraintSurrogates} }
#'      \item{nViolations}{ number of violated constraints }
#'      \item{maxViolation}{ maximum constraint violation. }
#'      \item{FEval}{  number of function evaluations in sequential optimizer. NA if it was a repair step }
#'      \item{Best}{  ever-best feasible objective value \code{fbest}. As long as there is 
#'            no feasible point, take among those with minimum number of violated constraints the
#'            one with minimum Fres. }
#'      \item{optimizer}{ e.g. "COBYLA"  }
#'      \item{optimizationTime}{  in sec}
#'      \item{conv}{ optimizer convergence code }
#'      \item{dist}{ distance of the current point (row of \code{cobra$A}) to the true solution 
#'            \code{cobra$solu} in rescaled space. If there is more than one solution, take the one
#'            which has the minimum distance element (since this is the solution to which the 
#'            current run converges). }
#'      \item{distOrig}{ same as \code{dist}, but in original space  }
#'      \item{XI}{  the DRC element used in the current iteration }
#'      \item{seed}{ the used seed in every run }
#'   }
#'
#'   The data frame \code{cobra$df2} contains one row per phase-II-iteration with columns 
#'   \describe{
#'      \item{iter}{ iteration index}
#'      \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predY} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{predVal}{   surrogate objective value + penalty }
#'      \item{predSolu}{   surrogate objective value at true solution (see \code{cobra$df$predSolu}) }
#'      \item{predSoluPenal}{   surrogate objective value + penalty at true solution (only diagnostics)}
#'      \item{sigmaD}{ the sigmaD element used in the current iteration (see \code{\link{cobraInit}})  }
#'      \item{penaF}{ penalty factor used in the current iteration (see \code{\link{cobraInit}}) }
#'      \item{XI}{  the DRC element used in the current iteration }
#'      \item{EPS}{ the current used margin for constraint function modeling (see \code{epsilonInit} in \code{\link{cobraInit}} ) }
#'   }
#'
#' @examples 
#' ## Initialize cobra. The problem to solve is the unconstrained sphere function sum(x^2).   
#'  
#' ## In version 1.1 and higher there is no need for defining a dummy 
#' ## constraint function for the unconstrained problems
#' d=2
#' fName="sphere"
#' cobra <- cobraInit(xStart=rep(5,d), fName=fName,
#'                    fn=function(x){c(obj=sum(x^2))},  
#'                    lower=rep(-10,d), upper=rep(10,d), feval=40)
#'                    
#' ## Run cobra optimizer
#' cobra <- cobraPhaseII(cobra)
#' 
#' ## The true solution is at solu = c(0,0)
#' ## where the true optimum is fn(solu)[1] = optim = 0
#' ## The solution found by SACOBRA:
#' print(getXbest(cobra))
#' print(getFbest(cobra))
#' 
#' ## Plot the resulting error (best-so-far feasible optimizer result - true optimum)
#' ## on a logarithmic scale:
#' optim = 0
#' plot(cobra$df$Best-optim,log="y",type="l",ylab="error",xlab="iteration",main=fName)
#' 
#' @seealso   \code{\link{cobraPhaseI}}, \code{\link{cobraInit}}
#' @author Wolfgang Konen, Samineh Bagheri, Patrick Koch, Cologne University of Applied Sciences
#' @export

##########################################################################################################
# Some rules about the COBRA-II-code:
#
# - cobra$df contains one row for each iteration, including initial points (see updateSaveCobra.R)
#   cobra$df2 contains one row for each phase-II iteration only (see updateSaveCobra.R)
# - cobra$PLOG is set by adFit (SACOBRA.R) and adFit is called at the start of each iteration
#   (see trainSurrogates.R)
# - cobra$solu is always in original input space. But when it is used (for diagnostics) in 
#   updateSaveCobra, then a local copy \code{solu} is made, and - if cobra$rescale==TRUE - 
#   \code{solu} is transformed to the rescaled space.
# - the rescaled space has the bounds [rep(cobra$newlower,d),rep(cobra$newupper,d)], (usually -1 and 1) 
# - cobra$xbest,cobra$fbest,cobra$ibest refer always to the same infill point (same iteration).
# - What is cobra$fbest before a feasible infill is found? - See \code{updateSaveCobra}:
#   If the new infill point is feasible, take its fn[1]-value (of course!). If it is not feasible,
#   leave it at the setting from cobraInit.R#400: from all points of the initial design
#   with minimum number of violated constraints, take the one with smallest Fres.
# 
##########################################################################################################

cobraPhaseII <- function(cobra){
  gc()
  verboseprint(cobra$verbose, important=FALSE,"There is at least one feasible point in the population Or PhaseI is skipped")
  verboseprint(2, important=TRUE,"PHASE II Started")
  phase<-"PHASE II"
  
  testit::assert("cobraPhaseII: cobra$fbest is NULL!",!is.null(cobra$fbest))
  testit::assert("cobraPhaseII: cobra$ibest is NULL!",!is.null(cobra$ibest))
  
  ########################################################################################################
  # STEP5: 
  # Initializing the parameters and
  #  Initialize the margin  and counters                                                                 #
  ########################################################################################################
  #cobra$RMSE<-c()
  #cobra$RMSEC<-NULL
  cobra$TM<-diag(cobra$dimension)
  newErr1<-0
  newErr2<-0
  err1<-c()
  err2<-c()
  ev1<- list()
  cobra$ApproxFrame<-list()
  cobra$selectedModel<-NULL
  cobra$MSXI<-c()
  cobra$WRITE_XI<-TRUE
  globalOptCounter<-1
  fn=cobra$fn
  dimension=cobra$dimension
  nConstraints <- cobra$nConstraints
  CHECKDIST=T       
  Tfeas <- cobra$Tfeas
  Tinfeas <- cobra$Tinfeas
  Cfeas<-0                # Starting Counters
  Cinfeas<-0
  EPS <- cobra$epsilonInit 
  num <- nrow(cobra$A)
  cobra$hesse<-diag(cobra$dimension)
  nRepair<- 0
  if(num==cobra$initDesPoints){
    ev1$predY = rep(NA,cobra$initDesPoints) # structure to store surrogate optimization results
    ev1$predVal = rep(NA,cobra$initDesPoints) 
    if(cobra$nConstraints!=0){
      ev1$predC = matrix(nrow=cobra$initDesPoints,ncol=cobra$nConstraints) # matrix to store predicted constraint values
      colnames(ev1$predC) <- paste("C",1:cobra$nConstraints,sep="") 
      ev1$feas <- sapply(1:nrow(cobra$Gres), FUN = function(i) !any(cobra$Gres[i,]>0)) # feasibility of initial design      
    }
    
    ev1$feasPred <- rep(NA,cobra$initDesPoints)
    ev1$optimizerConvergence = rep(1,cobra$initDesPoints) # vector to store optimizer convergence
    ev1$optimizationTime <- rep(0,cobra$initDesPoints)
    ev1$feval <- rep(NA,cobra$initDesPoints) # structure to store function evaluations on surrogate
    #cobra$fbestArray <- rep(cobra$fbest,cobra$initDesPoints) # /WK/ obsolete
  }else{
    ev1$predY<-cobra$df$predY
    ev1$predVal<-cobra$df$predVal
    #cobra$predC<-cobra$predC
    ev1$predC = cobra$predC
    ev1$feas<-cobra$df$feasible
    ev1$feasPred<-cobra$df$feasPred
    ev1$optimizerConvergence <- cobra$df$conv
    ev1$optimizationTime <- cobra$optimizationTime
    ev1$feval <- cobra$df$FEval
    #fbestArray <- cobra$fbestArray                     # /WK/ obsolete
  }
  attr(ev1,"state") <- "initialized"
  
  ## these variables are deprecated, we use now cobra$constraintSurrogates, 
  ## cobra$fitnessSurrogate* so that trainSurrogates can be in a separate 
  ## file trainSurrogates.R
  constraintSurrogates <- NULL # have them as variables on the global level of cobraPhaseII  
  fitnessSurrogate     <- NULL # such that inner trainSurrogates() can access them with "<<-"
  fitnessSurrogate1    <- NULL # built model according to the original value of the fitness values
  fitnessSurrogate2    <- NULL # built model according to the plog transformed of the fitness values
  
  constraintPrediction = NULL # actual constraint value prediction
  penaF <- cobra$penaF    
  sigmaD <- cobra$sigmaD;  
  cobra$important<-FALSE
  cobra$nCobyla <- 0;               # only for ISRESCOBY
  cobra$nIsres <- 0;                # only for ISRESCOBY
  #SB equality constraint handling 02.10.2015
  currentEps<-cobra$currentEps[1]
  
  ##------------------------------Trust Region variables-----------------------------|
  cobra$globalDevice<-NA
  cobra$localDevice<-NA
  cobra$TRcounter<-0
  cobra$TFR<-c() #min-max range of the fitness for the points in the trust region (gets updated in trustRegion)
  cobra$TSR<-c() #min-max range of the surrogate fitness for the points in the trust region (gets updated in trustRegion)
  cobra$Fratio<-c() #ratio of Frange nad Srange. SR/FR (when this ratio is verylarge it means model is oscilating)
  cobra$TRpop<-c()  #Size of the population in the trust region
  cobra$TRdelta<-c()
  cobra$TRapprox<-c()
  cobra$TRiter<-c()
  cobra$iA<-cobra$A
  cobra$tCenter<-cobra$xbest
  # -----------------------------------------------------------------------------------------------
  
  
  if (num >= cobra$feval) warning("num is after Phase I equal or larger than cobra$feval")
  
  # -----------------------------------------------------------------------------------------------
  # ---------------- helper functions cobraPhaseII ------------------------------------------------
  # -----------------------------------------------------------------------------------------------
  
  ### --- should later go into innerFuncs, but think about EPS and ro and fn
  
  fitFuncPenalRBF <- function(x) { 
    if(any(is.nan(x))){
      warning("fitFuncPenalRBF: x value is NaN, returning Inf")
      return(Inf)
    }
    
    y = interpRBF(x, cobra$fitnessSurrogate)
    if (cobra$trueFuncForSurrogates) y<-fn(x)[1]
    penalty<-0
    if(cobra$CONSTRAINED){
      constraintPrediction <-  interpRBF(x,cobra$constraintSurrogates) +EPS^2
      if (cobra$trueFuncForSurrogates) constraintPrediction <-  fn(x)[-1]+EPS^2
      violatedConstraints = which(constraintPrediction>0)
      penalty = sum(constraintPrediction[violatedConstraints]) 
    }
    
    penalty = penalty + distRequirement(x,cobra$fitnessSurrogate,cobra$ro)$sumViol*sigmaD[1]
    return(y + penalty*penaF[1])
  }
  
  
  
  distRequirement<- function(x,fitnessSurrogate,ro) {
    ed = ro - distLine(x,fitnessSurrogate$xp)
    violatedDist = which(ed>0)
    sumViol = sum(ed[violatedDist])
    return(list(ed=ed,   # vector of euclidean distances
                violatedDist=violatedDist,
                sumViol=sumViol))
  }

  # check whether the solution returned in subMin fulfills the distance requirement.
  # If not, increase sigmaD[1]. Return sigmaD.
  # 
  checkDistanceReq <- function(subMin,fitnessSurrogate,ro,sigmaD,CHECKDIST) {
    if (CHECKDIST) {
      sumViol = distRequirement(subMin$par,fitnessSurrogate,ro)$sumViol
      if (sumViol>0) {
        #
        # If distance requirement is not fulfilled, increase sigmaD[1] (the empirical penalty factor 
        # in fitFuncPenalRBF or subProb). This influences currently only NMKB and HJKB (not COBYLA),  
        # because sigmaD[1] is only used in fitFuncPenalRBF or subProb.
        #
        if(sigmaD[1]*sigmaD[2]<sigmaD[3]) {
          sigmaD[1] <- sigmaD[1]*sigmaD[2]
          verboseprint(cobra$verbose, important=FALSE,paste("***   Increasing sigmaD to: ",sigmaD[1],"at iteration",num,"  ***"))          
        }
      }
    }
    return(sigmaD)
  }
      
  # update cobra information (A, Fres, Gres)
  # update counters Cfeas, Cinfeas on the global level of cobraPhaseII
  #
  updateInfoAndCounters <- function(cobra,xNew,xNewEval,newNumViol,
                                    newMaxViol,trueMaxViol,phase,currentEps)
  {
    cobra$A<-rbind(cobra$A,xNew)
    cobra$TA<-rbind(cobra$TA,xNew)
    cobra$Fres <- c(cobra$Fres, xNewEval[1])
    #SB 02.10.2015 equality handling
    cobra$Gres = rbind(cobra$Gres,xNewEval[-1])
    cobra$currentEps<-c(cobra$currentEps,currentEps)
    cobra$numViol<-c(cobra$numViol,newNumViol)
    cobra$maxViol<-c(cobra$maxViol,newMaxViol)
    cobra$trueMaxViol<-c(cobra$trueMaxViol,trueMaxViol)

    
    cobra$phase<-c(cobra$phase,phase)
    if(nrow(cobra$A) %% cobra$verboseIter == 0){#important to print
      cobra$important=TRUE
    }else{
      cobra$important=FALSE
    }
    xNewIndex<-length(cobra$numViol)
    DEBUGequ<-FALSE
    if(cobra$equHandle$active && cobra$verbose==2)DEBUGequ=TRUE
    #if(newNumViol==0)browser(expr=(equMargin-newMaxViolOriginal <0))

    verboseprint(cobra$verbose, important=DEBUGequ,(sprintf("%s.[%d]: %f %f | %f | %f|[%f]" 
                  , phase ,nrow(cobra$A), cobra$A[xNewIndex,1] ,cobra$A[xNewIndex,2] , xNewEval[1] , newMaxViol,currentEps)))
    #SB: I added another print line which prints the best found solution after several iterations, if we do not have any interest for this the following lines can be commented
    realXbest<-inverseRescale(cobra$xbest,cobra)
    #realXbest<-sapply(1:length(cobra$xbest) , function(i){scales::rescale(cobra$xbest[i],from=c(cobra$newlower,cobra$newupper),to=c(cobra$lower[i],cobra$upper[i]))})
    
    if(cobra$equHandle$active){
     # browser()
      verboseprint(cobra$verbose, important=cobra$important,
                   sprintf("%s.[%d]: %f %f | %f | %f |[%f]", 
                           "Best Result" ,nrow(cobra$A), realXbest[1] ,realXbest[2] , cobra$fbest[1] , cobra$trueMaxViol [cobra$ibest],(currentEps)))
      
    }else{
      verboseprint(cobra$verbose, important=cobra$important,
                   sprintf("%s.[%d](%d): %f %f | %f | %f " ,
                           "Best Result" ,nrow(cobra$A),nrow(get("ARCHIVE",envir=intern.archive.env)), realXbest[1] ,realXbest[2] , cobra$fbest[1] , cobra$maxViol[cobra$ibest]))   
    }
    
    
    
    #browser()
    #check if the new point is feasible #Now the question is feasible on what???on shifted constraints
 
    #  browser(expr=cobra$numViol[xNewIndex]!=0)
    if(cobra$numViol[xNewIndex]==0){ 
      Cfeas <<- Cfeas+1
      Cinfeas <<- 0
    }else{
      Cinfeas <<- Cinfeas+1
      Cfeas <<- 0
    }
    return(cobra)
  }  
  
  
  # adjust margins (EPS)
  # may change EPS, currentEps, Cfeas, and Cinfeas on the global level of cobraPhaseII
  #
  adjustMargins <- function(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,epsMax,currentEps) {
    if(Cfeas >= Tfeas){
      EPS <<- EPS/2
      verboseprint(cobra$verbose, important=FALSE,sprintf("reducing epsilon to %f",EPS[1]))
      verboseprint(cobra$verbose, important=FALSE,sprintf("reducing equality margin to %f",currentEps))
      
      Cfeas <<- 0
    }
    
    if(Cinfeas>=Tinfeas){
      EPS <<- pmin(2*EPS,epsMax) 
      #currentEps <<- min(currentEps*(cobra$equHandle$inc),cobra$currentEps[1])
      verboseprint(cobra$verbose, important=FALSE,sprintf("increasing epsilon to %f",EPS[1]))
      verboseprint(cobra$verbose, important=FALSE,sprintf("increasing equality margin to %f",currentEps))
      
      Cinfeas <<- 0
    }
    
    if(cobra$equHandle$active)currentEps<<-modifyMu(Cfeas,Cinfeas,Tfeas,currentEps,cobra)
  }


  # -----------------------------------------------------------------------------------------------
  # ---------------- end helper functions cobraPhaseII --------------------------------------------
  # -----------------------------------------------------------------------------------------------
  
  ########################################################################################################
  # STEP6:                                                                                               #
  #  Improve the feasible point                                                                          #
  ########################################################################################################
  cobra$printP<-TRUE
  testit::assert("ev1 not initialized",attr(ev1,"state")=="initialized");
  # prior to entering while-loop, ev1 must have attribute state = "initialized"

  while(num < cobra$feval){


    ##########################################################
    # STEP6.1: UPDATE RBF MODEL for fitness and constraints  #
    ##########################################################
    gama<-cobra$XI[((globalOptCounter) %% length(cobra$XI))+1]#tempo

    if(cobra$MS$active && nrow(cobra$A) %% cobra$MS$freq==0){
      
      cobra<-selectModel(models=cobra$MS$models,
                        widths=cobra$MS$widths,cobra=cobra,freq=cobra$MS$freq,
                        slidingW=cobra$MS$slidingW,WinS=cobra$MS$WinS,
                        quant=cobra$MS$quant,gama)
      
    }
    
    #as long as no transformation of the fitness function is done, the models are simply generated as follows. 
     if(!cobra$TFlag) cobra <- trainSurrogates(cobra);  # side effect: cobra$constraintSurrogates, cobra$fitnessSurrogate
     
   # browser()
    #whitening transfomartion
    if(cobra$CA$active && any(cobra$CA$ITER<=nrow(cobra$A))){
      cobra$TFlag<-T  #transformation flag
      VALIDTRANSORFMATION<-T
      
      
      
      #surrogate of fitness function
      sfunc<-function(x){
        y<-predict.RBFinter(cobra$fitnessSurrogate,matrix(x,ncol=cobra$dimension))
        return(y)
      }
      
      #fitness function
      fitfunc<-function(x){
        y<-cobra$fn(x)
        return(y[1])
      }
      sortPop<-function(pop,Fres){
        temp<-cbind(pop,Fres)
        names(temp)<-c(sprintf(paste("x",c(1:ncol(pop)),sep="")),"Fres")
        temp<-temp[order(Fres),]
        return(list(pop=temp[,-ncol(temp)],Fres=temp[,ncol(temp)]))
      }
      
      #only every cobra$CA$ITER iterations the hessian matrix is updated
      if(any(cobra$CA$ITER==nrow(cobra$A))){

      switch(cobra$CA$HessianType,
             "real"={hesse<-numDeriv::hessian(fitfunc,x=cobra$xbest+runif(cobra$dimension,min=-0.001,max=0.001))}
            # "realp"={hesse<-pracma::hessian(f=fitfunc,x0=cobra$xbest+runif(cobra$dimension,min=-0.001,max=0.001))}, #deprecated because the results are not anybetter than the numDeriv version
            # "realp2"={hesse<-pracma::hessian_csd(f=fitfunc,x0=cobra$xbest+runif(cobra$dimension,min=-0.001,max=0.001))}, #deprecated because the results are not anybetter than the numDeriv version
             ,"surrogate"={hesse<-numDeriv::hessian(sfunc,x=cobra$xbest)}
             #,"cma"={hesse<-mycma_es(fn=fitfunc,par=cobra$xbest,lower = cobra$lower, upper=cobra$upper,
             #                       control=list(maxit=4*cobra$dimension,lambda=cobra$dimension+1,mu=cobra$dimension,sigma=1/nrow(cobra$A)))},
             #,"1+1cma"={hesse<-optim11(fn=fitfunc,xinit=cobra$xbest,budget=4*cobra$dimension^2+4*cobra$dimension)}, #this option is deprecated and the corresponding code is moved to the deprectaed folder
             #,"pca"={popc<-sortPop(cobra$A,cobra$Fres)$pop;EIG<-eigen(cov(popc[c(1:2*cobra$dimension),]));
             #             coefs<-EIG$values;
             #             coefs[coefs<1e-14]<-1e-14;
             #             coefs<-coefs/max(coefs);
             #             hesse<-sapply(c(1:ncol(EIG$vectors)),function(i)EIG$vectors[,i]*coefs[i]^-1) #still in testing phase and not available for SACOBRA 1.1 version
             )
        #calculating hesse^(-1/2)
        cobra$hesse<-hesse
        SVD<-svd(hesse)
        eps<-1e-25
        invD <- 1/sqrt(SVD$d)
        invD[abs(SVD$d/SVD$d[1])<eps] <- 0
        
        #M<-SVD$v%*%diag(invD)%*%t(SVD$v)
        M<-diag(invD)%*%t(SVD$v)
        cobra$TM<-M
      }
      
      
      #it is checked if the cobra$TM suggests a valid transformation
      #If the transformation is not valid then the VALIDTRANSORFMATION turns False and the transformation is cancelled
      if(all(is.na(cobra$TM))||all(cobra$TM==0.0)){
        cobra$TM<-diag(cobra$dimension)
        print(sprintf("no valid transformation is suggested in iteration %d",nrow(cobra$A)))
        VALIDTRANSORFMATION<-F
      }
      
    
    #-----only debug purposes---------#
    # realF is only used for diagnosis purposes and the real function calls are not stored in the archive    
     realF<-function(x){
       y<-cobra$fnNOarchive(x)[1]
       return(y)
     }   
     #the visHess function plots the hessian matrix error
     #deprecated
     #if(cobra$CA$visHessianError) visHess(OHesse=hesse,cobra)
    #-----only debug purposes---------#
     
     
     #Checking all the evaluated points to see if we can find a better solution
     Xs<-get("ARCHIVE",envir=intern.archive.env)
     Ys<-get("ARCHIVEY",envir=intern.archive.env)
     bestInd<-which(Ys==min(Ys))[1]
     XBEST<-Xs[bestInd,]
     myGrad<-XBEST-cobra$xbest
     #updating the best solution
     cobra$xbest<-XBEST
     cobra$xbest<- pmax(cobra$xbest,cobra$lower)   
     cobra$xbest<- pmin(cobra$xbest,cobra$upper)
     cobra$A[cobra$ibest,]<-cobra$xbest
     cobra$Fres[cobra$ibest]<-Ys[bestInd]
     cobra$fbest<-Ys[bestInd]
     
     # re-evaluating the current points after applying the M transformation
     # note that this step costs extra real function evaluations if cobra$CA$reEval=T
     #if(cobra$CA$reEval)
       AtoEval<-cobra$A
     #else
     #  AtoEval<-cobra$iA
     
     alpha<-cobra$CA$alpha
     tCenter<-cobra$xbest+alpha*myGrad #transformation center
     #popc<-sortPop(cobra$A,cobra$Fres)$pop
     #tCenter<-apply(popc[c(1:2*cobra$dimension),],2,mean)+alpha*myGrad
    # centerChange<-sqrt(sum((tCenter-cobra$tCenter)^2))
     
     if(VALIDTRANSORFMATION ){
            cobra$tCenter<-tCenter
            newA<-AtoEval-t(replicate(nrow(AtoEval),tCenter))
            newA<-newA%*%cobra$TM
            newA<-newA+t(replicate(nrow(AtoEval),tCenter))
            cobra$TA<-AtoEval
            cobra$TFres<-apply(newA,1,fitfunc) 
      
        #cobra$TA<-sel$AtoEval
        #cobra$TFres<-sel$TFres
       cobra <- trainSurrogates(cobra);# side effect: cobra$constraintSurrogates, cobra$fitnessSurrogate
       #evalModel(nPoints=10000,model=cobra$fitnessSurrogate,dimension=cobra$dimension,fn=realF,Tr=cobra$TM,Tc=cobra$tCenter) 
       #This function is deprecated and the R script is moved to the deprecated folder
     }else{ #if the validtransformation!=T then we turn off the transfomarion flag and we generate normal surrogates
       cobra <- trainSurrogates(cobra);# side effect: cobra$constraintSurrogates, cobra$fitnessSurrogate
         
     }


    }#end of model generation by means of whitening
    
    
    ##########################################################
    # STEP6.2: Determine Distance requirement                #
    ##########################################################
    #why -nRepair? - globalOptCounter counts only real optimization steps
    #gama<-cobra$XI[((globalOptCounter-nRepair) %% length(cobra$XI))+1]  #/SB/ using a counter only for the global optimization excluding repair and trust region
    gama<-cobra$XI[((globalOptCounter) %% length(cobra$XI))+1]  #/SB/ using a counter only for the global optimization excluding repair and trust region
    cobra$ro <- gama*cobra$l
    cobra$EPS<- EPS
    
    ##########################################################
    # STEP6.3: Optimize on surrogates                        #
    ##########################################################
    ptm <- proc.time()
    subMin <- list()
    # print_level = 2  # { 0 | 2 } print no or more diagnositc information --> only avail for nloptr, not for wrapper cobyla 
    verboseprint(cobra$verbose, important=FALSE,paste(cobra$seqOptimizer, " optimization on surrogate ..."))
    
    #### /SB/-Random Start algorithm
    if(cobra$sac$RS){
      cobra<-RandomStart(cobra)
      xStart<-cobra$xStart
      #if(any(cobra$xbest!=xStart)) cobra$noProgressCount<-0
    }else{
      xStart<-cobra$xbest
    }
   
    testit::assert("",all(xStart>=cobra$lower));
    testit::assert("",all(xStart<=cobra$upper));
    cobra <- checkIfCobraOptimizable(cobra);
   # browser()
    
  #cat("Starting optimizer ...\n");
  switch(cobra$seqOptimizer,
           
           #RANDOMLHS={subMin<-RS(fun=subProb,lb=cobra$lower, ub=cobra$upper, n=cobra$seqFeval)},
           COBYLA={ subMin<-nloptr::cobyla(xStart,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA,control=list(maxeval=cobra$seqFeval,xtol_rel=cobra$seqTol), cobra=cobra); subMin$feval=subMin$iter },
           ISRES ={ subMin<-isres2(xStart,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA, maxeval=cobra$seqFeval, cobra=cobra); subMin$feval=subMin$iter},
           HJKB = { subMin<-dfoptim::hjkb(xStart,fn=subProb,lower=cobra$lower,upper=cobra$upper,control=list(maxfeval=cobra$seqFeval), cobra=cobra) },
           NMKB = { subMin<-nmkb2(xStart,fn=subProb,lower=cobra$lower,upper=cobra$upper, control=list(maxfeval=cobra$seqFeval,tol=cobra$seqTol), cobra=cobra) },
           #ACTIVECMA  = { subMin <- ActiveOnePlusOneCMAES(xStart, subProb, length(cobra$xbest), opts=list(esname="ActiveCMAES", lb=cobra$lower, ub=cobra$upper, 
           #                maxFunEvals=cobra$seqFeval, 
           #                mu=cobra$seqMu, 
           #                lambda=cobra$seqLambda, 
           #                sigma=cobra$seqStepSize)); 
           #                subMin$convergence <- 1;
           #},
           #RANDOMSEARCH = { subMin <- randomSearch(cobra$xbest, fn=subProb2, lower=cobra$lower,upper=cobra$upper,control=list(maxfeval=cobra$seqFeval, sd=0.05))},
           ISRESCOBY = { subMin<-isresCobyla(xStart,fn=subProb2,hin=gCOBRA, cobra=cobra); subMin$feval=subMin$iter; cobra$nCobyla=subMin$nCobyla; cobra$nIsres=subMin$nIsres;},
           DEOPTIM={subMin<-DEoptim::DEoptim(fn=subProb,lower=cobra$lower,upper=cobra$upper,DEoptim::DEoptim.control(itermax=1000,trace=Inf,initialpop=xStart+matrix(runif(10*length(xStart)^2),nrow=10*length(xStart),ncol=length(xStart))), cobra=cobra); subMin$feval=subMin$optim$nfeval ; subMin$par=subMin$optim$bestmem;subMin$convergence=NA}
    )
    optimTime <- (proc.time()-ptm)[3]
    verboseprint(cobra$verbose, important=FALSE,paste(" finished (",subMin$feval,"iterations,",optimTime,"sec )"))
    #cat("Optimization time for", subMin$feval, " iterations:", optimTime, "seconds\n")  
    #cat("Predicted infill value:",subMin$value,"\n")
    
    if (cobra$DEBUG_XI) {
      # If cobra$DEBUG_XI==TRUE, then print xStart in every iteration 
      # and add later some extra debug info to cobra$df (columns fxStart, fxbest, errY, XI, RS)
      cat("** xStart =",xStart," **\n")
      cobra$xStart = xStart
      cobra$DEBUG_RS = (!all(xStart==cobra$xbest)) 
    }
    
    if (penaF[1]*penaF[2] < penaF[3])
      penaF[1] <- penaF[1]*penaF[2]
    
    cobra$penaF <- penaF
    cobra$sigmaD <- sigmaD <- checkDistanceReq(subMin,cobra$fitnessSurrogate,cobra$ro,sigmaD,CHECKDIST)    
    subMin$par<- pmax(subMin$par,cobra$lower)   
    subMin$par<- pmin(subMin$par,cobra$upper)
    #subMin$par[subMin$par==1]<-runif(length(which(subMin$par==1)),min=-0.9,max=0.9)
    #subMin$par[subMin$par==-1]<-runif(length(which(subMin$par==-1)),min=-0.9,max=0.9)
    
    xNew<-subMin$par
    attr(ev1,"state") <- "optimized"
    
    ##########################################################
    # STEP6.4: Evaluate real functions                       #
    ##########################################################
    #cobra$noProgressCount<-cobra$noProgressCount+1
    # cobra$noProgressCount is increased unconditionally here, it will be later reset to 0, if it 
    # turns out that the current iteration was successful (produced a new xBest). If there are
    # unsuccessful iterations in a row, then cobra$noProgressCount counts them.
    globalOptCounter<-globalOptCounter+1# this is a counter which counts all main iterates without repair or tr
    cobra$fe<-cobra$fe+1
    #
    #evaluate xNew on the real functions + do refine step (if cobra$equHandle$active) 
    ev1 <- evalReal(cobra,ev1,xNew,subMin$value,subMin$feval,subMin$convergence,optimTime,currentEps)
    cobra$TFres<-c(cobra$TFres,ev1$xNewEvalT[1])
    #calculate the p-Effect
    if((nrow(cobra$A)%%cobra$sac$onlineFreqPLOG==0 && cobra$sac$onlinePLOG) || nrow(cobra$A)==cobra$initDesPoints){
      cobra<-calcPEffect(cobra,ev1$xNew,ev1$xNewEval) 
    }
    
    
    ##########################################################
    # STEP6.5 & 6.6: Update Information and Counters         #
    ##########################################################
    cobra <- updateInfoAndCounters(cobra,ev1$xNew,ev1$xNewEval
                                   ,ev1$newNumViol,ev1$newMaxViol,ev1$trueMaxViol
                                   ,phase,currentEps)
    
    num<-nrow(cobra$A)
    cobra$predC <- ev1$predC
   
    
      
    ##########################################################
    # STEP6.8: Update and save cobra                         #
    ##########################################################
    # This includes the update of cobra$xbest, cobra$fbest, cobra$ibest.
    # In case of cobra$equHandle$active it will call updateCobraEqu in modifyEquCons.R
    # to do the job.
    cobra <- updateSaveCobra(cobra,ev1,subMin,sigmaD,penaF,gama,EPS,
                             fitFuncPenalRBF,distRequirement)
    ev1<-cobra$ev1
    
    # ---- WK-debug only ---
    DBG_F02=F
    if (DBG_F02) {
      #browser()
      d = cobra$dimension
      xmat=seq(-1,1,0.01)
      par(mfrow=c(2,2))
      for (k in 0:min(3,d-1)) {
        getOption("scipen")
        opt <- options("scipen" = -1)
        fmat=sapply(1:length(xmat),function(i){interpRBF(c(rep(-0.2,k),xmat[i],rep(-0.2,d-1-k)),cobra$fitnessSurrogate)})
        tmat=sapply(1:length(xmat),function(i){cobra$fnNOarchive(c(rep(-0.2,k),xmat[i],rep(-0.2,d-1-k)))[1]})
       # fomat=sapply(1:length(xmat),function(i){interpRBF(c(rep(1,k),xmat[i],rep(1,d-1-k)),cobra$fitnessSurrogate)})
       # tomat=sapply(1:length(xmat),function(i){cobra$fnNOarchive(c(rep(1,k),xmat[i],rep(1,d-1-k)))})
       # plot(xmat,fmat,type="l",main=sprintf("dimension=%s",k+1),ylim=c(0,max(c(fmat))),xlab=substitute(paste(x[i]),list(i=k+1)),ylab="f(x)")   # a plot along the x-axis, the narrow, flat 'valley'
        myylab<-substitute(paste("f(",x[i],")"),list(i=k+1))
        dim<-c(1:4)
        dim<-dim[-(k+1)]
        myylab<-substitute(paste("f(",x[i],", ",x[dim1],"=-1",", ",x[dim2],"=-1",", ",x[dim3],"=-1)"),list(i=k+1,dim1=dim[1],dim2=dim[2],dim3=dim[3]))
        #plot(xmat*5,fmat,type="l",ylim=c(0,max(fmat)),xlab=substitute(paste(x[i]),list(i=k+1)),ylab=myylab,lwd = 5, cex=2,cex.names=1.3,cex.axis=1.3,cex.lab=1.1)   # a plot along the x-axis, the narrow, flat 'valley'
        plot(xmat*5,fmat,type="l",xlab=substitute(paste(x[i]),list(i=k+1)),ylab=myylab,lwd = 5, cex=2,cex.names=1.3,cex.axis=1.3,cex.lab=1.1)   # a plot along the x-axis, the narrow, flat 'valley'
        lines(xmat*5,tmat,col="red",lwd = 2, cex=2)
        #lines(xmat,fomat,col="green")
        #lines(xmat,tomat,col="pink")
        
      }
      par(mfrow=c(1,1))
      legend(1, 5e+6, legend=c("real function", "surrogate"),
             col=c("red", "black"), lwd=c(2,5), cex=1,box.lwd = 0,box.col = "white")
    }
    # ---- WK-debug only ---
    
    ##########################################################
    # STEP6.7: Adjust Margins                                #
    ##########################################################
    adjustMargins(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,cobra$epsilonMax,currentEps)
    
   
    ##########################################################
    # STEP6.9: Repair Infeasible                             #
    ##########################################################
    if(cobra$repairInfeas==TRUE) 
    {
      xNewHasEverBestFitness = (cobra$Fres[length(cobra$numViol)] < cobra$fbest + cobra$ri$marFres) 
      if (!cobra$ri$repairOnlyFresBetter) xNewHasEverBestFitness = TRUE
      # /WK/ if repairOnlyFresBetter=T, repair only iterates with fitness < so-far-best-fitness  
      
      if( (ev1$newNumViol>=1) && 
          (xNewHasEverBestFitness) &&   
          (ev1$newMaxViol < cobra$ri$repairMargin) &&
          (num < cobra$feval) )            # /WK/ bug fix: no repair, if the last iteration would issue a repair
      {
        cobra$important=(num %% cobra$verboseIter ==0)
        #important to print

        # Build surrogate anew, based on current cobra$A, cobra$Gres
        # This is important for accurate constraint surrogates models near current infeasible point
        #SB 
        print("repair is called")
        cobra<-trainSurrogates(cobra) 
        repairInfeasible <- repairInfeasRI2
        if (cobra$ri$RIMODE==3) repairInfeasible <- repairChootinan
        #---this is the normal repairInfeasible call ---
        xNewRepaired<-repairInfeasible(ev1$xNew,ev1$xNewEval[-1], cobra$constraintSurrogates,cobra)
        #---this is the repairInfeasible call with checkit=TRUE (debug, extra printout) ---
        #xNewRepaired<-repairInfeasible(ev1$xNew,xNewEval[-1], cobra$constraintSurrogates,cobra,TRUE)
        #
        # /WK/ bug fix above: changed repairInfeasible(subMin$par,... to repairInfeasible(ev1$xNew,...
        
        nRepair <- nRepair +1 
        if(all(ev1$xNew==xNewRepaired)){
          print("cannot repair infeasible solution")
          attr(ev1,"state") <- "repaired"
        }
        else {
          xNew<-xNewRepaired
          cobra$fe<-cobra$fe+1
          attr(ev1,"state") <- "repairedWithSuccess"
          cobra$RSDONE<-c(cobra$RSDONE, "RP")
          
          ##########################################################
          # STEP R6.4: Evaluate real functions                     #
          ##########################################################
          ev1 <- evalReal(cobra,ev1,xNew,fitFuncPenalRBF(xNew), NA,NA,NA,currentEps)
          
          cobra$radi<-c(cobra$radi,cobra$radi[length(cobra$radi)])
          # /WK/ Is cobra$radi really needed? And even if so, is a prolongation of 
          #      cobra$radi correct within the repair step??
          
          ##########################################################
          # STEP R6.5 & 6.6: Update Information and Counters         #
          ##########################################################
          cobra <- updateInfoAndCounters(cobra,ev1$xNew,ev1$xNewEval
                                         ,ev1$newNumViol,ev1$newMaxViol,ev1$trueMaxViol
                                         ,phase,currentEps)
          
          num<-nrow(cobra$A)
          
          ##########################################################
          # STEP R6.8: Update and save cobra                       #
          ##########################################################
          cobra <- updateSaveCobra(cobra,ev1,subMin,sigmaD,penaF,gama,EPS,
                                   fitFuncPenalRBF,distRequirement)
          
         ##########################################################
         # STEP R6.7: Adjust Margins                              #
         ##########################################################
         adjustMargins(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,cobra$epsilonMax,currentEps)
         
        } # else of 'all(xNew==xNewRepaired)'
      } 
    } # if(cobra$repairInfeas==TRUE)
    #######################################################################
    
    
    
    if(   cobra$TrustRegion  
       && gama==0.0 && (num < cobra$feval)     # /WK/ gama==0.0 is experimental
       && cobra$maxViol[cobra$ibest] < 0.1)
    {   
      center <- switch(cobra$TRlist$center,
                       "xbest" = cobra$xbest,
                       "xnew" = xNew,
                       stop("type of TR center is not valid"))
      cobra<-trustRegion(cobra,center=center)   # side effect: cobra$fitnessSurrogateTR
      xNewTrustreg<-cobra$trustregX

      if(cobra$TRDONE){
        #TRcounter<-TRcounter+1 # /WK/ now inside trustRegion
        if(any(cobra$xbest!=xNewTrustreg)){ 
          xNew<-xNewTrustreg
          xNew <- pmax(xNew,cobra$lower)   
          xNew <- pmin(xNew,cobra$upper)  
          cobra$fe<-cobra$fe+1
          attr(ev1,"state") <- "TRWithSuccess"
          
          ##########################################################
          # STEP TRA6.4: Evaluate real functions                     #
          ##########################################################
          ev1 <- evalReal(cobra,ev1,xNew,fitFuncPenalRBF(xNew),NA,NA,NA,currentEps
                          ,fitnessSurrogate=cobra$fitnessSurrogateTR)

          cobra$radi<-c(cobra$radi,cobra$radi[length(cobra$radi)])
          newPredY = tail(ev1$predY,1)
          if(cobra$DEBUG_TR){cat("[",nrow(cobra$A),"]TR DEBUG ::approx error=",abs(ev1$xNewEval[1]-newPredY),"\n")
                             cobra$TRapprox<-c(cobra$TRapprox,abs(ev1$xNewEval[1]-newPredY)) 
                             cobra$TRiter<-c(cobra$TRiter,nrow(cobra$A))}
          
          ##########################################################
          # STEP TRA6.5 & 6.6: Update Information and Counters         #
          ##########################################################
          cobra <- updateInfoAndCounters(cobra,ev1$xNew,ev1$xNewEval
                                         ,ev1$newNumViol,ev1$newMaxViol,ev1$trueMaxViol
                                         ,phase,currentEps)
          if(cobra$DEBUG_TR && all(cobra$xbest==ev1$xNew)){
            cat("[",nrow(cobra$A),"]TR DEBUG ::TR improved the best objective"
                ,cobra$fbestArray[length(cobra$fbestArray)-1]-cobra$fbestArray[length(cobra$fbestArray)],"\n")
          }
          
          num<-nrow(cobra$A)
          ##########################################################
          # STEP TRA6.8: Update and save cobra                     #
          ##########################################################
          cobra <- updateSaveCobra(cobra,ev1,subMin,sigmaD,penaF,gama,EPS,
                                   fitFuncPenalRBF,distRequirement,
                                   fitnessSurrogate=cobra$fitnessSurrogateTR)
          ##########################################################
          # STEP TRA6.7: Adjust Margins                            #
          ##########################################################
          adjustMargins(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,cobra$epsilonMax,currentEps)
          
        }else{
          warning("trust region spotted the old point")
          attr(ev1,"state") <- "TR"
        } # if-else(any(cobra$xbest!=xNewTrustreg))
        
      } # if(cobra$TRDONE)
      
     #cobra<-adaptRadi(cobra)
    } # if(cobra$TrustRegion...) 
    

    #17.10.2017 SB only for diagnosis purposes
    #if(cobra$LocalExpensiveSearch && any(cobra$LES$Iter==num)){
     # cobra<-LocalExpensiveSearch(cobra)
   # } 
    
 } # while(num)
  
  #
  # NEW /WK/2016/01: We invalidate (set to NA) all iterates in cobra$df$Best before the 1st feasible point
  #
  indFeas <- which(cobra$df$feasible) 
  if (length(indFeas)==0) {
    warning("No feasible solutions found in cobraPhaseII!");    
  } else {
    # invalidate iterates before the 1st feasible point
    cobra$df$Best[1:(indFeas[1]-1)] <- NA
  }
  
  cobra$optimizationTime <- ev1$optimizationTime;
  cobra$predC <- ev1$predC

  cobra$adjustMargins<-adjustMargins; # this is to make adjustMargins (and all functions in its 
                                    # environment = all inner functions of cobraPhaseII)
                                    # accessible from the outside
                                    # (e.g. via 
                                    #    environment(cobra$adjustMargins)$fitFuncPenalRBF 
                                    # )

  # reset some 'inner' variables of cobra so that the list returned is simpler:
  if (!cobra$saveSurrogates) {
    cobra$constraintSurrogates <- NULL;
    cobra$fitnessSurrogate <- NULL
    cobra$constraintSurrogatesTR <- NULL;
    cobra$fitnessSurrogateTR <- NULL
  }
  cobra$fitnessSurrogate1 <- NULL
  cobra$fitnessSurrogate2 <- NULL
  cobra$printP <- NULL

  return(cobra)
}
