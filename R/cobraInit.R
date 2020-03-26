#
#Samineh Bagheri, Patrick Koch, Wolfgang Konen
#Cologne University of Applied Sciences
#
#April,2014 - May,2015
#cobraInitial.R
#

#DEBUG_RBF=list(active=FALSE,overlayTrueZ=FALSE,DO_SNAPSHOT=FALSE,every=2)

######################################################################################
# cobraInit
#
#' Initial phase for SACOBRA optimizer
#'
#' In this phase the important parameters are set and the initial design population are evaluated on the real function. The problem to solve is: 
#' \deqn{ \mbox{Minimize}\quad  f(\vec{x}) , \vec{x} \in [\vec{a},\vec{b}] \subset \mathbf{R}^d }
#' \deqn{ \mbox{subject to}\quad g_i(\vec{x}) \le 0, i=1,\ldots,m    }
#' \deqn{ \mbox{~~~~~~~~~~}\quad\quad h_j(\vec{x}) = 0, j=1,\ldots,r.   }
#'
#Detail:
#' If \code{epsilonInit} or \code{epsilonMax} are NULL on input, then \code{cobra$epsilonInit}
#' and \code{cobra$epsilonMax},  resp., are set to \code{0.005*l} where \code{l} is the smallest 
#' side of the search box.
#' 
#' Note that the parameters \code{penaF}, \code{sigmaD}, \code{constraintHandling} are only 
#' relevant for penalty-based internal optimizers \code{\link[dfoptim]{nmkb}} or HJKB. They are NOT relevant  
#' for default optimizer \code{\link[nloptr]{cobyla}}.
#' 
#' Although the software was originally designed to handle only constrained optimization problems, 
#' it can also address unconstrained optimization problems  
#' 
#' How to code which constraint is equality constraint? - Function \code{fn} should return 
#' an \eqn{(1+m+r)}-dimensional vector with named elements. The first element is the objective, the 
#' other elements are the constraints. All equality constraints should carry the name \code{equ}. 
#' (Yes, it is possible that multiple elements of a vector have the same name.) 
#'
#' @param xStart        a vector of dimension \code{d} containing the starting point for the optimization problem
#' @param fn            objective and constraint functions: \code{fn} is a function accepting 
#'                      a \code{d}-dimensional vector \eqn{\vec{x}} and returning an \eqn{(1+m+r)}-dimensional
#'                      vector \code{c(}\eqn{f,g_1,\ldots,g_m,h_1,\ldots,h_r}\code{)}  
#' @param fName         the results of \code{\link{cobraPhaseII}} are saved to \code{<fname>.Rdata}
#' @param lower         lower bound \eqn{\vec{a}} of search space, same dimension as \code{xStart}
#' @param upper         upper bound \eqn{\vec{b}} of search space, same dimension as \code{xStart}
#' @param feval         maximum number of function evaluations
#' @param initDesign    ["LHS"] one out of ["RANDOM","LHS","BIASED","OPTIMIZED","OPTBIASED"]
#' @param initDesPoints [\code{2*d+1}] number of initial points, must be smaller than \code{feval}
#' @param initDesOptP   [NULL] only for initDesign=="OPTBIASED": number of points for the "OPT"
#'                        phase. If NULL, take initDesPoints.
#' @param initBias      [0.005] bias for normal distribution in "OPTBIASED" and "BIASED"
#' @param skipPhaseI    [TRUE] if TRUE, then skip \code{\link{cobraPhaseI}}
#' @param seqOptimizer  ["COBYLA"] string defining the optimization method for COBRA phases I 
#'                        and II, one out of ["COBYLA","ISRES","HJKB","NMKB","ISRESCOBY"]
#' @param seqFeval      [1000] maximum number of function evaluations on the surrogate model
#' @param seqTol        [1e-6] convergence tolerance for sequential optimizer, see param \code{tol} 
#'                      in \code{\link[dfoptim]{nmkb}} or param \code{control$xtol_rel} 
#'                      in \code{\link[nloptr]{cobyla}}
#' @param ptail         [TRUE] TRUE: with, FALSE: without polynomial tail in \code{trainRBF}
#' @param squares      [TRUE] set to TRUE for including the second order polynomials in building the fitness and constraint surrogates in \code{trainRBF}
#' @param XI            [DRCL] magic parameters for the distance requirement cycle (DRC)
#' @param epsilonInit   [NULL] initial constant added to each constraint to maintain a certain margin to boundary
#' @param epsilonMax    [NULL] maximum for constant added to each constraint 
#' @param cobraSeed     [42] seed for random number generator
#' @param conTol        [0.0] constraint violation tolerance
#' @param repairInfeas  [FALSE] if TRUE, trigger the repair of appropriate infeasible solutions
#' @param ri            [\code{\link{defaultRI}()}] list with other parameters for 
#'                      \code{\link{repairInfeasRI2}}
#' @param saveSurrogates [FALSE] if TRUE, then \code{\link{cobraPhaseII}} returns the last surrogate models in
#'                      cobra$fitnessSurrogate and cobra$constraintSurrogates
#' @param saveIntermediate [FALSE] if TRUE, then \code{\link{cobraPhaseII}} saves intermediate results
#'                      in dir 'results/' (create it, if necessary)
#' @param RBFmodel      ["cubic"] a string for the type of the RBF model, "cubic", "Gaussian" or "MQ"
#' @param RBFwidth      [-1] only relevant for Gaussian RBF model. Determines the width \eqn{\sigma}. 
#'                      For more details see parameter \code{width} in \code{\link{trainGaussRBF}} in \code{RBFinter.R}. 
#' @param widthFactor   [1.0] only relevant for Gaussian RBF model. Additional constant 
#'                      factor applied to each width \eqn{\sigma} 
#' @param GaussRule     ["One"] only relevant for Gaussian RBF model, see \code{\link{trainGaussRBF}}                                     
#' @param RBFrho        [0.0] experimental: 0: interpolating, > 0, approximating (spline-like) Gaussian RBFs
#' @param trueFuncForSurrogates  [FALSE] if TRUE, use the true (constraint & fitness) functions
#'                      instead of surrogates (only for debug analysis)
#' @param equHandle     [\code{\link{defaultEquMu}()}] list with of parameters for            
#'                      equality constraint handling described in \code{\link{defaultEquMu}()}. equHandle$active is set to TRUE by default.
#' @param rescale       [TRUE] if TRUE, transform the input space from \code{[lower,upper]} 
#'                      to hypercube \code{[newlower,newupper]^d}
#' @param newlower      [-1] lower bound of each rescaled input space dimension, if \code{rescale==TRUE}
#' @param newupper      [+1] upper bound of each rescaled input space dimension, if \code{rescale==TRUE}
#' @param solu          [NULL] the best-known solution (only for diagnostics). This is normally a 
#'                      vector of length d. If there are multiple solutions, it is a matrix with d
#'                      columns (each row is a solution). If NULL, then the current best point
#'                      will be used in \code{\link{cobraPhaseII}}. 
#'                      \code{solu} is given in original input space.
#' @param TrustRegion   [FALSE] if TRUE, perform trust region algorithm \code{\link{trustRegion}}. 
#' @param TRlist        [\code{\link{defaultTR}()}] a list of parameters, needed only 
#'                      in case \code{TrustRegion==TRUE}.  
#' @param sac           [\code{\link{defaultSAC}(DOSAC)}] list with other parameters for SACOBRA.  
#' @param DOSAC         [1] set one out of [0|1|2]. \cr
#'                      0: COBRA-R settings, \cr 1: SACOBRA settings, \cr 2: SACOBRA settings with fewer parameters. \cr 
#'                      The precise settings are documented in \code{\link{defaultSAC}}.
#' @param penaF         [c(3,1.7,3e5)] parameters for dynamic penalty factor (fct subProb in 
#'                      \code{\link{cobraPhaseII}}): \code{c(start,augment,max)}, only relevant \code{if seqOptimizer==HJKB} or \code{seqOptimizer==NMKB}
#' @param sigmaD        [c(3,2.0,100)] parameters for dynamic distance factor (fct subProb in 
#'                      \code{\link{cobraPhaseII}}): \code{c(start,augment,max)}, , only relevant \code{if seqOptimizer==HJKB} or \code{seqOptimizer==NMKB}
#' @param constraintHandling ["DEFAULT"] (other choices: "JOINESHOUCK", "SMITHTATE", "COIT", "BAECKKHURI";
#'                      experimental, only relevant \code{if seqOptimizer==HJKB} or \code{seqOptimizer==NMKB}
#'                      see the code in function \code{subProb} in \code{\link{cobraPhaseII}})          
#' @param MS            [\code{\link{defaultMS}()}] list of online model selection parameters described in \code{\link{defaultMS}}. 
#' If \code{MS$active = TRUE} then the type of RBF models for each function will be selected automatically and the \code{RBFmodel} parameter becomes irrelevant.
#' @param DEBUG_XI      [FALSE] if TRUE, then print in \code{\link{cobraPhaseII}} extra debug information: 
#'                      \code{xStart} in every iteration to console and add some extra debug 
#'                      columns to \code{cobra$df}
#' @param DEBUG_RBF     [\code{defaultDebugRBF()}] list with settings for visualization RBF (only for \code{d==2}) 
#' @param DEBUG_TRU     [FALSE] visualize trust-region RBF (only for dimension==2) 
#' @param DEBUG_TR      [FALSE] prints information about trust region status and visualisation for \code{d==2} (coming soon)
#' @param DEBUG_RS      [FALSE] prints the RS probability in each iteration in the console
#' @param verbose       [1] set one out of [0|1|2], how much output to print
#' @param verboseIter   [10] an interegr value. Printing the summarized results after each \code{verboseIter} iterations.
#' @param conditioningAnalysis [\code{\link{defaultCA}()}] A list with setting for the objective function conditioning analysis and online whitening
#' @return \code{cobra}, an object of class COBRA, this is a (long) list containing most
#'   of the argument settings (see above) and in addition (among others):
#'      \item{\code{A}}{ (feval x dim)-matrix containing the initial design points in input .  
#'            space. If rescale==TRUE, all points are in  \strong{rescaled} input space. }
#'      \item{\code{Fres}}{ a vector of the objective values of the initial design points }
#'      \item{\code{Gres}}{ a matrix of the constraint values of the initial design points }
#'      \item{\code{nConstraints}}{ the total number \eqn{m+r} of constraints } 
#'      \item{\code{Tfeas}}{ the threshhold parameter for the number of consecutive iterations 
#'            that yield feasible solutions before margin epsilon is reduced }
#'      \item{\code{Tinfeas}}{ the threshhold parameter for the number of consecutive iterations 
#'            that yield infeasible solutions before margin epsilon is increased }
#'      \item{\code{numViol}}{ number of constraint violations }
#'      \item{\code{maxViol}}{ maximum constraint violation}
#'      \item{\code{trueMaxViol}}{ maximum constraint violation}
#'      \item{\code{trustregX}}{A vector of all refined solutions generated by trust region algorithm (see \code{trustRegion})}
#'      
#'   
#'   Note that \code{cobra$Fres}, \code{cobra$fbest}, \code{cobra$fbestArray} and similar contain 
#'   always the objective values of the orignial function \code{cobra$fn[1]}. (The surrogate models 
#'   may be trained on a \code{\link{plog}}-transformed version of this function.)
#' 
#' @examples 
#' ## Initialize cobra. The problem to solve is the sphere function sum(x^2)    
#' ## with the equality constraint that the solution is on a circle with 
#' ## radius 2 and center at c(1,0).
#' d=2
#' fName="onCircle"
#' cobra <- cobraInit(xStart=rep(5,d), fName=fName,
#'                    fn=function(x){c(obj=sum(x^2),equ=(x[1]-1)^2+(x[2]-0)^2-4)},  
#'                    lower=rep(-10,d), upper=rep(10,d), feval=40)
#'                    
#' ## Run sacobra optimizer
#' cobra <- cobraPhaseII(cobra)
#' 
#' ## The true solution is at solu = c(-1,0) (the point on the circle closest 
#' ## to the origin) where the true optimum is fn(solu)[1] = optim = 1
#' ## The solution found by SACOBRA:
#' print(getXbest(cobra))
#' print(getFbest(cobra))
#' 
#' ## Plot the resulting error (best-so-far feasible optimizer result - true optimum)
#' ## on a logarithmic scale:
#' optim = 1
#' plot(abs(cobra$df$Best-optim),log="y",type="l",ylab="error",xlab="iteration",main=fName)
#' 
#' @seealso   \code{\link{startCobra}}, \code{\link{cobraPhaseI}}, \code{\link{cobraPhaseII}}
#' @author Wolfgang Konen, Samineh Bagheri, Patrick Koch, Cologne University of Applied Sciences
#' @export
#' 
#' 
#' 
######################################################################################
cobraInit <- function(xStart, fn, fName, lower, upper, feval, 
                      initDesign="LHS",
                      initDesPoints=2*length(xStart)+1, initDesOptP=NULL, initBias=0.005,
                      skipPhaseI=TRUE,
                      seqOptimizer="COBYLA", seqFeval=1000, seqTol=1e-6, 
                      ptail=TRUE, squares=TRUE, conTol=0.0,
                      DOSAC=1, sac=defaultSAC(DOSAC), 
                      repairInfeas=FALSE, ri=defaultRI(),
                      RBFmodel="cubic", RBFwidth=-1,GaussRule="One",widthFactor=1.0,RBFrho=0.0,MS=defaultMS(),
                      equHandle=defaultEquMu(),
                      rescale=TRUE,newlower=-1,newupper=1,  
                      XI=DRCL, 
                      TrustRegion=FALSE,TRlist=defaultTR(),
                      conditioningAnalysis=defaultCA(),
                      penaF=c(3.0, 1.7, 3e5), sigmaD=c(3.0,2.0,100), constraintHandling="DEFAULT",
                      verbose=1,verboseIter=10,
                      DEBUG_RBF=defaultDebugRBF(), DEBUG_TR=FALSE, 
                      DEBUG_TRU=FALSE, DEBUG_RS=FALSE, DEBUG_XI=FALSE, 
                      trueFuncForSurrogates=FALSE, 
                      saveIntermediate=FALSE, saveSurrogates=FALSE,
                      epsilonInit=NULL, epsilonMax=NULL, solu=NULL,
                      cobraSeed=42 )
{
  #if(length(xStart)>2)vis<-FALSE
  originalfn<- fn
  originalL <- lower
  originalU <- upper
  phase<-"init"

  testit::assert("cobraInit: xStart contains NaNs",all(!is.nan(xStart)));
  testit::assert("cobraInit: lower<upper violated",all(lower<upper));
    
  dimension<-length(xStart)         # number of parameters
  #browser()
  if(rescale){
    lb<-rep(newlower,dimension)
    up<-rep(newupper,dimension)
    xStart<-sapply(1:dimension , function(i){scales::rescale(xStart[i],to=c(lb[i],up[i]),from=c(originalL[i],originalU[i]))
    })
    if(!is.null(solu))solu<-sapply(1:dimension , function(i){scales::rescale(solu[i],to=c(lb[i],up[i]),from=c(originalL[i],originalU[i]))
    })
    fn<-rescaleWrapper(fn,originalL,originalU,dimension,newlower,newupper)
    lower<-lb
    upper<-up
  }
  #
  
  
  testit::assert("cobraInit: Too many init design points", initDesPoints<feval)
  
  CONSTRAINED=T
  xStartEval<-fn(xStart)
  nConstraints<-length(xStartEval)-1
  if(nConstraints==0)CONSTRAINED<-F
  assert("This version does not support conditioning analysis for constrained problems ",!CONSTRAINED || !conditioningAnalysis$active)
 if(!CONSTRAINED)testit::assert("cobraInit: This version does not support trust Region functionality for unconstrained Problems", !TrustRegion )
  
  # old version
  # testit::assert("cobraInit: There should be at least one constraint!",nConstraints>0)
  # NEW[27.09.2017] the software accepts unconstarined problems
  testit::assert("cobraInit: nConstraints cannot be smaller than 0",nConstraints>=0)

  if(!CONSTRAINED)verboseprint(verbose=verbose, important=TRUE,"An unconstrained problem is being addressed")
  
  
  l <- min(upper - lower) # length of smallest side of search space
  if (is.null(epsilonInit)) epsilonInit<- 0.005*l
  if (is.null(epsilonMax))  epsilonMax <- 2*0.005*l
  if (is.null(initDesOptP)) initDesOptP <- initDesPoints
  

  
  set.seed(cobraSeed)
  dimension<-length(xStart)         # number of parameters
  iteration<-0  
  # We just have one starting point
  xStart <- xStart
  
  I<-matrix()
  Gres <- matrix()
  Fres <- c()
  randomInitialPoints <- 10 # additional random points, required for optimized design to avoid crashing of SVD in the next phases
  
  #adding archiving function to fn

  assign('ARCHIVE', NULL, envir=intern.archive.env)
  assign('ARCHIVEY', NULL, envir=intern.archive.env)
  fnNOarchive<-fn
  fn<-function(x){

    y<-fnNOarchive(x)
    ARCHIVE<-rbind(get("ARCHIVE",envir=intern.archive.env),as.numeric(x))
    assign("ARCHIVE", ARCHIVE, envir = intern.archive.env)

    ARCHIVEY<-rbind(get("ARCHIVEY",envir=intern.archive.env),as.numeric(y[1]))
    assign("ARCHIVEY", ARCHIVEY, envir = intern.archive.env)
    return(y)
  }
 
  # helper for switch(initDesign): evaluate random solutions in I with fn
  randomResultsFactory <- function(I,fn,dimension) {
    sapply(1:nrow(I), function(i){
      verboseprint(verbose=0,important=FALSE,sprintf("iteration %03d: ",i))
      x<-I[i,]    
      res<-fn(x)
      return(res)
    })
  }
  
  # helper for switch(initDesign): clip all entries in I[,j] to be between lower[j] and upper[j] 
  clipLowerUpper <- function(I,dimension) {
    I <- t(unlist(sapply(1:nrow(I), FUN=function(i)pmax(I[i,],lower))))
    I <- t(unlist(sapply(1:nrow(I), FUN=function(i)pmin(I[i,],upper))))
    return(I)    
  }
  
  cat("\n")
  sw=switch(initDesign,
         "RANDOM" = {
           set.seed(cobraSeed)
        
           I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)stats::runif(dimension,lower,upper))))    
                     
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)
           
           randomResults<-randomResultsFactory(I,fn,dimension)
           # update structures for random solutions
           # update structures for random solutions
           #WK: In order to adapt the code to address unconstraint problems
           if(nConstraints==0){
             Gres<-NULL
             Fres <- as.vector(randomResults)                                                          
             
           }else{
             Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint
             colnames(Gres)<-rownames(randomResults)[-1]
             Fres <- as.vector(randomResults[1,])                                                          
           }
           names(Fres) <- NULL
           colnames(I)=NULL
           rownames(I)=NULL
           "ok"
         },
         
         "LHS" = {      # latin hypercube sampling + xStart
           set.seed(cobraSeed)
           I <- lhs::randomLHS(initDesPoints-1,dimension)    # LHS with values uniformly in [0,1]
           #I <- lhs::randomLHS(initDesPoints,dimension)    # LHS with values uniformly in [0,1]
           for (k in 1:ncol(I)) {
             I[,k] <- lower[k] + I[,k]*(upper[k]-lower[k])
           }
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)

           randomResults<-randomResultsFactory(I,fn,dimension)
           #browser()
           # update structures for random solutions
           #SB: In order to adapt the code to address unconstraint problems
           if(nConstraints==0){
             Gres<-NULL
             Fres <- as.vector(randomResults)                                                          
             
           }else{
             Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint
             colnames(Gres)<-rownames(randomResults)[-1]
             Fres <- as.vector(randomResults[1,])                                                          
           }
           names(Fres) <- NULL
           colnames(I)=NULL
           rownames(I)=NULL
           "ok"
         },
         
         "BIASED" = {  # biased initial population based on starting point
           I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)rnorm(dimension,as.vector(xStart),initBias))))
           names(xStart)=NULL
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)
           
           randomResults<-randomResultsFactory(I,fn,dimension)
           # update structures for random solutions
           #WK: In order to adapt the code to address unconstraint problems
           if(nConstraints==0){
             Gres<-NULL
             Fres <- as.vector(randomResults)                                                          
             
           }else{
             Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint
             colnames(Gres)<-rownames(randomResults)[-1]
             Fres <- as.vector(randomResults[1,])                                                          
           }
           colnames(Gres)<-rownames(randomResults)[-1]
           names(Fres) <- NULL
           "ok"
         },
         
         "BIASEDINF" = {  # infeasible initial population based on starting point
           I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)rnorm(dimension,as.vector(xStart),initBias))))
           names(xStart)=NULL
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)
           
           randomResults<-randomResultsFactory(I,fn,dimension)
           # update structures for random solutions
           #WK: In order to adapt the code to address unconstraint problems
           if(nConstraints==0){
             Gres<-NULL
             Fres <- as.vector(randomResults)                                                          
             
           }else{
             Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint
             colnames(Gres)<-rownames(randomResults)[-1]
             Fres <- as.vector(randomResults[1,])                                                          
           }
           colnames(Gres)<-rownames(randomResults)[-1]
           names(Fres) <- NULL
           #browser()
           infeasibleSolution = sapply(1:initDesPoints, FUN=function(i)any(Gres[i,]>0))
           ind = which(infeasibleSolution)
           I <- I[ind,]
           Gres <- as.matrix(Gres[ind,])
           Fres <- Fres[ind]
           nInfeasible <- nrow(I)
           
           while(nInfeasible < initDesPoints){
             x <- runif(dimension,0,1)
             I <- rbind(I,x)
             nInfeasible <- nInfeasible + 1
             #evaluate random solution on fn
             res<-fn(x) 
             Gres <- rbind(Gres,res[2:length(res)])
             Fres <- c(Fres,res[1])
           }
           Gres <- as.matrix(Gres)
           colnames(I)=NULL
           rownames(I)=NULL
           names(Fres) <- NULL
           "ok"
         },
         
         "OPTIMIZED" =  { # optimized (HJKB) initial points    
           # start Hooke & Jeeves initial search
           initialHookeJeeves <- initialHjkb(xStart,fn=fn,lower=lower,upper=upper,control=list(maxfeval=initDesPoints) )
           I <- as.matrix(initialHookeJeeves$historyData$xArchive)
           duplicatedIndices <- which(duplicated(I))
           if(length(duplicatedIndices)!=0){
             I <- I[-duplicatedIndices,]
           }
           
           I <- clipLowerUpper(I,dimension)
           colnames(I)=NULL
           rownames(I)=NULL
           Gres <- as.matrix(initialHookeJeeves$historyData$constraintArchive)
           if(length(duplicatedIndices)!=0){
           Gres <- Gres[-duplicatedIndices,]
           }
           Fres <- initialHookeJeeves$historyData$yArchive
           if(length(duplicatedIndices)!=0){
             Fres <- Fres[-duplicatedIndices]
           }

           names(Fres) = NULL
           #TODO probably need to reset initial design size, because Hooke&Jeeves exceeds budget:
           initDesPoints <- length(Fres)
           "ok"           
         },
            
         "OPTCOBYLA"=, "OPTBIASED"= {  # COBYLA-optimized initial points + BIASED design
              # start COBYLA initial search
              #
              #resetSoluArchive()    # this was necessary for fnArchive_OLD.R
              #setArchiveFunc(fn)
              fnArchiveF <- fnArchiveFactory(fn);
              initialCobyla <- nloptr::cobyla(xStart,fn=function(x){fnArchiveF(x)[1]},lower=lower,upper=upper
                                      , hin=function(x){-fn(x)[-1]} 
                                      , control=list(maxeval=initDesOptP,xtol_rel = 1e-9))
              I <- (environment(fnArchiveF))$getSoluArchive()
              I <- clipLowerUpper(I,dimension)
              
              
              randomResults<-randomResultsFactory(I,fn,dimension)
              # update structures for random solutions
              #WK: In order to adapt the code to address unconstraint problems
              if(nConstraints==0){
                Gres<-NULL
                Fres <- as.vector(randomResults)                                                          
                
              }else{
                Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint
                colnames(Gres)<-rownames(randomResults)[-1]
                Fres <- as.vector(randomResults[1,])                                                          
              }
              
              if (initDesign=="OPTBIASED") {
                # find the best feasible point (if any, else the best infeasible point) ...
                feas<-apply(Gres,1,function(x)all(x<=0))
                if (any(feas==TRUE)) Fres[!feas] <- Inf
                xStart <- I[which.min(Fres),]
                #
                # ... and construct initDesPoints 'BIASED' design points around it
                I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)rnorm(dimension,as.vector(xStart),initBias))))
                names(xStart)=NULL
                I <- rbind(I,xStart)
                I <- clipLowerUpper(I,dimension)
                randomResults<-randomResultsFactory(I,fn,dimension)
                # update structures for random solutions
                #WK: In order to adapt the code to address unconstraint problems
                if(nConstraints==0){
                  Gres<-NULL
                  Fres <- as.vector(randomResults)                                                          
                  
                }else{
                  Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint
                  colnames(Gres)<-rownames(randomResults)[-1]
                  Fres <- as.vector(randomResults[1,])                                                          
                }
              }

              colnames(Gres)<-rownames(randomResults)[-1]
              names(Fres) <- NULL
              #testit::assert("",initDesPoints==length(Fres))
              #TODO probably need to reset initial design size, because COBYLA varies budget:
              initDesPoints <- length(Fres)
              "ok"           
            },
            
         "InvalidInitDesign"
  ) # switch
  cat("\n")
  testit::assert(sprintf("cobraInit: Wrong value %s for initDesign",initDesign),
                 sw!="InvalidInitDesign")
  # parameters for cycling distance requirement
  #teta<-teta   #distance requirement phaseI
  XI<-XI #c(0.01, 0.001, 0.0005)      #distance requirement phaseII
  
  Tfeas<-floor(2*sqrt(dimension) )  # The threshhold parameter for the number of consecutive iterations that yield feasible solution before the margin is reduced
  Tinfeas<-floor(2*sqrt(dimension)) # The threshold parameter for the number of consecutive iterations that yield infeasible solutions before the margin is increased
  Nmax<-feval
  
  
  fe<-nrow(I) #number of function evaluations
  fbest<-c() # best function value found
  xbest<-c() # best point found # a numeric vector with the length of dimension
  fbestArray<-c()
  xbestArray<-matrix()
    
  ########################################################################################################
  # Update structures                                                                              #
  ########################################################################################################
  #testit::assert("Gres is not a matrix!",is.matrix(Gres) || nConstraints!=0)
 
  #Added by SB 26.03.2020
  if (utils::packageVersion("nloptr") <= "1.2.2.1"){
    options(nloptr.show.inequality.warning = FALSE)
  }
  
  
  A<-I  # A contains all evaluated points
  n<-nrow(A)
  cobra<-list(fn=fn,
              fnNOarchive=fnNOarchive,
            xStart=xbest,
            fName=fName,
            dimension=dimension, 
            nConstraints=nConstraints, 
            lower=lower,
            upper=upper,
            newlower=newlower,
            newupper=newupper,
            originalL=originalL,
            originalU=originalU,
            originalfn=originalfn,
            rescale=rescale,
            feval=feval, 
            A=A, 
            #fbestArray=fbestArray, 
            #xbestArray=xbestArray, 
            #xbest=xbest, 
            #fbest=fbest,
            #ibest=ibest,
            Fres=Fres, 
            Gres=Gres, 
           # numViol=numViol,
           # maxViol=maxViol,
            epsilonInit=epsilonInit, 
            epsilonMax=epsilonMax, 
            XI=XI, 
            #drFactor=drFactor,
            Tinfeas=Tinfeas,
            Tfeas=Tfeas,
            iteration=iteration, 
            initDesPoints=initDesPoints,
            penaF=penaF,
            cobraSeed=cobraSeed,
            conTol=conTol,
            constraintHandling=constraintHandling,
            equHandle=equHandle,
            l=l, 
            repairInfeas=repairInfeas,
            fe=fe,
            ri=ri,
            radi=c(),
            RBFmodel=RBFmodel,
            RBFwidth=RBFwidth,
            RBFrho=RBFrho,
            RULE=GaussRule,
            widthFactor=widthFactor,
            sigmaD=sigmaD,
            ptail=ptail,
            squares=squares,
            saveIntermediate=saveIntermediate,
            saveSurrogates=saveSurrogates,
            skipPhaseI=skipPhaseI,
            solu=solu,
            seqOptimizer=seqOptimizer,
            seqFeval=seqFeval,
            seqTol=seqTol,
            sac=sac,
            trueFuncForSurrogates=trueFuncForSurrogates,
            TrustRegion=TrustRegion,
            TRlist=TRlist,
            TRDONE=FALSE,
            TRind=c(),
            trustregX=c(),
            fCount=0,
            sCount=0,
            DOSAC=DOSAC,
            PLOG=FALSE,
            pShift=0,
            pEffect=NA,
            MS=MS,
            noProgressCount=0,
            DEBUG_XI=DEBUG_XI,
            DEBUG_RBF=DEBUG_RBF,
            DEBUG_TRU=DEBUG_TRU,
            DEBUG_TR=DEBUG_TR,
            DEBUG_RS=DEBUG_RS,
            verbose=verbose,
            verboseIter=verboseIter,
            CA=conditioningAnalysis,
            TFlag=F,
            TM=diag(dimension),
            phase=rep(phase,initDesPoints),
            CONSTRAINED=CONSTRAINED
            )
  cobra$equIndex<-which(colnames(Gres)=="equ")
  ## /SB/ 21.09.2015
  cobra$DEBUG_RBF <- setOpts(cobra$DEBUG_RBF,defaultDebugRBF())
  #cobra$LES <- setOpts(cobra$LES,defaultLES(cobra))
  cobra$CA <- setOpts(cobra$CA,defaultCA())
 if(is.character(cobra$CA$ITER)) if(cobra$CA$ITER=="all")cobra$CA$ITER<-seq(initDesPoints,feval,1)
  
  
  ## /SB/ 12.01.2016 Initializing the model selection parameters
  cobra$MS <- setOpts(cobra$MS,defaultMS())
  ## /WK/ if cobra$ri
  cobra$ri <- setOpts(cobra$ri,defaultRI())
  ## /SB/: extension to SACOBRA
  cobra$TRlist<-setOpts(cobra$TRlist,defaultTR())
  cobra$radi<-rep(cobra$TRlist$radiInit,initDesPoints)
 
 
  if(DOSAC>0) {     
    verboseprint(cobra$verbose, important=FALSE,"Parameter and function adjustment phase")
    cobra$sac<-setOpts(cobra$sac,defaultSAC(DOSAC))
    cobra$pEffect<-cobra$sac$pEffectInit

    if(cobra$sac$aDRC){
      if(length(XI)!=length(DRCL)){
        warning("XI is different from default (DRCL), but sac$aDRC==TRUE, so XI will be set by automatic DRC adjustment!")  
      }else if(any(XI!=DRCL)){
        warning("XI is different from default (DRCL), but sac$aDRC==TRUE, so XI will be set by automatic DRC adjustment!")
      }
      
      verboseprint(cobra$verbose,important=FALSE,"adjusting DRC")
      DRC<-adDRC(max(cobra$Fres),min(cobra$Fres))
      cobra$XI<-DRC
    }
    
    # --- adFit is now called in *each* iteration of cobraPhaseII (adaptive plog) ---
    cobra$GRfact<-1
    cobra$finMarginCoef<-1
    if(cobra$sac$aCF && nConstraints!=0 ){
      verboseprint(cobra$verbose,important=FALSE,"adjusting constraint functions")
      cobra<-adCon(cobra)
      #cobra$fn<-fn  
    }
    cobra$RSDONE<-rep(NA,initDesPoints)
  } # (DOSAC)
  
  
  
  
  
  
  numViol<-sapply(1:initDesPoints , function(i){ # number of initial Constraint violations
    return(sum(cobra$Gres[i,]>0))
  })
  maxViol<-sapply(1:initDesPoints , function(i){
    #y<-max(0,cobra$Gres[i,])
    y<-max(0,cobra$Gres[i,])
    
    return(y)
  })
  trueMaxViol<-maxViol
  
  
  #SB: related to handling equality 02.10.2015
  currentEps<-NA
  equIndex<-which(colnames(Gres)=="equ")
  if(length(equIndex)==0){
    cobra$equHandle$active=FALSE
  }else if(equHandle$active){
    cobra$equIndex<-equIndex
    cobra$equHandle <- setOpts(equHandle,defaultEquMu())
    tempG<-cobra$Gres
    tempG[,cobra$equIndex]<-abs(tempG[,cobra$equIndex])
    trueMaxViol<-sapply(1:initDesPoints , function(i){
      y<-max(0,tempG[i,]*cobra$GRfact)
     # y<-max(0,tempG[i,])
      
      return(y)
    })
    
    switch(cobra$equHandle$initType,
           useGrange={currentEps<-cobra$GrangeEqu},
           TAV={tempG[tempG<0]<-0;currentEps<-median(apply(tempG,1,sum))},
           TMV={currentEps<-median(maxViol)},
           EMV={currentEps<-median(apply(tempG[,cobra$equIndex],1,max))}
    ) 
    currentEps<-max(currentEps,cobra$equHandle$equEpsFinal)
    cobra$currentEps<-rep(currentEps,initDesPoints)
    maxViol<-sapply(1:initDesPoints , function(i){
      y<-max(0,tempG[i,]-currentEps)
      return(y)
    })
    numViol<-sapply(1:initDesPoints , function(i){ # number of initial Constraint violations
      return(sum(tempG[i,]-currentEps>0))
    })
  }
  

  
  cobra$numViol=numViol
  cobra$maxViol=maxViol
  cobra$trueMaxViol=trueMaxViol
  
  
  if(0 %in% numViol){
    fbest<-min(Fres[which(numViol==0)])
    #xbest <- I[which.min(fbest),]        # /WK/2015-09-03/ this was a bug
    xbest <- I[which(Fres==fbest)[1],]
    ibest <- which(Fres==fbest)[1]
  }else{
    # if there is no feasibe point yet, take one from the points with min number of violated constraints:
    minNumIndex<-which(numViol==min(numViol))
    # /WK/ the folllowing lines can select pretty big fbest values which lead to a very bad pShift:
    #index<-minNumIndex[which(maxViol[minNumIndex]==min(maxViol[minNumIndex]))]
    #fbest <- Fres[index[1]]
    #xbest <- A[index[1],]
    #ibest <- index[1]
    # /WK/ alternatve: select among the min-num-viol point the one with smallest Fres
    FresMin <- Fres[minNumIndex]
    ind <- which.min(FresMin)
    index <- minNumIndex[ind]
    
    fbest <- Fres[index[1]]
    xbest <- A[index[1],]
    ibest <- index[1]
  }
  fbestArray<-rep(fbest,initDesPoints)
  xbestArray<-xbest
  
  for(i in c(2:n)){
    xbestArray<-rbind(xbestArray,xbest)
  }
  cobra$fbest<-fbest
  cobra$xbest<-xbest
  cobra$ibest<-ibest
  cobra$fbestArray<-fbestArray
  cobra$xbestArray<-xbestArray
  cat("*** Starting run with seed",cobraSeed,"***\n")
  print("Initialization is done")
  attr(cobra,"state") <- "initialized"
  class(cobra) <- c("COBRA","list")
  return(cobra)
}

verboseprint<-function(verbose, important, message){
  if(verbose!=0){ # if verbose == 0 do not print anything
                  # if verbose == 1 print only important messages 
                  # if verbose == 2 print everything
    if((verbose==2) || ((verbose==1) && (important==TRUE)) ){
      print(message)
    }
  }
}
verbosecat<-function(verbose, message, important=FALSE){
  if(verbose!=0){ # if verbose== 0 do not cat anything
                  # if verbose ==1 cat only important messages 
                  # if verbose ==2 cat everything
    if((verbose==2) || ((verbose==1) && (important==TRUE)) ){
      cat(message)
    }
  }
}


#' 
#' 
#' 
