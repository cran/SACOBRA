#innerFuncs.R
###Wolfgang Konen
###Cologne University of Applied Science
###April/2016

checkIfCobraOptimizable <- function(cobra){
  testit::assert("cobra$ro is not set to a numeric",is.numeric(cobra$ro));
  testit::assert("cobra$EPS is not set to a numeric",is.numeric(cobra$EPS));
  testit::assert("cobra$dimension is not set to a numeric",is.numeric(cobra$dimension));
  testit::assert("cobra$nConstraints is not set to a numeric",is.numeric(cobra$nConstraints));
  testit::assert("cobra$A is not a matrix",is.matrix(cobra$A));
  if(cobra$CONSTRAINED)testit::assert("cobra$constraintSurrogates is not an RBFinter object",inherits(cobra$constraintSurrogates,"RBFinter"));
  testit::assert("cobra$fitnessSurrogate is not an RBFinter object",inherits(cobra$fitnessSurrogate,"RBFinter"));
  if (cobra$trueFuncForSurrogates)
    testit::assert("cobra$fn is not a function",is.function(cobra$fn));
  if (cobra$seqOptimizer=="ISRESCOBYLA") {
    testit::assert("cobra$nCobyla is not set to a numeric",is.numeric(cobra$nCobyla));
    testit::assert("cobra$nIsres is not set to a numeric",is.numeric(cobra$nIsres));
  }
  attr(cobra,"state") <- "optimizable";
  return(cobra)
}

# surrogate penalty function for unconstrained optimization methods  - PHASE I
subProbPhaseI<-function(x,cobra){ 
  if(any(is.nan(x))){
    warning("subProbPhaseI: x value is NaN, returning Inf")
    return(Inf)
  }
  distance <- distLine(x,cobra$A)
  subC<-pmax((cobra$ro-distance),0)
  penalty1<-sum(subC)
  num<-nrow(cobra$A)

  constraintPrediction <- calcConstrPred(x,cobra);
  
  maxViolation <- sapply(1:length(constraintPrediction), FUN=function(i)max(0,constraintPrediction[i]))
  #penalty2 <- sum(maxViolation)
  # Constraint handling: Coit et al. (1996)
  NFT <- rep(0.1,length(maxViolation))
  kappa <- 2
  
  penalty2 <- (length(maxViolation)*sum(maxViolation)*num)*sum(sapply(1:length(maxViolation),FUN=function(i)(maxViolation[i]/NFT[i])^kappa))
  # y<-interpRBF(x, fitnessSurrogate) + (penalty1*cobra$sigmaD[1] + penalty2)*cobra$penaF[1]
  #y<- sum(max(interpRBF(x,constraintSurrogates),0)^2) + (penalty1*cobra$sigmaD[1] + penalty2)*cobra$penaF[1]
  
  #f <- interpRBF(x, cobra$fitnessSurrogate)      # /WK/2016-05/ Bug fix: f is NOT part of subprobPhaseI
  #y <- f  + penalty2 #+ penalty1*cobra$sigmaD[1] #
  y <- penalty2 #+ penalty1*cobra$sigmaD[1]
  return(y)
}

# surrogate evaluation of 'f' for constrained optimization methods  - PHASE I
subProb2PhaseI <- function(x,cobra){

  if(any(is.nan(x))){
    warning("subProb2PhaseI: x value is NaN, returning Inf")
    return(Inf)
  }
  #y<-predict.RBFinter(constraintSurrogates,matrix(x,ncol=dimension))
  
  constraintPrediction <- calcConstrPred(x,cobra);
  
  y<- sum(max(constraintPrediction,0)^2)
  return(y)
}



# surrogate evaluation of '\vec{g}' for constrained optimization methods  - PHASE I
gCOBRAPhaseI <- function(x,cobra) {

  h <- c()
  distance <- distLine(x,cobra$A)
  subC<-pmax((cobra$ro-distance),0)
  h[1] <- sum(subC)

  constraintPrediction <- calcConstrPred(x,cobra);

  h <- (-1.0)*c(h[1], constraintPrediction) # TODO -1* ... is required for COBYLA constraints, maybe also for other optimizers? 
  return(h)
}

# surrogate penalty function for unconstrained optimization methods  - PHASE II
subProb<-function(x,cobra){ 
  #x<-x%*%cobra$TM
  #x<-as.numeric(x)
 # browser()
  
  if(any(is.nan(x))){
    warning("subProb: x value is NaN, returning Inf")
    return(Inf)
  }
  #cat(">>> Entering subProb\n")
  distance <- distLine(x,cobra$A)
  subC<-pmax((cobra$ro-distance),0)
  penalty1<-sum(subC)
  penalty2<-0
  if(cobra$CONSTRAINED){
    constraintPrediction <- calcConstrPred(x,cobra);
    #if (cobra$trueFuncForSurrogates) constraintPrediction <-  cobra$fn(x)[-1]+cobra$EPS^2
    maxViolation <- sapply(1:length(constraintPrediction), FUN=function(i)max(0,constraintPrediction[i]))
    penalty2 <- sum(maxViolation)
  }
  
  y <- subProbConstraintHandling(x,cobra,penalty1,penalty2,maxViolation,cobra$sigmaD,cobra$penaF)  
  if (cobra$trueFuncForSurrogates) y<- cobra$fn(x)[1] + (penalty1*cobra$sigmaD[1]*100 + penalty2)*cobra$penaF[1]
  #cat(">>SubProb: ", interpRBF(x, cobra$fitnessSurrogate), " / ", (C * cobra$feval)^alpha * sum(maxViolation)^beta , " ||", (C * cobra$feval)^alpha, " / ", sum(maxViolation)^beta ,"\n")
  #cat("<<< leaving subProb\n")
  return(y)
} # subProb()

# surrogate evaluation of 'f' for constraint optimization methods  - PHASE II
subProb2 <- function(x,cobra){

  if(any(is.nan(x))){
    warning("subProb2: x value is NaN, returning Inf")
    return(Inf)
  }
  
  if (cobra$trueFuncForSurrogates) {
    y<-cobra$fn(x)[1]
  } else {
    y<-predict.RBFinter(cobra$fitnessSurrogate,matrix(x,ncol=cobra$dimension))
  }
  
  return(y)
}

# surrogate evaluation of '\vec{g}' for constrained optimization methods - PHASE II
gCOBRA <- function(x,cobra) {

  if(any(is.nan(x))){
    warning("gCOBRA: x value is NaN, returning Inf")
    return(c(Inf,rep(Inf,cobra$nConstraints)))
  }
  
  if(cobra$CONSTRAINED)constraintPrediction <- calcConstrPred(x,cobra);
  h <- c()
  distance <- distLine(x,cobra$A)
  subC<-pmax((cobra$ro-distance),0)
  # h[1] <- sum(subC)*cobra$drFactor
  h[1] <- sum(subC)
  
  
  DBG=FALSE
  if (DBG & h[1]>0) {
    cat("gCOBRA: ",h,max(constraintPrediction),"\n")
    if (h<770) browser()
  }
  
  if(cobra$CONSTRAINED){
    h <- (-1.0)*c(h[1], constraintPrediction) # TODO -1* ... is required for COBYLA constraints, maybe also for other optimizers? 
  }else{
    h<--h
  }
  #cat("<<< leaving interpRBF\n")
  return(h)
}



### --- snip 

isresCobyla <- function(xStart,fn=subProb2,hin=gCOBRA, cobra) {
  #maxeval=cobra$seqFeval; 
  #subMin$feval=subMin$iter;
  subMin1 <- isres2(xStart,fn=fn,lower=cobra$lower,upper=cobra$upper,hin=hin, maxeval=cobra$seqFeval, cobra=cobra); 
  subMin2 <- nloptr::cobyla(xStart,fn=fn,lower=cobra$lower,upper=cobra$upper,hin=hin, control=list(maxeval=cobra$seqFeval,xtol_rel=cobra$seqTol), cobra=cobra); 
  if (subMin1$value < subMin2$value) {
    subMin1$nIsres <- cobra$nIsres+1;
    subMin1$nCobyla <- cobra$nCobyla;
    return (subMin1);
  } else {
    subMin2$nCobyla <- cobra$nCobyla+1;
    subMin1$nIsres <- cobra$nIsres;
    return (subMin2);
  }
}


# -----------------------------------------------------------------------------------------------
# ----------------  helper functions for subprob*, gCOBRA*  -------------------------------------
# -----------------------------------------------------------------------------------------------

# calculate constraintPrediction:
# -- this includes the proper handling of cobra$equHandle$active and 
# -- cobra$trueFunctionForSurrogates
# Known callers:
#   subprobPhaseI, subprob2PhaseI, gCOBRAPhaseI, 
#   subprob, gCOBRA (phase II)
#
calcConstrPred <- function(x,cobra) {
  if(cobra$equHandle$active){
    currentEps<-cobra$currentEps[length(cobra$currentEps)]
    if (cobra$trueFuncForSurrogates){   
      constraintPrediction1<- cobra$fn(x)[-1]
    } else {
      constraintPrediction1<- interpRBF(x,cobra$constraintSurrogates)
    }
    constraintPrediction1[cobra$equIndex]<- constraintPrediction1[cobra$equIndex]-currentEps  #this is supposed to create h(x)-mu
   # constraintPrediction2<- -constraintPrediction1[cobra$equIndex]-2*currentEps #this is supposed to be -h(x)-mu #why 2*currentEps?
    constraintPrediction2<- -constraintPrediction1[cobra$equIndex]-currentEps #this is supposed to be -h(x)-mu 
    
    constraintPrediction<- c(constraintPrediction1,constraintPrediction2)+cobra$EPS^2
    
  } else  {  # i.e. if(!cobra$equHandle$active)
    if (cobra$trueFuncForSurrogates) {
      constraintPrediction <-   cobra$fn(x)[-1]+cobra$EPS^2
    } else {
      constraintPrediction <-  interpRBF(x,cobra$constraintSurrogates)+cobra$EPS^2
    }
    
  }
  #/WK/ Bug fix: the above lines for cobra$trueFuncForSurrogates==TRUE were missing before
  
  return(constraintPrediction)
}

# constraint handling for suProb, the surrogate penalty function for unconstrained 
# optimization methods in cobraPhaseII.R
subProbConstraintHandling <- function(x,cobra,penalty1,penalty2,maxViolation,sigmaD,penaF)
{
  switch(cobra$constraintHandling,
         
         # J.A. Joines and C.R. Houck:  On the use of non-stationary penalty functions to solve 
         # nonlinear constrained optimization problems with GA's. In Proceedings of the First IEEE 
         # Conference on Evolutionary Computation, p. 579-584, (1994)
         JOINESHOUCK={
           C=0.5
           alpha=2
           beta=2
           y<-interpRBF(x, cobra$fitnessSurrogate) + (C * cobra$feval)^alpha * sum(maxViolation)^beta 
         },
         
         # A.E. Smith and D.M. Tate: Genetic optimization using a penalty function. In Proceedings of the 
         # Fifth International Conference on Genetic Algorithms p. 499-505, (1993)
         SMITHTATE={
           #TODO
           lambda = 1
           beta1 = 3
           beta2 = 2
           kappa = 1
           d = sum(maxViolation)^kappa
           y <- interpRBF(x, cobra$fitnessSurrogate) + penalty
         },
         
         # D.W. Coit, A.E. Smith and D.M. Tate: Adaptive penalty methods for genetic optimization of 
         # constrained combinatorial problems. In INFORMS Journal on Computing 8(2), (1996)
         COIT={
           fFeas <- cobra$fbest
           fAll <- min(cobra$Fres)
           NFT <- rep(0.05, cobra$nConstraints)
           kappa <- 1
           y <- interpRBF(x, cobra$fitnessSurrogate) + (fFeas - fAll) * sum((maxViolation/NFT)^kappa)
         },
         
         # T. Baeck and S. Khuri: An evolutionary heuristic for the maximum independent set problem 
         # In Proceedings of the First IEEE Conference on Evolutionary Computation, p. 531-535, (1994)
         BAECKKHURI={
           #TODO
           K=10^9
           p = cobra$nConstraints
           s <- length(cobra$numViol)
           penalty <- K - s*(K/cobra$nConstraints)
           y <- interpRBF(x, cobra$fitnessSurrogate) + penalty
         },
         
         DEFAULT={
           y<-interpRBF(x, cobra$fitnessSurrogate) + (penalty1*sigmaD[1] + penalty2)*penaF[1]
         }
  )
  
  #PKDEBUG=TRUE|FALSE
  #if(PKDEBUG){# Joines and Houck 1994
  #  C=5
  #  alpha=2
  #  beta=2
  #  y<-interpRBF(x, cobra$fitnessSurrogate) + (C * cobra$feval)^alpha * sum(maxViolation)^beta 
  #}else{y<-interpRBF(x, cobra$fitnessSurrogate) + (penalty1*sigmaD[1] + penalty2)*penaF[1]}
  
  return(y)
}
