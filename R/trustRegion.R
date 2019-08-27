#trustRegion.R
#'
#' Performs trust region refinement
#'
#' If \code{cobra$TrustRegion==TRUE}  (see \code{\link{cobraInit}}), then the \code{trustRegion} functionality is applied
#' every iteration in order to refine the best solution so far.
#' This function builds a local model around the best solution and runs a local search in the trust region 
#' to refine the best solution and find a better solution in the neighborhood.
#'
#' @param cobra     an object of class \code{cobra}, which is basically a list  (see \code{\link{cobraInit}})
#' @param center    [cobra$xbest] the center of the trust region
#' @return the modified \code{cobra} with new/updated elements
#'      \item{TRDONE}{ logical, is \code{TRUE} if there are more than d+1 points in the trusted 
#'            region and thus surrogates can be trained. Otherwise \code{FALSE}.}
#'      \item{trustregX}{ if \code{TRDONE==TRUE} the refined solution from the trust-region call,
#'            otherwise \code{NA} }
#'    If \code{TRDONE==TRUE} the relevant lists and counters  (\code{A,Fres,df,...}) 
#'    of \code{cobra} will be updated in \code{\link{cobraPhaseII}} as well.
#'    
#' @author Samineh Bagheri (\email{samineh.bagheri@@th-koeln.de})

trustRegion<-function(cobra, center=cobra$xbest){
  delta<-cobra$radi[length(cobra$radi)]
  TRDONE<-FALSE
  trustregX<-rep(NA,length(cobra$xbest))
  
  while(delta <= cobra$TRlist$radiMax && TRDONE==FALSE){
    TRlower<-pmax(center-delta,cobra$lower)
    TRupper<-pmin(center+delta,cobra$upper)
    ind<-c()
    
    switch(cobra$TRlist$shape,
           "cube" ={ind<-which(sapply(1:nrow(cobra$A),indCheck<-function(i){
             all((TRupper-cobra$A[i,] >= 0 ) & (cobra$A[i,]-TRlower >= 0))
           })) },
           "sphere"={radius<-sqrt(ncol(cobra$A))*delta;
                     distances<-as.matrix(stats::dist(rbind(center,cobra$A)))[-1,1];
                     ind<-which(distances<radius);
                     ind<-ind-1},
           "hypercross"={ ind<-which(sapply(1:nrow(cobra$A),indCheck<-function(i){
             (TRupper-cobra$A[i,] >= 0 ) && (cobra$A[i,]-TRlower >= 0)
           })) 
           #(only for hypercross case)
           TRlower<-apply(cobra$A[ind,],2,min)
           TRupper<-apply(cobra$A[ind,],2,max)
           
           },
           stop("the selected shape for trust region is not valid")
    ) # switch
    dimension<-ncol(cobra$A)
    if(length(ind) < dimension+1 ){ # for building an RBF we need at least d+1 points, if we have less points then building 
      trustregX<-rep(NA,dimension)            # a local model in the defined trust region is not possible
      #print("Trust Region cannot be performed ")
      cobra$TRDONE<-FALSE
      if(cobra$DEBUG_TR)cat("[",nrow(cobra$A),"]TR DEBUG ::TRpop=",length(ind),"\n")
      delta<-2*delta
    }else{ 
      cobra$TRDONE<-TRUE
      TRDONE<-TRUE
      cobra$RSDONE<-c(cobra$RSDONE,"TR")
      #build the local model over the trust region
    }
  } # while()  
  
  if(cobra$TRDONE){
    ############################Training the local model###############################################
    cat(paste(">> Training", cobra$TRlist$shape,"Trust Region [delta=",delta,", n=",length(ind), "]" ,cobra$RBFmodel,"surrogates","...\n" ))
    TRA<-cobra$A[ind,]
    TRGres<-cobra$Gres[ind,]
    TRFres<-cobra$Fres[ind]
    if(cobra$DEBUG_TR){
    TRFR<-(max(TRFres)-min(TRFres)) #FR feature for the Trust region area
    cobra$TFR<-c(cobra$TFR,TRFR)
    cat("[",nrow(cobra$A),"]TR DEBUG ::TRFR=",TRFR,"\n")}
   
    sw=switch(cobra$RBFmodel,
              "cubic" =     {constraintSurrogatesTR <- trainCubicRBF(TRA,TRGres,ptail=cobra$ptail,squares=cobra$squares)
                             fitnessSurrogateTR <- trainCubicRBF(TRA,TRFres,ptail=cobra$ptail,squares=cobra$squares)},
              "Gaussian"=   {constraintSurrogatesTR <- trainGaussRBF(TRA,TRGres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,widthFactor=cobra$widthFactor);
                             fitnessSurrogateTR <- trainGaussRBF(TRA,TRFres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,widthFactor=cobra$widthFactor)},
              "MQ"=   {constraintSurrogatesTR <- trainMQRBF(TRA,TRGres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,widthFactor=cobra$widthFactor);
              fitnessSurrogateTR <- trainMQRBF(TRA,TRFres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,widthFactor=cobra$widthFactor)}
    )
    
    #---generate random points in the trust region area
    if(cobra$DEBUG_TR){
     popSize<-1000
     RandomPopTR<-lhs::randomLHS(popSize,ncol(cobra$A))
     TRS<-apply(RandomPopTR,1,function(x){
       return(interpRBF(x, fitnessSurrogateTR))
     })
     TRSR<-max(TRS)-min(TRS)   #range of the surrogate fitness in the area of trust region, 
                               #measured by popSize points in the trust region
     cobra$TSR<-c(cobra$TSR,TRSR)
     if(cobra$DEBUG_TR)cat("[",nrow(cobra$A),"]TR DEBUG ::SR/FR=",TRSR/TRFR,"\n")
     cobra$Fratio<-c(cobra$Fratio,TRSR/TRFR)
     cobra$TRpop<-c(cobra$TRpop,length(ind))
     cobra$TRdelta<-c(cobra$TRdelta,delta)
    }
    #trainTRsurrogates
    ############################Training the local model finished###############################################
    cobra$fitnessSurrogateTR<-fitnessSurrogateTR
    cobra$constraintSurrogatesTR<-constraintSurrogatesTR
    
    if (cobra$DEBUG_TRU) {
      if (cobra$dimension!=2) stop("cobra$DEBUG_TRU currently only for d=2")
      if (cobra$rescale==F) stop("cobra$DEBUG_TRU currently only for cobra$rescale==TRUE")
      
      ##this will be deleted if it is not going to be used
      rescaler<-function(x,TOL,TOU,FROML,FROMU){
        x<-sapply(1:length(x) , function(i){scales::rescale(x[i],to=c(TOL[i],TOU[i]),from=c(FROML[i],FROMU[i]))
        })

        return(x)
      }
      TRlowerO<-rescaler(TRlower,TOL=cobra$originalL,TOU=cobra$originalU,FROML=cobra$lower,FROMU=cobra$upper)
      TRupperO<-rescaler(TRupper,TOL=cobra$originalL,TOU=cobra$originalU,FROML=cobra$lower,FROMU=cobra$upper)
      
      N=100

        
      X=outer(rep(1,N),seq(TRlower[1],TRupper[1],length=N))  
      Y=t(outer(rep(1,N),seq(TRlower[2],TRupper[2],length=N)))   
      #Xr=outer(rep(1,N),seq(cobra$newlower,cobra$newupper,length=N))  
      #Yr=t(Xr)   
      newWin <- is.null(cobra$TrWinExists)
      localModel <- drawSurrogate3d(fitnessSurrogateTR,cobra$fName,X,Y,X,Y,newWindow=newWin,
                            lower=TRlower,upper=TRupper,device=cobra$localDevice,title="local",cobra$A)
      ZS<-localModel$ZS
      cobra$localDevice<-localModel$devNum
      cobra$TrWinExists <- TRUE

      TRFres2 <- TRFres;
      rgl::points3d(cbind(TRA,TRFres2),color="white",size=10)
      #       rngX = TRupper[1] - TRlower[1]
      #       rngY = TRupper[2] - TRlower[2]
      #       rngZ = max(max(TRFres2)-min(TRFres2),max(ZS)-min(ZS))
      #       rgl::aspect3d(1/rngX,1/rngY,1/rngZ)
      
      approxErr <- NA
      VISUALIZE_TRUE_Z <- TRUE
      if (VISUALIZE_TRUE_Z) {
        # this branch will cost some time
        fitFunc <- function(x) {
          z=cobra$fn(x)
          return(z[1])
        }
        Z = X*0
        for (i in 1:N)
          for (j in 1:N)
            Z[i,j] = fitFunc(c(X[i,j],Y[i,j]))      
        drawSurface3d(X,Y,Z)
        approxErr = sqrt(sum((Z-ZS)^2))/(N*N)
      }

      newpnt <- TRA[nrow(TRA),]
      cat(sprintf("N=%d, new pnt  at (%7.2f,%7.2f), surrogate range = %7.2f, approx err=%7.2f\n",
                  fitnessSurrogateTR$npts,newpnt[1],newpnt[2],max(ZS)-min(ZS),approxErr))
      wm <- which.min(ZS)
      minpnt <- c(X[wm],Y[wm])
      cat(sprintf("N=%d, surr min at (%7.2f,%7.2f)\n",
                  fitnessSurrogateTR$npts,minpnt[1],minpnt[2]))
      
      DO_SNAPSHOTS=F
      if (DO_SNAPSHOTS) {
        if (fitnessSurrogateTR$npts %% 2 == 0) {
          if (fitnessSurrogateTR$npts==6) browser()
          if (!file.exists("images.d")) dir.create("images.d")
          rgl::rgl.snapshot(sprintf("images.d/%s-TR-%03d.png",cobra$fName,
                                    cobra$fitnessSurrogate$npts))       
        }
      }
    } # if (cobra$DEBUG_TRU)
    
    cobra$TRcounter <- cobra$TRcounter + 1
    cobra$TRgama<-cobra$XI[((cobra$TRcounter) %% length(cobra$XI))+1] 
    #cobra$TRgama<-0
    ro<-cobra$TRgama*cobra$l  
    
    subProbTR2 <- function(x){
      y<-predict.RBFinter(fitnessSurrogateTR,matrix(x,ncol=dimension))
      return(y)
    }
    
    gCOBRATR <- function(x) {
      h <- c()
      distance <- distLine(x,cobra$A[ind,])
      distCent <- sqrt(sum((x-center)^2))
      subC<-pmax((ro-distance),0)
     # h[1] <- sum(subC)*cobra$drFactor
      h[1] <- sum(subC)
      
      h[2] <- max(delta-center,0)       # /WK/ new constraint, experimental: sphere around center
      #cat(">>> Entering interpRBF\n")
      constraintPredictionTR <- interpRBF(x,constraintSurrogatesTR)+cobra$epsilonInit^2
      #h <- (-1.0)*c(h[1],h[2], constraintPredictionTR) # TODO -1* ... is required for COBYLA constraints, maybe also for other optimizers? 
      h <- (-1.0)*c(h[1], constraintPredictionTR) 
      #cat("<<< leaving interpRBF\n")
      return(h)
    }
    
    #generate a random starting point in the trust region
    #xStart<-runif(ncol(cobra$A),min=TRlower,max=TRupper)
   # assert("center of TR is not on the current best point",all(cobra$xbest==center))
    xStart<-center
    switch(cobra$seqOptimizer,
           COBYLA={ subMin<-nloptr::cobyla(xStart,fn=subProbTR2,lower=TRlower,upper=TRupper,hin=gCOBRATR,control=list(maxeval=cobra$seqFeval,xtol_rel=cobra$seqTol)); subMin$feval=subMin$iter },
           ISRES ={ subMin<-isres2(xStart,fn=subProbTR2,lower=TRlower,upper=TRupper,hin=gCOBRATR, maxeval=cobra$seqFeval); subMin$feval=subMin$iter}
           #RANDOMSEARCH = { subMin <- randomSearch(xStart, fn=subProb2, lower=TRlower,upper=TRupper,control=list(maxfeval=cobra$seqFeval, sd=0.05))}
    )

    trustregX<-subMin$par
    if(cobra$DEBUG_TR)cat("[",nrow(cobra$A),"]TR DEBUG ::Move length=",as.numeric(stats::dist(rbind(cobra$xbest,trustregX))),"\n")

  }
  cobra$trustregX<-trustregX
if(is.na(trustregX[1])){
  output<-NA
}else{
  output<-cobra$fn(trustregX)[1]
}

if(cobra$DEBUG_TR)print(sprintf("%s.[%d]: %f %f | %f" 
        , "Refined x " ,nrow(cobra$A), trustregX[1] ,trustregX[2] , output ))
  return(cobra)
}

