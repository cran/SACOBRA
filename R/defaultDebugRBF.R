#defaultDebugRBF
#'
#' Default settings for debug visualization RBF (only for \code{d==2})
#' 
#' Sets default values for debug visualization RBF  of SACOBRA.  
#'
#' @return  \code{DEBUG_RBF} a list of the follwing elements:
#'    \item{active}{If set to TRUE then \code{\link{debugVisualizeRBF}} is called every \code{DEBUG_RBF$every} iterations}
#'    \item{overlayTrueZ}{If set to TRUE overlay the true objective function}
#'    \item{DO_SNAPSHOT}{do rgl.snapshot every \code{DEBUG_RBF$every} iteration and store it \cr in \code{sprintf("images.d/\%s-\%03d.png",cobra$fname,npts)}}
#'    \item{every}{Frequency of calling the \code{\link{debugVisualizeRBF}} function}
#' @seealso   \code{\link{debugVisualizeRBF}}
#' @export
defaultDebugRBF<-function(){
  DEBUG_RBF=list(active=FALSE,
                 overlayTrueZ=FALSE,
                 DO_SNAPSHOT=FALSE,
                 every=2)
  return(DEBUG_RBF)  
}


# debugVisualizeRBF
#'
#' Optional visualization of surrogate models for 2-dimensional problems, d=2.
#' 
#' @param \code{cobra}  parameter list, we need here 
#'     \describe{
#'       \item{\code{DEBUG_RBF$overlayTrueZ}}{=T: overlay the true objective function}
#'       \item{\code{DEBUG_RBF$DO_SNAPSHOT}}{=T: do rgl.snapshot every \code{DEBUG_RBF$every} iteration and store it in \code{sprintf("images.d/\%s-\%03d.png",cobra$fname,npts)} }
#'       \item{\code{DEBUG_RBF$every}}{  see above}}
#' @param fitnessSurrogate  The surrogate model for the objective function
#' @param A                 The whole population of solutions which is a matrix of n times d, where n is the number of solutions and d is the dimension of problem, here d=2
#' @param Fres              The objective value on the evaluated solutions \code{A}
#' @return \code{cobra}, an object of class COBRA, with modified items:
#'    \item{\code{TrueZ}}{True objective values of the evaluated points used for viusalization}
#'    \item{\code{TrueFeas}}{A vector of values 0: for infeasible evaluated points or 1: for feasible points.}
#'    \item{\code{globalDevice}}{}
# Detail:
#' This function is called only if cobra$DEBUG_RBF$active==T.
#' An assertion fires if cobra$dimension!=2 or if cobra$rescale==F.
#' @keywords internal   
debugVisualizeRBF <- function(cobra,fitnessSurrogate,A,Fres){
  if (cobra$dimension!=2) stop("cobra$DEBUG_RBF currently only for d=2")
  if (cobra$rescale==F) stop("cobra$DEBUG_RBF currently only for cobra$rescale==TRUE")
  
  N=100
  X=outer(rep(1,N),seq(cobra$originalL[1],cobra$originalU[1],length=N))  
  Y=t(outer(rep(1,N),seq(cobra$originalL[2],cobra$originalU[2],length=N)))   
  Xr=outer(rep(1,N),seq(cobra$newlower,cobra$newupper,length=N))  
  Yr=t(Xr)   
  
  if (is.null(cobra$TrueZ)) {
    # we first set Z <- NA for infeasible points, take the max w/o NA and 
    # then set all NA's in Z to this max value
    fitFunc <- function(x) {
      z <- cobra$originalfn(x)
      if (any(z[-1]>0)) z[1] <- NA      # tag infeasible points with NA
      return(z[1])
    }
    Z = X*0
    for (i in 1:N)
      for (j in 1:N)
        Z[i,j] = fitFunc(c(X[i,j],Y[i,j]))  
    cobra$TrueFeas = X*0
    cobra$TrueFeas[!is.na(Z)] <- 1
    infeasVal <- max(Z,na.rm=TRUE)
    Z[is.na(Z)] <- infeasVal
    cat("infeasVal: ", infeasVal,"\n")
    cobra$TrueZ <- Z
  }
  
  newWin <- (cobra$initDesPoints==nrow(A))
  overlayTrueZ <- cobra$DEBUG_RBF$overlayTrueZ
  if(cobra$trueFuncForSurrogates)fitnessSurrogate<-cobra$fn
  globalModel <- drawSurrogate3d(fitnessSurrogate,cobra$fName,X,Y,Xr,Yr,newWindow=newWin, 
                                 TrueZ=cobra$TrueZ,overlayTrueZ=overlayTrueZ,
                                 lower=cobra$originalL,upper=cobra$originalU,
                                 device=cobra$globalDevice,title="global",A=cobra$A)
  ZS<-globalModel$ZS
  cobra$globalDevice<-globalModel$devNum
  origBest<-t(sapply(1:length(cobra$xbest),function(i)inverseRescale(cobra$xbest[i],cobra)))
  
  if(class(fitnessSurrogate)!="function"){
  origA <- t(sapply(1:nrow(A),function(i)inverseRescale(A[i,],cobra)))
  CENTROIDS<-t(sapply(1:nrow(fitnessSurrogate$xp),function(i)inverseRescale(fitnessSurrogate$xp[i,],cobra)))
  Fres2 <- Fres;
  rgl::points3d(cbind(origA,Fres2),color="white",size=10)
  }
  # rgl::points3d(cbind(origBest,cobra$fbest),color="yellow",size=15)
   solu<-t(sapply(1:length(cobra$solu),function(i)inverseRescale(cobra$solu[i],cobra)))
   rgl::points3d(cbind(solu,1e+6),color="red",size=15)
   
  approxErr = sqrt(sum((cobra$TrueZ-ZS)^2))/(N*N)
  newpnt <- inverseRescale(A[nrow(A),],cobra)
  cat(sprintf("N=%d, new pnt  at (%7.2f,%7.2f), surrogate range = %7.2f, approx err=%7.2f\n",
              nrow(A),newpnt[1],newpnt[2],max(ZS)-min(ZS),approxErr))
  ZS[cobra$TrueFeas==0] <- max(ZS)
  wm <- which.min(ZS)
  minpnt <- c(X[wm],Y[wm])
  cat(sprintf("N=%d, surr min at (%7.2f,%7.2f)\n",
              nrow(A),minpnt[1],minpnt[2]))
  # only diagnostics, needed for cobra$df & cobra$df2 /WK/
  testit::assert("cobra$solu is missing", !is.null(cobra$solu))
  solu <- cobra$solu; 
  if (cobra$rescale) 
    if (is.matrix(solu)) {
      solu <- t(sapply(1:nrow(solu),function(i){ forwardRescale(solu[i,],cobra)}))
    } else {
      solu <- forwardRescale(solu,cobra);
    }
  # now solu is always in *rescaled* input space
  dL <- distLine(solu,rbind(newpnt,minpnt))
  print(c(dL,XI=cobra$df$XI[nrow(A)]))
  #DO_SNAPSHOT=T; every=2;
  if (cobra$DEBUG_RBF$DO_SNAPSHOT) {
    if (nrow(cobra$A)==cobra$initDesPoints)     #|| nrow(cobra$A)==40
      browser();   # possibility to adjust RGL-plot before first snapshot
    snapShot3d(fitnessSurrogate,cobra$DEBUG_RBF$every,sprintf("%s-03b",cobra$fName),A);
  }
  return(cobra)
} # debugVisualizeRBF