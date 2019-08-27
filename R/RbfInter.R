
#Wolfgang Konen
#Cologne University of Applied Science
#RbfInter.R
#
# RBF interpolation model
#
 
#' Euclidean distance of \code{x} to all \code{xp}
#'
#' Euclidean distance of \code{x} to a line of points \code{xp}
#'
#' \code{distLine} is up to 40x faster than using \code{\link[stats]{dist}} and taking only the first row or 
#' column of the distance matrix returned. 
#'
#' @param x       vector of dimension d 
#' @param xp      n points \eqn{x_i} of dimension d are arranged in (n x d) matrix \code{xp}.
#'                If \code{xp} is a vector, it is interpreted as (n x 1) matrix, i.e. d=1. 
#'
#' @return        vector of length n, the Euclidean distances
#' 
distLine <- function(x,xp) {
  if (is.vector(xp)) xp <- as.matrix(xp)
  z = outer(rep(1,nrow(xp)),x) - xp
  z = sqrt(rowSums(z*z))
}

# helper for train...RBF function
calcRHS <- function(U,d2) {
  if (length(d2)!=0) {
    if (is.vector(U))
      rhs = c(U,rep(0,d2))
    else if (is.matrix(U)) 
      rhs = rbind(U,matrix(0,d2,ncol(U)))
    else 
      stop("U is neither vector nor matrix!")
  } else {
    rhs = U
    if (!is.vector(U) && !is.matrix(U))
      stop("U is neither vector nor matrix!")
  }
  #browser()
  return (rhs);
}  



#Stable calculation of  matrix M by using singular value decomposition svd
svdInv <- function(M) {
  eps = 1e-14
  s = svd(M)
  invD <- 1/s$d
  invD[abs(s$d/s$d[1])<eps] <- 0
  invM <- s$v %*% diag(invD) %*% t(s$u)
  
  return(invM)
}


#Internal function train RBF function
#After Phi matrix being built in trainCubic, trainGaussian and trainMQ functions, 
#it will be passed to this function to become augmented and find the RBF model parameters 
#by inverting the augmented matrix  
# 
trainRBF<-function(phi,U,ptail,squares,xp,type,DEBUG2,width=NA,rho=0){

  npts=nrow(xp)
  #AUGMENTED<-T
  d2<-NULL   #size of the polynomial tail
  pMat<-NULL #polynomial matrix
  regMat<-NULL
  
  if(!ptail && !squares){
    M = phi
  }else{
    if (ptail)  pMat=cbind(e=rep(1.0,npts),xp)          # linear tail LH(1,x1,x2,...)
    if (squares) pMat=cbind(pMat,xp*xp)     # ... plus direct squares x1^2, x2^2, ...
    
    d2<-ncol(pMat)
    nMat=matrix(0,d2,d2)                # matrix of zeros
    M<-cbind(phi,pMat)
    
      QQ=t(pMat)
      M<-rbind(M,cbind(QQ,nMat)) 
      M<-M+diag(nrow(M))*rho^2
  }
  
  invM = svdInv(M)
  rhs <- calcRHS(U,d2) 
  coef = invM %*% rhs
  ### 
  ### The linear-equation solution with 'solve' is more accurate when matrix M has full rank.
  ### But it crashes for rank-deficient or near-singular matrices M.
  ### This version only for vector-like input U:
  ###
  #coef2 = solve(M,rhs)
  #print(c(max(abs(M %*% coef - rhs)),max(abs(M %*% coef2 - rhs)),max(abs(rhs))))

  # this check does not yet work as expected, leave DEBUG==FALSE /WK/
  DEBUG=FALSE
  if (DEBUG) {
    inv2 = MASS::ginv(M)
    coef2 = inv2 %*% rhs
    testit::assert("rows and cols invM and inv2 do not match",nrow(invM)==nrow(inv2),ncol(invM)==ncol(inv2))
    if (nrow(invM)<=ncol(invM)) 
      testit::assert("invM-check 1 failed",max(abs(invM %*% M - diag(rep(1,nrow(invM)))))<1e-5)
    else
      testit::assert("invM-check 2 failed",max(abs(M %*% invM - diag(rep(1,ncol(invM)))))<1e-5)    
  }
 
  rbf.model = list(coef=coef
                   ,xp=xp
                   ,d=d2
                   ,npts=nrow(xp)
                   ,ptail=ptail
                   ,squares=squares
                   ,type=type,width=width)  
  if (DEBUG2) {
    rbf.model$M = M
    rbf.model$rhs = rhs
  }
  class(rbf.model) <- c("RBFinter","list")
  return(rbf.model)
}






#----------------------------------------------------------------------------------
#' Fit cubic RBF interpolation to training data X for d>1.
#'
#' The model at a point \eqn{z=(z_1,...,z_d)} is fitted using n sample points \eqn{x_1, ..., x_n} 
#' \cr
#'    \deqn{ s(z) = \lambda_1*\Phi(||z-x_1||)+... +\lambda_n*\Phi(||z-x_n||)
#'                  + c_0 + c_1*z_1 + ... + c_d*z_d  }
#' \cr
#' where \eqn{\Phi(r)=r^3} denotes the cubic radial basis function. The coefficients \eqn{\lambda_1, 
#' ..., \lambda_n, c_0, c_1, ..., c_d} are determined by this training procedure.\cr
#' This is for the default case \code{squares==FALSE}. In case \code{squares==TRUE} 
#' there are d additional pure square terms and the model is
#' \cr
#'    \deqn{ s_{sq}(z) = s(z) + c_{d+1}*z_1^2 + ... + c_{d+d}*z_d^2 } 
#' \cr
#' In case \code{ptail==FALSE} the polynomial tail (all coefficients \eqn{c_i}) is omitted completely.
#'
#' The linear equation system is solved via SVD inversion. Near-zero elements 
#' in the diagonal matrix \eqn{D} are set to zero in \eqn{D^{-1}}. This is numerically stable 
#' for rank-deficient systems.
#'
#' @param xp      n points \eqn{x_i} of dimension d are arranged in (n x d) matrix \code{xp}
#' @param U       vector of length n, containing samples \eqn{f(x_i)} of 
#'                the scalar function \eqn{f} to be fitted \cr
#'                - or - \cr
#'                (n x m) matrix, where each column 1,...,m contains one vector of samples
#'                \eqn{f_j(x_i)} for the m'th model, j=1,...,m
#' @param squares [FALSE] flag, see description
#' @param ptail   [TRUE] flag, see description
#' @param rho     [0.0] experimental: 0: interpolating, >0, approximating (spline-like) 
#'                Gaussian RBFs
#' @param DEBUG2  [FALSE] if TRUE, save \code{M} and \code{rhs} on return value
#' @param width   [NA] non relevant for the parameter-free cubic RBF              
#'                
#' @return \code{rbf.model},  an object of class \code{RBFinter}, which is basically a list 
#' with elements:
#'      \item{coef}{  (n+d+1 x m) matrix holding in column m the coefficients for the m'th 
#'                    model:      \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d}.  
#'                    In case \code{squares==TRUE} it is an (n+2d+1 x m) matrix holding  
#'                    additionally the coefficients \eqn{c_{d+1}, ..., c_{d+d}}.}                    
#'      \item{xp}{  matrix xp   }
#'      \item{d}{  size of the polynomial tail. If \code{length(d)==0} it means no polynomial tail will be used for the model. In case of ptail==T && squares==F d will be dimension+1 and in case of ptail==T && squares==T d will be 2*dimension+1 }
#'      \item{npts}{  number n of points \eqn{x_i} }
#'      \item{ptail}{TRUE or FALSE (see description)}
#'      \item{squares}{  TRUE or FALSE (see description)  }
#'      \item{type}{  "CUBIC"}
#'      \item{width}{NA, irrelevant for the parameter-free cubic RBF}
#'      
#' @seealso   \code{\link{trainGaussRBF}}, \code{\link{trainMQRBF}} \code{\link{predict.RBFinter}}, \code{\link{interpRBF}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@th-koeln.de}), Samineh Bagheri (\email{samineh.bagheri@@th-koeln.de})
trainCubicRBF <- function(xp, U, ptail=TRUE, squares=FALSE, rho=0.0,DEBUG2=FALSE,width=NA) {

  edist=as.matrix(stats::dist(xp,upper=T))         # euclidean distance matrix
  phi=edist*edist*edist                     # cubic RBF (npts x npts matrix)
                                            
  rbf.model<-trainRBF(phi,U,ptail,squares,xp,type="CUBIC",DEBUG2,width=NA,rho=rho)
  
  return(rbf.model)
}

#----------------------------------------------------------------------------------
#' Fit Gaussian RBF model to training data for d>1.
#' 
#' The model for a point \eqn{z=(z_1,...,z_d)} is fitted using n sample points \eqn{x_1, ..., x_n} 
#' \cr
#'    \deqn{ s(z) = \lambda_1*\Phi(||z-x_1||)+... +\lambda_n*\Phi(||z-x_n||)
#'                  + c_0 + c_1*z_1 + ... + c_d*z_d  }
#' \cr    
#' where \eqn{\Phi(r)=\exp(-r^2/(2*\sigma^2))} denotes the Gaussian radial basis function with width
#' \eqn{\sigma}. The coefficients \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d} are determined 
#' by this training procedure.\cr
#' This is for the default case \code{squares==FALSE}. In case \code{squares==TRUE} 
#' there are d additional pure square terms and the model is
#' \cr
#'    \deqn{ s_{sq}(z) = s(z) + c_{d+1}*z_1^2 + ... + c_{d+d}*z_d^2 } 
#' In case \code{ptail==FALSE} the polynomial tail (all coefficients \eqn{c_i}) is omitted completely.
#'  
#' The linear equation system is solved via SVD inversion. Near-zero elements 
#' in the diagonal matrix \eqn{D} are set to zero in \eqn{D^{-1}}. This makes  
#' rank-deficient systems numerically stable.
#' 
#' @param xp      n points \eqn{x_i} of dimension d are arranged in (n x d) matrix \code{xp}
#' @param U       vector of length n, containing samples \eqn{u(x_i)} of 
#'                the scalar function \eqn{u} to be fitted \cr
#'                - or - \cr
#'                (n x m) matrix, where each column 1,...,m contains one vector of samples
#'                \eqn{u_j(x_i)} for the m'th model, j=1,...,m
#' @param squares [FALSE] flag, see 'Description'
#' @param ptail   [TRUE] flag, see description
#' @param width   [-1] either a positive real value which is the constant width \eqn{\sigma} for all 
#'                Gaussians in all iterations, or -1. If -1, the appropriate width \eqn{\sigma} is 
#'                calculated anew in each iteration with one of the rules \code{RULE},
#'                based on the distribution of data points \code{xp}.              
#' @param RULE    ["One"] one out of ["One" | "Two" | "Three"], different rules for automatic 
#'                estimation of width \eqn{\sigma}. Only relevant if \code{width = -1},   
#' @param widthFactor [1.0] additional constant factor applied to each width \eqn{\sigma} 
#' @param rho     [0.0] experimental: 0.0: interpolating, >0.0, approximating (spline-like) 
#'                Gaussian RBFs
#' @param DEBUG2  [FALSE] if TRUE, save \code{M} and \code{rhs} on return value
#' @return \code{rbf.model},  an object of class \code{RBFinter}, which is basically a list 
#'      with elements:
#'      \item{coef}{  (n+d+1 x m) matrix holding in column m the coefficients for the m'th 
#'                    model:      \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d}.  
#'                    In case \code{squares==TRUE} it is an (n+2d+1 x m) matrix holding  
#'                    additionally the coefficients \eqn{c_{d+1}, ..., c_{d+d}}.}
#'      \item{xp}{ matrix xp   }
#'      \item{d}{ size of the polynomial tail. If \code{length(d)==0} it means no polynomial tail will be used for the model. In case of ptail==T && squares==F d will be dimension+1 and in case of ptail==T && squares==T d will be 2*dimension+1 }
#'      \item{npts}{  number n of points \eqn{x_i} }
#'      \item{ptail}{TRUE or FALSE (see description)}
#'      \item{squares}{  TRUE or FALSE (see description)}
#'      \item{width}{ the calculated width \eqn{\sigma} }
#'      \item{type}{  "GAUSS"}
#'
#' @seealso   \code{\link{trainCubicRBF}}, \code{\link{predict.RBFinter}}, \code{\link{interpRBF}}
#' @author Wolfgang Konen, Samineh Bagheri
#'          
#----------------------------------------------------------------------------------
trainGaussRBF <- function(xp, U, ptail=TRUE, squares=FALSE, 
                          width,RULE="One",widthFactor=1.0, rho=0.0, DEBUG2=F) {
  #browser()
  if(width==-1){#When no width is given by user then width is selected automatically
    #RULE<-"One"
    #Automatic width adaptation is in testing phase therefore we would like to test several rules
    #from literature and then select teh most relevant ones
    
    #>> rule number One:
    #The rule is taken from Benoudjit,2002 (But originally from Haykin,1999)
    #It works based on  finding a compromise between locality and smoothness
    
    #>> rule number Two:
    #The rule is taken from Benoudjit,2002 (But originally from Moody and Darken,1989)
    #a vector is returned as width and width[i] is the width factor for the ith RBF
    
    #>> rule number Three:
    # This rule is taken from Sun and Jin,2013
    # The idea is basically adapting the width to the smallest interval of two points in each coordiante

    switch(RULE,
           "One" =     width<-(max(stats::dist(xp)))/sqrt(2*nrow(xp)),
           "Two" = {
              k<-2
              width<-sapply(1:nrow(xp),function(i){
                centroid<-matrix(xp[i,],nrow=1)
                ncentroid<-xp[-i,]
                #browser()
                #index<-FNN::knnx.index(ncentroid,centroid, k=k) #finding the index of k nearest neigbors
                mindists<-FNN::knnx.dist(ncentroid,centroid, k=k)
                #mindists<-mindists^2
                y<-(1/k)*(sqrt(sum(mindists*mindists)))
                return(y)
              })},
           "Three" = {width<-c()
                      for(i in 1:ncol(xp)){
                        interval<-max(xp[,i])-min(xp[,i])
                      width<-min(width,interval)  
             
           }}
           )
    verboseprint(verbose=0,important=FALSE,paste("Automatic adjustment of RBF width=",width[1]))
    
  }
 

  edist=as.matrix(stats::dist(xp,upper=T))
  width = width*widthFactor
  if (length(width)==1) {                   # /WK/ bug fix: former code was not correct for the case 
    wmat=width^2                            # that width may be a vector of length nrow(xp)
  } else {
    w = outer(rep(1,length(width)),width)
    wmat = w*w                              # this choice leads to non-symmetric phi, but it 
                                            # has smallest error on RbfFit2.R (multiply each width by 
                                            # factor 3 - 5)
    #wmat = w*t(w)                          # this is (wmat)_ij = width_i * width_j, it
                                            # guarantees that phi will be a symmetric matrix,
                                            # but it has bigger error  
  }
  phi <- exp(-0.5*(edist*edist)/wmat)       # Gaussian RBF with width=width
                                            # /WK/ experimental: rho>0 means spline approximation
  # phi <- phi - diag(npts)*npts*rho          # /WK/ instead of exact interpolating RBF (rho=0)

  rbf.model<-trainRBF(phi,U,ptail,squares,xp,type="GAUSS",DEBUG2,width=width,rho=rho)
  
  return(rbf.model)
}

#----------------------------------------------------------------------------------
#' Fit multiquadric RBF model to training data for d>1.
#'
#' The model for a point \eqn{z=(z_1,...,z_d)} is fitted using n sample points \eqn{x_1, ..., x_n} 
#' \cr
#'    \deqn{ s(z) = \lambda_1*\Phi(||z-x_1||)+... +\lambda_n*\Phi(||z-x_n||)
#'                  + c_0 + c_1*z_1 + ... + c_d*z_d  }
#' \cr    
#' where \eqn{\Phi(r)=\sqrt{(1+(r/\sigma)^2)}} denotes the multiquadrics radial basis function with width
#' \eqn{\sigma}. The coefficients \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d} are determined 
#' by this training procedure.\cr
#' This is for the default case \code{squares==FALSE}. In case \code{squares==TRUE} 
#' there are d additional pure square terms and the model is
#' \cr
#'    \deqn{ s_{sq}(z) = s(z) + c_{d+1}*z_1^2 + ... + c_{d+d}*z_d^2 } 
#' In case \code{ptail==FALSE} the polynomial tail (all coefficients \eqn{c_i}) is omitted completely.
#'
#' The linear equation system is solved via SVD inversion. Near-zero elements 
#' in the diagonal matrix \eqn{D} are set to zero in \eqn{D^{-1}}. This makes  
#' rank-deficient systems numerically stable.
#'
#' @param xp      n points \eqn{x_i} of dimension d are arranged in (n x d) matrix \code{xp}
#' @param U       vector of length n, containing samples \eqn{u(x_i)} of 
#'                the scalar function \eqn{u} to be fitted \cr
#'                - or - \cr
#'                (n x m) matrix, where each column 1,...,m contains one vector of samples
#'                \eqn{u_j(x_i)} for the m'th model, j=1,...,m
#' @param squares [FALSE] flag, see 'Description'
#' @param ptail   [TRUE] flag, see description
#' @param width   [-1] either a positive real value which is the constant width \eqn{\sigma} for all 
#'                Gaussians in all iterations, or -1. If -1, the appropriate width \eqn{\sigma} is 
#'                calculated anew in each iteration with one of the rules \code{RULE},
#'                based on the distribution of data points \code{xp}.              
#' @param RULE    ["One"] one out of ["One" | "Two" | "Three"], different rules for automatic 
#'                estimation of width \eqn{\sigma}. Only relevant if \code{width = -1},   
#' @param widthFactor [1.0] additional constant factor applied to each width \eqn{\sigma} 
#' @param rho     [0.0] experimental: 0.0: interpolating, >0.0, approximating (spline-like) 
#'                Gaussian RBFs
#' @param DEBUG2  [FALSE] if TRUE, save \code{M} and \code{rhs} on return value
#' @return \code{rbf.model},  an object of class \code{RBFinter}, which is basically a list 
#' with elements:
#'      \item{coef}{  (n+d+1 x m) matrix holding in column m the coefficients for the m'th 
#'                    model:      \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d}.  
#'                    In case \code{squares==TRUE} it is an (n+2d+1 x m) matrix holding  
#'                    additionally the coefficients \eqn{c_{d+1}, ..., c_{d+d}}.}
#'      \item{xp}{    matrix xp   }
#'      \item{d}{ size of the polynomial tail. If \code{length(d)==0} it means no polynomial tail will be used for the model. In case of ptail==T && squares==F d will be dimension+1 and in case of ptail==T && squares==T d will be 2*dimension+1 }
#'      \item{npts}{  number n of points \eqn{x_i} }
#'      \item{ptail}{TRUE or FALSE (see description)}
#'      \item{squares}{  TRUE or FALSE (see description)}
#'      \item{width}{ the calculated width \eqn{\sigma} }
#'      \item{type}{  "MQ"}
#'      
#' @seealso   \code{\link{trainCubicRBF}}, \code{\link{predict.RBFinter}}, \code{\link{interpRBF}}
#' @author Wolfgang Konen, Samineh Bagheri
#'          
#----------------------------------------------------------------------------------
trainMQRBF <- function(xp, U, ptail=TRUE, squares=FALSE, 
                       width,RULE="One",widthFactor=1.0, rho=0.0, DEBUG2=F) {
  if(width==-1){#When no width is given by user then width is selected automatically
    #RULE<-"One"
    #Automatic width adaptation is in testing phase therefore we would like to test several rules
    #from literature and then select teh most relevant ones
    
    #>> rule number One:
    #The rule is taken from Benoudjit,2002 (But originally from Haykin,1999)
    #It works based on  finding a compromise between locality and smoothness
    
    #>> rule number Two:
    #The rule is taken from Benoudjit,2002 (But originally from Moody and Darken,1989)
    #a vector is returned as width and width[i] is the width factor for the ith RBF
    
    #>> rule number Three:
    # This rule is taken from Sun and Jin,2013
    # The idea is basically adapting the width to the smallest interval of two points in each coordiante
    
    switch(RULE,
           "One" =     width<-(max(stats::dist(xp)))/sqrt(2*nrow(xp)),
           "Two" = {
             k<-2
             width<-sapply(1:nrow(xp),function(i){
               centroid<-matrix(xp[i,],nrow=1)
               ncentroid<-xp[-i,]
               #browser()
               #index<-FNN::knnx.index(ncentroid,centroid, k=k) #finding the index of k nearest neigbors
               mindists<-FNN::knnx.dist(ncentroid,centroid, k=k)
               #mindists<-mindists^2
               y<-(1/k)*(sqrt(sum(mindists*mindists)))
               return(y)
             })},
           "Three" = {width<-c()
           for(i in 1:ncol(xp)){
             interval<-max(xp[,i])-min(xp[,i])
             width<-min(width,interval)  
             
           }}
    )
    verboseprint(verbose=0,important=FALSE,paste("Automatic adjustment of RBF width=",width[1]))
    
  }
  

  edist=as.matrix(stats::dist(xp,upper=T))         # euclidean distance matrix

  width = width*widthFactor
  if (length(width)==1) {                   # /WK/ bug fix: former code was not correct for the case 
    wmat=width^2                            # that width may be a vector of length nrow(xp)
  } else {
    w = outer(rep(1,length(width)),width)
    wmat = w*w                              # this choice leads to non-symmetric phi, but it 
    # has smallest error on RbfFit2.R (multiply each width by 
    # factor 3 - 5)
    #wmat = w*t(w)                          # this is (wmat)_ij = width_i * width_j, it
    # guarantees that phi will be a symmetric matrix,
    # but it has bigger error  
  }
   phi <- sqrt(1+edist*edist/wmat)           # multiquadric RBF

  rbf.model<-trainRBF(phi,U,ptail,squares,xp,type="MQ",DEBUG2,width=width,rho=rho)
  
  return(rbf.model)
}




#----------------------------------------------------------------------------------
#' Apply the trained cubic, MQ or Gaussian RBF interpolation to new data for d>1.
#' 
#' @param x         vector holding a point of dimension d
#' @param rbf.model trained RBF model (or set of models), see \code{\link{trainCubicRBF}} 
#'                  or \code{\link{trainGaussRBF}}
#'                
#' @return          value \eqn{s(\vec{x})} of the trained model at \eqn{\vec{x}} \cr
#'                  - or - \cr
#'                  vector \eqn{s_j( \vec{x})} with values for all trained models \eqn{j=1,...,m} at \eqn{\vec{x}}
#' 
#' @seealso   \code{\link{trainCubicRBF}},  \code{\link{trainMQRBF}},  \code{\link{trainGaussRBF}}, \code{\link{predict.RBFinter}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@th-koeln.de})
#----------------------------------------------------------------------------------
interpRBF <- function(x,rbf.model) {
  #testit::assert("non-conform",length(x)==ncol(rbf.model$xp))
  testit::assert("rbf.model is not of class RBFinter",any(class(rbf.model)=="RBFinter"))
  if (length(x)!=ncol(rbf.model$xp)) {
    cat("problem in interpRBF\n")
    #browser()
    stop("Problem in interpRBF")
  }
  ed = distLine(x,rbf.model$xp)# euclidean distance of x to all rbf.model$xp, 
                                      # this is up to 40x faster than dist() !!
                                      # ed is a vector of length nrow(xp)
  
  val<-c()
  if(length(rbf.model$width)>1){
    for(counter in c(1:length(rbf.model$width))){
      cwidth<-rbf.model$width[counter]
      type<-rbf.model$type[counter]
      switch(type,
             "CUBIC" = {
               ph = ed*ed*ed           
             },
             "GAUSS" = {
               ph = exp(-0.5*(ed/cwidth)^2)  # this works correctly even in case where width is a vector of length nrow(xp)
             },
             "MQ" = {
               ph =sqrt((ed^2)/(cwidth^2)+1)
             })

      if (rbf.model$ptail) {
        if (rbf.model$squares) 
          lhs = c(ph,1,x,x*x)
        else 
          lhs = c(ph,1,x)
      } else {
        lhs = ph
      }  
      val = c(val,lhs %*% rbf.model$coef[,counter] )
    }
    
  }else{
    switch(rbf.model$type,
           "CUBIC" = {
             ph = ed*ed*ed           
           },
           "GAUSS" = {
             ph = exp(-0.5*(ed/rbf.model$width)^2)  # this works correctly even in case where width is a vector of length nrow(xp)
           },
           "MQ" = {
             ph =sqrt((ed^2)/(rbf.model$width^2)+1)
           })
    
    if (rbf.model$ptail) {
      if (rbf.model$squares) 
        lhs = c(ph,1,x,x*x)
      else 
        lhs = c(ph,1,x)
    } else {
      lhs = ph
    }   
    #browser()
    val = as.vector(lhs %*% rbf.model$coef) 
  }
  
  return (val)
}


#----------------------------------------------------------------------------------
#' Apply cubic or Gaussian or MQ RBF interpolation
#' 
#' Apply cubic or Gaussian or MQ RBF interpolation to a set of new data points for d>1.
#' 
#' @param rbf.model trained RBF model (or set of models), see \code{\link{trainCubicRBF}} 
#'                  or \code{\link{trainGaussRBF}}
#' @param newdata   matrix or data frame with d columns. Each row contains a data point 
#'                  \eqn{x_i,\ i=1,\ldots,n}
#' @param ...       (not used)
#'                
#' @return          vector of model responses \eqn{s(x_i)}, one element for each data point \eqn{x_i} \cr
#'                  - or - \cr
#'                  if \code{rbf.model} is a set of \code{m} models, a \code{(n x m)}-matrix 
#'                  containing in each row the response \eqn{s_j(x_i)} of all models 
#'                  \eqn{j = 1,\ldots,m}  to \eqn{x_i}
#'  
#' @seealso   \code{\link{trainCubicRBF}}, \code{\link{trainGaussRBF}}, \code{\link{interpRBF}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@th-koeln.de})
#----------------------------------------------------------------------------------
predict.RBFinter <- function(rbf.model,newdata,...) {
  val = t(sapply(1:nrow(newdata),function(i)interpRBF(as.numeric(newdata[i,]),rbf.model)))
  # sapply binds the result vectors from interpRBF() by default **column-wise** together.
  # Therefore we use t() to return in the ith *row* of val 
  # the results for the ith row in newdata (row-wise).
 
#### 
#### This alternative version, which does the job of interpRBF right here (in-place),
#### is not really faster than the one sapply-line above (!)#
####
#   x = as.matrix(newdata)
#   nx = nrow(x)
#   ed = t(sapply(1:nx,function(i)(dist(rbind(x[i,],rbf.model$xp)))[1:rbf.model$npts] ))
#   ph = ed*ed*ed
#   if (rbf.model$squares) 
#     lhs = cbind(ph,rep(1,nx),x,x*x)
#   else 
#     lhs = cbind(ph,rep(1,nx),x)
#   val = as.vector(lhs %*% rbf.model$coef)
  
  # if there is only one model in rbf.model:
  if (ncol(rbf.model$coef)==1) val = as.vector(val)
  return(val)
}
