# splitMaStrings
#

splitMaStrings<-function(x){
  y<-list()
  for(i in 1:length(x)){
    if(is.na(pmatch("Gaussian",x[i]))){
      
      temp<-strsplit(x[i],"MQ")
      temp[[1]][1]<-"MQ"
      if(is.na(pmatch("MQ",x[i]))){
        temp[[1]][1]<-"cubic"
      }
      
    }else{
      temp<-strsplit(x[i],"Gaussian")
      temp[[1]][1]<-"Gaussian"
      
    }
    y[i]<-temp
  }
  return(y)
}

#' Training surrogates for objective and constraint functions 
#' 
#' This function is called during the second phase of SACOBRA \code{\link{cobraPhaseII}} to train the surrogate models for objective and constarint functions 
#' of the optimization problem passed to the SACOBRA framework. This function takes all the so-far evalauted points stored in \code{cobra$A} and builds
#' 
#'
#' @param cobra an object of class COBRA, this is a (long) list containing all settings
#'        from \code{\link{cobraPhaseII}}
#'
#' @return \code{cobra}, an object of class COBRA, 
#'    with the following elements modified:
#'      \item{\code{constraintSurrogates}}{  surrogate models for all constraints }
#'      \item{\code{fitnessSurrogate}}{  surrogate model for objective function  }
#'      \item{\code{fitnessSurrogate1}}{ surrogate model for objective function,  w/o \code{\link{plog}} }
#'      \item{\code{fitnessSurrogate2}}{ surrogate model for objective function, with \code{\link{plog}} }
#'
#' @keywords internal   
trainSurrogates <- function(cobra) {
  verboseprint(cobra$verbose,important=FALSE,paste(">> Training" ,cobra$RBFmodel,"surrogates","..."))
  ptm <- proc.time()
  cobra$Fres <- as.vector(cobra$Fres)
  A<-cobra$A <- as.matrix(cobra$A)

  #added option of adaptive plog
  ind <- 1:length(cobra$Fres) #SB: do we need this parameter at all?
  if(cobra$sac$aFF){
    # print("adjusting fitness function")
    cobra<-adFit(cobra,ind)
    Fres<-cobra$SurrogateInput
  }else{
    Fres=cobra$Fres[ind]
  }
  
  if(cobra$TFlag){
    Fres<-cobra$TFres
    A<-cobra$TA
  }
  
  if(cobra$DOSAC >0){
    if(cobra$PLOG[length(cobra$PLOG)] && cobra$printP){
      verboseprint(cobra$verbose, important=TRUE,paste("PLOG transformation is done ( iter=",nrow(cobra$A),")"))
    }   
  } 
  cobra$printP<-FALSE
  
  
  
  if (cobra$CONSTRAINED) {
    cobra$Gres <- as.matrix(cobra$Gres)
    Gres=cobra$Gres[ind,]  
    if(cobra$MS$apply==F || cobra$MS$active==F){
      sw=switch(cobra$RBFmodel,
                "cubic" =     {cobra$constraintSurrogates <- trainCubicRBF(A,Gres,ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)
                cobra$fitnessSurrogate <- trainCubicRBF(A,Fres,ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)},
                "Gaussian"=   {cobra$constraintSurrogates <- trainGaussRBF(A,Gres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor);
                cobra$fitnessSurrogate <- trainGaussRBF(A,Fres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor)},
                "MQ" = {cobra$constraintSurrogates <- trainMQRBF(A,Gres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor);
                cobra$fitnessSurrogate <- trainMQRBF(A,Fres,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor)},
                "InvalidRBFmodel"
      )
    }else{
      
      lastSelectedModel<-cobra$selectedModel[nrow(cobra$selectedModel),]
      lastSelectedModel<-splitMaStrings(lastSelectedModel)
      cons<-length(lastSelectedModel)-1
      #objective model
      sw=switch(lastSelectedModel[[1]][1],
             "cubic"={fitnessModel<-"cubic";
             cobra$fitnessSurrogate <- trainCubicRBF(A,Fres,ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)},
             "Gaussian"={fitnessModel<-"Gauss";
             fwidth<-as.numeric(lastSelectedModel[[1]][2]);
             cobra$fitnessSurrogate <- trainGaussRBF(A,Fres,1,
                                                     ptail=cobra$ptail,squares=cobra$squares,
                                                     RULE=cobra$RULE,
                                                     rho=cobra$RBFrho,
                                                     widthFactor=fwidth)},
             "MQ"={fitnessModel<-"MQ";
             fwidth<-as.numeric(lastSelectedModel[[1]][2]);
             cobra$fitnessSurrogate <- trainMQRBF(A,Fres,1,
                                                  ptail=cobra$ptail,squares=cobra$squares,
                                                  RULE=cobra$RULE,
                                                  rho=cobra$RBFrho,
                                                  widthFactor=fwidth)})
      
      
      COEF<-NULL
      types<-c()
      widths<-c()
      Gres<-as.matrix(Gres)
      for(con in c(1:cons)){
        
        sw=switch(lastSelectedModel[[con+1]][1],
               "cubic"={conModel<-"CUBIC";
               model <- trainCubicRBF(A,Gres[,con],ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho);
               COEF<-cbind(COEF,model$coef);
               cwidth<-0},
               "MQ"={conModel<-"MQ";
               cwidth<-as.numeric(lastSelectedModel[[con+1]][2]);
               model <- trainMQRBF(A,Gres[,con],1,
                                   ptail=cobra$ptail,squares=cobra$squares,
                                   RULE=cobra$RULE,
                                   rho=cobra$RBFrho,
                                   widthFactor=cwidth);
               COEF<-cbind(COEF,model$coef)},
               "Gaussian"={conModel<-"GAUSS";
               cwidth<-as.numeric(lastSelectedModel[[con+1]][2]);
               model <- trainGaussRBF(A,Gres[,con],1,
                                      ptail=cobra$ptail,squares=cobra$squares,
                                      RULE=cobra$RULE,
                                      rho=cobra$RBFrho,
                                      widthFactor=cwidth);
               COEF<-cbind(COEF,model$coef)}
        )
        
        types<-c(types,conModel)
        widths<-c(widths,cwidth)
      }     
      
      rbf.model<-list(coef=COEF
                      ,xp=A
                      ,d=ncol(A)
                      ,npts=nrow(A)
                      ,squares=T
                      ,ptail=cobra$ptail
                      ,width=widths
                      ,type=types
                      #,AUGMENTED=model$AUGMENTED
                      
      )  
      class(rbf.model) <- c("RBFinter","list")
      cobra$constraintSurrogates<-rbf.model
      
    }

  } #end of if cobra$CONSTRAINED
  
  if(!cobra$CONSTRAINED){
    
    if(cobra$MS$apply==T && cobra$MS$active==T){
      lastSelectedModel<-cobra$selectedModel[nrow(cobra$selectedModel),]
     # print(lastSelectedModel)
      lastSelectedModel<-splitMaStrings(lastSelectedModel)
      objectiveKernel<- lastSelectedModel[[1]][1]
      objw<-as.numeric(lastSelectedModel[[1]][2])
    }else{
      objectiveKernel<-cobra$RBFmodel
      objw<-cobra$widthFactor
    }
    #objective model
    sw=switch(objectiveKernel,
              "cubic"={fitnessModel<-"cubic";
              cobra$fitnessSurrogate <- trainCubicRBF(A,Fres,ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)},
              "Gaussian"={fitnessModel<-"Gauss";
              cobra$fitnessSurrogate <- trainGaussRBF(A,Fres,1,
                                                      ptail=cobra$ptail,squares=cobra$squares,
                                                      RULE=cobra$RULE,
                                                      rho=cobra$RBFrho,
                                                      widthFactor=objw)},
              "MQ"={fitnessModel<-"MQ";
              cobra$fitnessSurrogate <- trainMQRBF(A,Fres,1,
                                                   ptail=cobra$ptail,squares=cobra$squares,
                                                   RULE=cobra$RULE,
                                                   rho=cobra$RBFrho,
                                                   widthFactor=objw)})
  }
  
  if (cobra$DEBUG_RBF$active){ 
    print(nrow(A))
    cobra <- debugVisualizeRBF(cobra,cobra$fitnessSurrogate,A,Fres)        # see defaultDebugRBF.R
  }

  #SB: added the possibilty to measure p-effect after every 10 iterations
  if((cobra$sac$onlinePLOG && nrow(cobra$A)%%cobra$sac$onlineFreqPLOG==0 )|| nrow(cobra$A)==cobra$initDesPoints){
    
   # Fres1<-cobra$Fres             #without plog
   # Fres2<-sapply(cobra$Fres,plog)#with plog
    
    #adapted to the hessian calculation
    Fres1<-Fres             #without plog
    Fres2<-sapply(Fres,plog)#with plog
    
    #two models are built after each onlineFreqPLOG iterations:
    #                                            fitnessSurrogate1-> 
    sw=switch(cobra$RBFmodel,
              "cubic" =     {
                cobra$fitnessSurrogate1 <- trainCubicRBF(A,Fres1,ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)
                cobra$fitnessSurrogate2 <- trainCubicRBF(A,Fres2,ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)},
              "Gaussian"=   {
                cobra$fitnessSurrogate1 <- trainGaussRBF(A,Fres1,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor)
                cobra$fitnessSurrogate2 <- trainGaussRBF(A,Fres2,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor)},
              "MQ" = {cobra$fitnessSurrogate1 <- trainMQRBF(A,Fres1,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor);
              cobra$fitnessSurrogate2 <- trainMQRBF(A,Fres2,cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,rho=cobra$RBFrho,widthFactor=cobra$widthFactor)},
              "InvalidRBFmodel"
    ) 
    
  }
  testit::assert(sprintf("Wrong value %s for cobra$RBFmodel",sw),sw!="InvalidRBFmodel")
  DO_ASSERT=F
  if (DO_ASSERT) {
    #might need adjust due to rescale /WK/  
    conFunc <- {function(x)cobra$fn(x)[-1];}
    Gres = t(sapply(1:nrow(cobra$A),function(i){conFunc(cobra$A[i,])}))
    testit::assert("Gres-assertion failed",all(Gres==cobra$Gres))
    testit::assert("cobra$A-assertion failed",all(cobra$A==cobra$constraintSurrogates$xp))
    Gres = t(sapply(1:nrow(cobra$A),function(i){interpRBF(cobra$A[i,],cobra$constraintSurrogates)}))
    for (i in 1:ncol(cobra$Gres)) {
      z = (Gres[,i]-cobra$Gres[,i]) / (max(Gres[,i])-min(Gres[,i]))
      if (max(abs((z))>1e-5)) {
        verboseprint(cobra$verbose, important=FALSE,paste("interpRBF(..,cobra$constraintSurrogates)-assertion failed for constraint",i,""))
        print(max(abs(z))) #do we need this?
      }        
    }
    cat("All assertions passed\n")
  }
  
  verboseprint(cobra$verbose, important=FALSE,paste(" finished (",(proc.time() - ptm)[3],"sec )"))
  return(cobra);
  
} # trainSurrogates()


