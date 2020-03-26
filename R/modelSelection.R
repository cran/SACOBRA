


#' Default settings for the model-selection part of SACOBRA.
#'
#' Sets default values for the model-selection part \code{cobra$MS} of SACOBRA.  \cr
#' It is shown that different types of RBFs can deliver different qualites in modeling different functions. 
#' Using the online model selection functionality boosted the overall performace of SACOBRA on a large set of constrained problems.
#' The algorithm trains every function (objective and constraints) with a given pool of models including different RBF types and width parameters.
#' The type of model which performs the best in the last iterations \code{WinS} will be selected for each function. 
#' The quality of the models are determined by different measures of approximation error in each iteration
#' \deqn{f(\vec{x}_{new})-s(\vec{x}_{new})}
#' 
#
#Details:
#' With the call \code{\link{setOpts}(MS,defaultMS())} it is possible to extend a partial list 
#' \code{MS} to a list containing all \code{MS}-elements (the missing ones are taken from 
#' \code{defaultMS()}).\cr \cr
#' \strong{NOTE:} Because of common crash observation, it is not recommended to include Gaussian model in the set of models especially for problems which require more than 100 function evaluations.
#'
#' @return MS,  a list with the following elements 
#'      \item{active}{[F] If set to TRUE then \code{selectModel} calculates the best model for each constraint(s)/objective function }
#'      \item{models}{[c("cubic","MQ")] a set of model types that will be used to build the pool of models. Three types of RBF are implemneted "cubic", "Gaussian" and "MQ" (multiquadric). 
#' Users can select one or combination of these models. Users can select a set of "width parameters" for "MQ" and "Gaussian" by setting \code{widths} parameter. }
#'      \item{widths}{[c(0.01,0.1,1,10)] a set of values for width parameter of RBF models. Only relevant if \code{models} include "Gaussian" or "MQ". } 
#'      \item{freq}{[1] controls how often \code{selectModel} is called. 
#'      In every \code{freq} iterations all the selected models are trained for all constraint/objective function(s)}
#'      \item{slidingW}{[T] when set to FALSE it uses the information taken from all the past iterations to asses the quality of the models. 
#'      When set to TRUE, activates the sliding window functionality and it takes the information of the last \code{WinS} iterations (see \code{WinS}).}
#'      \item{WinS}{[1] size of the sliding window }
#'      \item{quant}{[3] 3: median, 2:0.25, 4:0.75. The measure used to compare the quality of the model in the last window.}
#'      \item{apply}{[T] if set to FALSE then the selected models are not used during the optimization. Only for debugging purposes. }
#'      \item{considerXI}{[F] If set to T then a subset of the approximation errors which are related to the current (DRC element) 
#'      are considered to make the model selection decision}
#' @seealso   \code{\link{setOpts}}
#' @author Samineh Bagheri 
#' @export
defaultMS<-function(){
  MS<-list(  active=FALSE,
             models=c("cubic","MQ"),
             widths=c(0.01,0.1,1,10),
             freq=1,
             slidingW=F,
             WinS=10,
             quant=3, #3: median, 2: 25%, 4:75%
             apply=T,
             considerXI=F

  )
  return(MS)
}


# selectModel is ...
selectModel<-function(models=c("cubic","Gaussian","MQ"),
                      widths=c(0.01,0.1,1,10),cobra,freq=1,
                      slidingW,WinS,quant,gama){

  cobra$MSXI<-c(cobra$MSXI,gama)
  if(cobra$MS$considerXI){
    acceptedIndices<-which(cobra$MSXI==gama)
  }else{
    acceptedIndices<-c(1:length(cobra$MSXI))
  }
  
  num<-nrow(cobra$A)
  mycolNames<-c("model","objective")
  Fres<-cobra$Fres
  A<-cobra$A
  test<-A[num,]
  trainPoints<-c(1:(num-1)) 
  testRes<-Fres[num]
  if(cobra$CONSTRAINED){
    mycolNames<-c(mycolNames,paste("g",c(1:ncol(cobra$Gres)),sep=""))
    Gres<-cobra$Gres
    testRes<-c(Fres[num],Gres[num,])
  }
  
    

  

  counter0<-0
  
    if(any(models=="cubic")){
      counter0<-counter0+1
      fitModel<-c()
      
      fitModel <- trainCubicRBF(A[trainPoints,],Fres[trainPoints],ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)
      fitPred <- interpRBF(test,fitModel)
      pred<-fitPred
      
      if(cobra$CONSTRAINED){
        conModel<-c()
        conModel <- trainCubicRBF(A[trainPoints,],Gres[trainPoints,],ptail=cobra$ptail,squares=cobra$squares,rho=cobra$RBFrho)
        conPred <- interpRBF(test,conModel) 
        pred<-c(pred,conPred)
      }
      
      
      
      newApproxErr<-abs(testRes-pred)
      cobra$ApproxFrame$cubic<-rbind(cobra$ApproxFrame$cubic,newApproxErr)
    }
counter1<-counter0
      if(any(models=="Gaussian")){
      
      for(width in widths){
        counter1<-counter1+1
        newApproxErr<-NULL
        conModel<-c()
        fitModel<-c()
        if(cobra$RBFwidth<0)
          W<-cobra$RBFwidth
        else 
          W<-1
        
        fitModel <- trainGaussRBF(A[trainPoints,],Fres[trainPoints],W,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,
                                  rho=cobra$RBFrho,widthFactor=width) 
        
        fitPred <- interpRBF(test,fitModel)
        pred<-fitPred
        
        if(cobra$CONSTRAINED){
          conModel <- trainGaussRBF(A[trainPoints,],Gres[trainPoints,],W,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,
                                  rho=cobra$RBFrho,widthFactor=width)
          conPred <- interpRBF(test,conModel)
          pred<-c(pred,conPred)
          
        }
        newApproxErr<-abs(testRes-pred)
        #index<-which(width==widths)+1
        newApproxErrF<-t(data.frame(newApproxErr))
        if(length(cobra$ApproxFrame)<counter1){
          cobra$ApproxFrame<-c(cobra$ApproxFrame,list(newApproxErrF))
        }else{
          cobra$ApproxFrame[[counter1]]<-rbind(cobra$ApproxFrame[[counter1]],newApproxErrF)
          
        }
      }
      
      names(cobra$ApproxFrame)[c((counter0+1):counter1)]<-paste("Gauss",widths)
    }

  counter<-counter1
  if(any(models=="MQ")){
    for(width in widths){
      counter<-counter+1
      newApproxErr<-NULL
      conModel<-c()
      fitModel<-c()

      fitModel <- trainMQRBF(A[trainPoints,],Fres[trainPoints],cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,
                                rho=cobra$RBFrho,widthFactor=width) 
      fitPred <- interpRBF(test,fitModel)
      pred<-fitPred

      if(cobra$CONSTRAINED){
        conModel <- trainMQRBF(A[trainPoints,],Gres[trainPoints,],cobra$RBFwidth,ptail=cobra$ptail,squares=cobra$squares,RULE=cobra$RULE,
                               rho=cobra$RBFrho,widthFactor=width) 
        conPred <- interpRBF(test,conModel)
        pred<-c(pred,conPred)
        
      }
      
      
      newApproxErr<-abs(testRes-pred)
      
      index<-which(width==widths)+1
      newApproxErrF<-t(data.frame(newApproxErr))
      if(length(cobra$ApproxFrame)<counter){
        cobra$ApproxFrame<-c(cobra$ApproxFrame,list(newApproxErrF))
      }else{
        cobra$ApproxFrame[[counter]]<-rbind(cobra$ApproxFrame[[counter]],newApproxErrF)
        
      }
    }
    names(cobra$ApproxFrame)[c((counter1+1):counter)]<-paste("MQ",widths)
  }
  
  #frameR<-nrow(cobra$ApproxFrame[[1]][acceptedIndices,])
  frameR<-length(acceptedIndices)
  MedianFrame<-NULL
  
  if(any(models=="cubic")){
  if(slidingW && frameR>WinS){
       MedianFrame<-data.frame(model="cubic",t(apply(data.frame(data.frame(cobra$ApproxFrame$cubic[acceptedIndices,])[c(max(frameR-WinS,1):frameR),]),2,FUN=function(x){
         return(quantile(x)[quant])
         })))
  }else{
    if(is.matrix(cobra$ApproxFrame$cubic[acceptedIndices,]))
      tempFrame<-cobra$ApproxFrame$cubic[acceptedIndices,]
    
    else
      tempFrame<-data.frame(matrix(cobra$ApproxFrame$cubic[acceptedIndices,],nrow=length(acceptedIndices),ncol=cobra$nConstraints+1))
      #tempFrame<-data.frame(t(cobra$ApproxFrame$cubic[acceptedIndices,]))
    
    MedianFrame<-data.frame(model="cubic",t(apply(tempFrame,2,FUN=function(x){
      return(quantile(x)[quant])
    })))
    
  }
    names(MedianFrame)<-mycolNames
    
  }

  

  if(any(models=="Gaussian")){ 
    for(width in widths){
    #index<-which(width==widths)+1
    index<-which(names(cobra$ApproxFrame)==paste("Gauss",width))
    
    if(slidingW && frameR>WinS){
      MedianFrame<-rbind(MedianFrame,
                         data.frame(model=paste("Gaussian",width,sep=""),
                                    t(apply(
                                      data.frame(data.frame(cobra$ApproxFrame[[index]][acceptedIndices,])[c(max(frameR-WinS,1):frameR),])
                                      ,2,FUN=function(x){
                                        return(quantile(x)[quant])
                                      })))
      ) 
    }else{
      
      if(is.matrix(cobra$ApproxFrame[[index]][acceptedIndices,]))
        tempFrame<-cobra$ApproxFrame[[index]][acceptedIndices,]
      else
        tempFrame<-data.frame(matrix(cobra$ApproxFrame$cubic[acceptedIndices,],nrow=length(acceptedIndices),ncol=cobra$nConstraints+1))
      #  tempFrame<-data.frame(t(cobra$ApproxFrame[[index]][acceptedIndices,]))
      
      MedianFrame<-rbind(MedianFrame,
                         data.frame(model=paste("Gaussian",width,sep=""),t(apply(tempFrame,2,FUN=function(x){
                           return(quantile(x)[quant])
                         })))
      )   
    }
  }
  }
  

  if(any(models=="MQ")){
    for(width in widths){
     # index<-which(width==widths)+1
      index<-which(names(cobra$ApproxFrame)==paste("MQ",width))
      if(slidingW && frameR>WinS){
        newTemp<-data.frame(model=paste("MQ",width,sep=""),
                            t(apply(
                              data.frame(data.frame(cobra$ApproxFrame[[index]][acceptedIndices,])[c(max(frameR-WinS,1):frameR),])
                              ,2,FUN=function(x){
                                return(quantile(x)[quant])
                              })))
        names(newTemp)<-mycolNames
        MedianFrame<-rbind(MedianFrame,newTemp) 
      }else{
        
        if(is.matrix(cobra$ApproxFrame[[index]][acceptedIndices,]))
          tempFrame<-cobra$ApproxFrame[[index]][acceptedIndices,]
        else
          tempFrame<-data.frame(matrix(cobra$ApproxFrame$cubic[acceptedIndices,],nrow=length(acceptedIndices),ncol=cobra$nConstraints+1))
          #tempFrame<-data.frame(t(cobra$ApproxFrame[[index]][acceptedIndices,]))
        
        
        newTemp<-data.frame(model=paste("MQ",width,sep=""),t(apply(tempFrame,2,FUN=function(x){
          return(quantile(x)[quant])
        })))
        names(newTemp)<-mycolNames

        MedianFrame<-rbind(MedianFrame,newTemp)   
      }
    } 
  }
  ModelName<-as.character(MedianFrame$model)
  
  #MedianFrame is a n times m data.frame, where n is the number of models in the pool of models and m is 2+nConstarints. 
  #The first columnn constains strings (The name of the RBF types) the second column is the approximate error for the objective function
  #The rest of columsn are the representing the approximate error for constraint functions
  testit::assert("MedianFrame dimension is wrong:",ncol(MedianFrame)==cobra$nConstraints+2)
  
  
  cobra$selectedModel<-rbind(cobra$selectedModel,apply(data.frame(MedianFrame[,-1]),2,function(x){
      index<-which(x==min(x))
      index<-index[1]
      y<-ModelName[index]
       return(y)
             }))
  #cobra$selectedModel is dataframe with 1+nConstarints columns
  print(paste(cobra$selectedModel[nrow(cobra$selectedModel),]))
  return(cobra)
}