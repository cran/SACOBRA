debugModel<-function(type="RMSE",nPoints=1000,cobra){
 # browser()
 dimension<-ncol(cobra$A)
 RandomPoints<-t(as.matrix(sapply(1:nPoints, FUN=function(i)stats::runif(dimension,-1,1))))
 realEval<-apply(RandomPoints,1,cobra$fn)
 realEvalObj<-realEval[1,]
 fitnessSurrogate<-cobra$fitnessSurrogate
 modelEval<-apply(RandomPoints,1,getPredY0,fitnessSurrogate=fitnessSurrogate,cobra=cobra)
 RMSE<-sum((realEvalObj-modelEval)^2)/length(modelEval)
# print(paste("RMSE (obj):",RMSE))
 
 realEvalCon<-realEval[-1,]
 constraintSurrogates<-cobra$constraintSurrogates
# browser()
 modelCon<-apply(RandomPoints,1,interpRBF,rbf.model=constraintSurrogates)
 diff<-abs(modelCon-realEvalCon)
 if(cobra$nConstraints==1){
   RMSEC<-sqrt(mean(sum(diff^2)))
 }else{
   RMSEC<-apply(diff,1,function(x){
     sqrt(mean( sum(x^2)))
   })  
 }
 
# browser()
 print(RMSEC)
 cobra$RMSE<-c(cobra$RMSE,RMSE)
 cobra$RMSEC<-rbind(cobra$RMSEC,RMSEC)
 return(cobra)
}


defaultMODELD<-function(){
  DM<-list(active=FALSE,
           type="RMSE",
           nPoints=200,
           freq=10)
  return(DM)  
}