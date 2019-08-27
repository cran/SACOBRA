#defaultSAC.R
#' 
#' Default settings for the SACOBRA part of SACOBRA.
#' 
#' Sets suitable defaults for the SACOBRA part of SACOBRA. \cr
#' With the call \code{\link{setOpts}(mySAC,defaultSAC())} it is possible to extend a partial list 
#' \code{mySAC} to a list containing all \code{sac}-elements (the missing ones are taken from 
#' \code{defaultSAC()}).
#' 
#' For backward compatibility, a logical DOSAC (deprecated) is mapped from FALSE to 0 and 
#' from TRUE to 1.
#' 
#' @param DOSAC  [0|1|2] with default 1.\cr
#'                0: COBRA-R settings (turn off SACOBRA), \cr
#'                1: SACOBRA settings, \cr
#'                2: SACOBRA settings with fewer parameters and more online adjustements (aFF and aCF are done parameter free).
#'                        
#'
#' @return a list with the following elements (the values in parantheses \code{[|]} are the values 
#'    for \code{DOSAC=[0|1|2]}):
#'      \item{RS}{  flag for random start algorithm \code{[FALSE|TRUE|TRUE]}  }
#'      \item{RStype}{type of the function to calculate probability to start the internal optimizer 
#'          with a random starting point\code{[NA|"SIGMOID"|"CONSTANT"]} (see function 
#'          \code{\link{RandomStart}} in \code{SACOBRA.R})}
#'      \item{RSmax}{maximum probability of a random start when \code{RStype=="SIGMOID"} 
#'          (see \code{\link{RandomStart}} in \code{SACOBRA.R}). If \code{RStype=="CONSTANT"} then  
#'          random start is done with a constant probability determined from mean(c(RSmax,RSmin)) 
#'          \code{[NA|0.3|0.3]}}
#'      \item{RSmin}{minimum probability of a random start when \code{RStype=="SIGMOID"}  
#'          (see \code{\link{RandomStart}} in \code{SACOBRA.R}) \code{[NA|0.05|0.05]} }
#'      \item{RSAUTO}{If TRUE then in every iteration where the fraction of feasible points in the population 
#'          is smaller than 0.05, the RS probability is set to 0.3. \code{[FALSE|FALSE|TRUE]}  }
#'      \item{aDRC}{  flag for automatic DRC adjustment \code{[FALSE|TRUE|TRUE]}  }
#'      \item{aFF}{  flag for automatic objective function transformation \code{[FALSE|TRUE|TRUE]}  }
#'      \item{aCF}{  flag for automatic constraint function transformation \code{[FALSE|TRUE|TRUE]}  }
#'      \item{TFRange}{  threshold, if \code{FRange} is larger than \code{TFRange}, then apply automatic objective
#'          function transformation (see \code{\link{plog}}). \code{[Inf|1e+05|-1]}  }
#'      \item{TGR}{  threshold, if \code{GRatio} is larger than \code{TGR}, then apply automatic 
#'          constraint function transformation. \code{GRatio} is the ratio "largest \code{GRange} /
#'          smallest \code{GRange}" where \code{GRange} is the min-max range of a specific constraint.
#'          If \code{TGR < 1}, then the transformation is always performed.  \code{[Inf|1e+03|-1]}.}
#'      \item{Cs}{  If  \code{Cs} iterations in a row do not improve the ever-best feasible solution, 
#'          then perform a random restart. \code{[10|10|10]}  }
#'      \item{adaptivePLOG}{  (experimental) flag for objective function transformation with \code{\link{plog}}, 
#'          where the parameter \code{pShift} is adapted during iterations. 
#'          \code{[FALSE|FALSE|FALSE]} }
#'      \item{onlinePLOG}{ flag for online decision making wether use plog or not according to p-effect \code{\link{plog}}. 
#'          \code{[FALSE|FALSE|TRUE]} }
#'      \item{pEffectInit}{Initial pEffect value when using onlinePLOG. If pEffectInit >= 2 then the initial model is built after plog transformation.
#'          \code{[NA|NA|2]}}
#'
#' @seealso   \code{\link{cobraInit}}, \code{\link{cobraPhaseII}}
#'
#' @author Samineh Bagheri, Cologne University of Applied Sciences
#' @export
#'
defaultSAC<-function(DOSAC=1){
  if (is.logical(DOSAC)) {
    if(DOSAC==TRUE) DOSAC=1
    if(DOSAC==FALSE) DOSAC=0
  }
  if (!(DOSAC %in% c(0,1,2))) 
    stop("DOSAC has to be either 0, 1 or 2.")
  
  if(DOSAC==1){
    sac<-list(RS=TRUE,
              RStype="SIGMOID", # "CONSTANT",  "SIGMOID"
              RSAUTO=FALSE,
              RSmax=0.3, #maximum probability of a random start
              RSmin=0.05, #minimum probability of a random start
              aDRC=TRUE,
              aFF=TRUE,
              TFRange=1e+05,
              aCF=TRUE,
              TGR=1e+03,   #GR threshold to perform constraint modification, 
              #if TGR < 1 transformation is always performed,
              # Inf value means transfomration is never performed
              # a positive value laarger than 1 means: 
              #the transfomartion is only performed for problems with a GR larger than TGR
              conPLOG=FALSE, #perform plog transformation for all constraints (In testing phase) 
              conFitPLOG=FALSE, #perform plog transformation for all constraints and fitness (In testing phase)
              adaptivePLOG=FALSE, #the effectivity of adaptivePLOG is not fully proved therefore I keep it as FALSE for now
              onlinePLOG=FALSE,
              onlineFreqPLOG=10, # number of iterations in a row which after that the online DOPLOG check is done
              pEffectInit=0,
              minMaxNormal=F,
              onlineMinMax=F,
              Cs=10)
  }
  if(DOSAC==2){
    sac<-list(RS=TRUE,
              RStype="CONSTANT", #"CONSTANT",  "SIGMOID"
              RSAUTO=TRUE,
              RSmax=0.3, #maximum probability of a random start
              RSmin=0.05, #minimum probability of a random start
              aDRC=TRUE,
              aFF=TRUE,
              TFRange=-1,
              aCF=TRUE,
              TGR=-1,   #GR threshold to perform constraint modification, 
              #if TGR < 1 transformation is always performed,
              # Inf value means transfomration is never performed
              # a positive value laarger than 1 means: 
              #the transfomartion is only performed for problems with a GR larger than TGR
              conPLOG=FALSE, #perform plog transformation for all constraints (In testing phase) 
              conFitPLOG=FALSE, #perform plog transformation for all constraints and fitness (In testing phase)
              adaptivePLOG=FALSE, #the effectivity of adaptivePLOG is not fully proved therefore I keep it as FALSE for now
              onlinePLOG=TRUE,
              onlineFreqPLOG=10, # number of iterations in a row which after that the online DOPLOG check is done
              pEffectInit=3,
              minMaxNormal=F,
              onlineMinMax=F,
              Cs=10)
  }
  
  if(DOSAC==0){
    sac<-list(RS=FALSE,
              RStype="SIGMOID", #"CONSTANT",  "SIGMOID"
              RSAUTO=FALSE,
              RSmax=0.3, #maximum probability of a random start
              RSmin=0.05, #minimum probability of a random start
              aDRC=FALSE,
              aFF=FALSE,
              TFRange=Inf,
              aCF=FALSE,
              TGR=Inf,   #GR threshold to perform constraint modification, 
              #if TGR < 1 transformation is always performed,
              # Inf value means transfomration is never performed
              # a positive value laarger than 1 means: 
              #the transfomartion is only performed for problems with a GR larger than TGR 
              conPLOG=FALSE, #perform plog transformation for all constraints (In testing phase)
              conFitPLOG=FALSE, #perform plog transformation for all constraints and fitness (In testing phase)
              adaptivePLOG=FALSE,
              onlinePLOG=FALSE,
              onlineFreqPLOG=10,
              pEffectInit=0,
              minMaxNormal=F,
              onlineMinMax=F,
              Cs=10)
  }
return(sac)  
}

# setOpts
#' 
#' Merge the options from a partial list and the default list
#'
#' @param opts       a partial list of options
#' @param defaultOpt a list with default values for every element
#'
#' @return a list combined from \code{opts} and \code{defaultOpt} where every available element 
#'    in \code{opts} overrides the default. For the rest of the elements the value from \code{defaultOpt}
#'    is taken. \cr
#'    A warning is issued for every element appearing in \code{opts} but not in \code{defaultOpt}
#'    
#' @seealso   \code{\link{defaultRI}}, \code{\link{defaultSAC}}, \code{\link{defaultTR}}, \code{\link{defaultEquMu}}
#'
#' @author Samineh Bagheri, Wolfgang Konen, Cologne University of Applied Sciences
#' @export
#'
setOpts<-function(opts,defaultOpt){
  
  setting<-defaultOpt
  
  if(methods::hasArg(opts)){
    matching<-intersect(names(opts),names(defaultOpt))
    setting[matching] <- opts[matching]
    
    notMatching <- setdiff(names(opts),
                           names(defaultOpt))
    if(length(notMatching)!=0) warning(paste("The
            following arguments are ignored: ", notMatching))
    
  }
  
  
  return(setting)
}