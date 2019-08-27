######################################################################################
# multiRunPlot_2
#
#' Plot the results from multipe COBRA runs.
#'
#' Plot for each run one black curve 'error vs. iterations' and aggregate the mean curve (red) and 
#' the median curve (green) of all runs. 
#' DIFFERENCE to \code{\link{multiRunPlot}}: 'error' is the distance of the ever-best feasible point 
#' in input space to the true solution \code{solu}.
#'
#' Print some diagnostic information: final median & mean error, percentage of runs which meet the
#' target (only if \code{optim} is available)). \cr
#'
#' @param dfAll      the data frame of all runs, obtained with \code{\link{multiCOBRA}} or loaded 
#'                   from .Rdata file
#' @param solu       the true solution in input space of the problem (only for diagnostics).
#' @param fName      ["multiRun"] the name of the .Rdata file, printed as subtitle
#' @param main       [""] the name of the problem (e.g. "G01 problem"), printed as title  
#' @param xlim       the x limits
#' @param ylim       the y limits
#' @param ylog       [TRUE]  logarithmic y-axis
#' @param xlog       [FALSE] logarithmic x-axis
#' @param target     [0.05] a single run meets the target, if the final error is smaller than \code{target}
#' @param plotPDF    [FALSE] if TRUE, plot to \code{<fName>.pdf}
#' @param subPDF     [NULL] optional subdirectory where .pdf should go
#' @param legendWhere  ["topright"]
#' @param absErr     [FALSE] if TRUE, plot abs(error) instead of error.
#'
#' @return \code{z3}, a vector containing for each run the ever-best feasible objective value
#'
#' @seealso   \code{\link{multiRunPlot}}, \code{\link{multiCOBRA}}, \code{\link{cobraPhaseII}}
#' @author Wolfgang Konen, Samineh Bagheri, Cologne University of Applied Sciences
#' @export
#' 
multiRunPlot_2 <- function(dfAll,solu,fName="multiRun",main="",
                           xlim=NULL,ylim=c(1e-05,1e+04),ylog=TRUE, 
                           xlog=FALSE,target=0.05, plotPDF=FALSE, subPDF=NULL,
                           legendWhere="topright",absErr=FALSE ) 
{
  if (absErr) {
    absFunc = abs
  } else {
    absFunc = function(x)x;
  }
  runList = unique(dfAll$run)
  ffcMin = min(dfAll$ffc)
  ffcMax = max(dfAll$ffc)
  ffcChk = sapply(runList,function(i)dfAll$ffc[max(which(dfAll$run==i))])
  if (any(ffcChk!=ffcMax)) {
    warning(sprintf("Not all runs terminate after same func evals: %s",paste(ffcChk,collapse=",")))
    ffcMax = min(ffcChk)    # only symptomatic cure for z3 below
  }
  ylab="dist to solu"
  if (is.null(xlim)) xlim = c(ffcMin,ffcMax)
  df = dfAll[dfAll$run==runList[1],]
  slog = ""
  if (xlog) slog=paste(slog,"x",sep="")
  if (ylog) slog=paste(slog,"y",sep="")
  pdfname = paste("results.d",sub(".Rdata",".pdf",fName),sep="/")
  if (!is.null(subPDF)) pdfname=paste(subPDF,pdfname,sep="/")
  if (plotPDF) {
    grDevices::pdf(pdfname,width=7,height=5.24)
    graphics::par(cex.axis=1.3)
    graphics::par(cex.lab=1.3)
    graphics::par(cex.sub=0.9)
  } else {
    graphics::par(cex.axis=1.0)
    graphics::par(cex.lab=1.0)
    graphics::par(cex.sub=0.7)    
  }
  
  if (is.matrix(solu)) {
    dimension=ncol(solu)
  } else {
    dimension=length(solu)
  }
  
  mk_Xbest <- function(df,dimension) {
    dfXbest = df[,(ncol(df)-dimension+1):ncol(df)]
    Best <- df$everBestFeas
    for (k in 1:nrow(df)) {
      if (!is.na(Best[k])) {
        if (!is.na(Best[k-1])) {
          if (Best[k]==Best[k-1]) dfXbest[k,]=dfXbest[k-1,]
        }
      }
    }
    return(dfXbest)
  }
  #browser()
  
  distAll = NULL
  for (i in runList) {
    df = dfAll[dfAll$run==runList[i],]
    
    # set dfAll$everBestFeas of each run to NA until the first feasible iterate arrives 
    ind <- which(df$feas==TRUE)
    if (length(ind)>0) {
      ind=min(ind)
      if (ind>1 & ind!=Inf)
        df$everBestFeas[1:(ind-1)] <- NA;   
      dfAll$everBestFeas[dfAll$run==runList[i]] <- df$everBestFeas;
    } else {
      cat(sprintf("NOTE: No feasible solution in run %d\n",i))
      dfAll$everBestFeas[dfAll$run==runList[i]] <- NA;      
    }   

    dfXbest <- mk_Xbest(df,dimension)
    if (is.matrix(solu)) {
      distS = NULL
      for (k in 1:nrow(solu)) {
        distSolu<- (dfXbest - outer(rep(1,nrow(df)),solu[k,]))^2
        distSolu<- sqrt(rowSums(distSolu))
        distS <- rbind(distS,distSolu)
        distSolu <- apply(distS,2,min)
      }
    } else {
      distSolu<- (dfXbest - outer(rep(1,nrow(df)),solu))^2
      distSolu<- sqrt(rowSums(distSolu))
    }
    distAll <- rbind(distAll,data.frame(distSolu=distSolu,run=rep(i,nrow(df))))    
    
    if (i==runList[1]) {
      graphics::plot(df$ffc,distSolu,log=slog,type="l",main=main,sub=fName
                     ,xlim=xlim,ylim=ylim,ylab=ylab,xlab="iterations") 
    }
    df$distSolu <- distSolu
    graphics::lines(df$ffc,df$distSolu)
  }
  
  graphics::legend(legendWhere,c("mean","median"),lwd=2,col=c("red","green"))
  #browser()
  indInfeas = which(dfAll$feas==FALSE)
  z = stats::aggregate(distAll$distSolu,list(dfAll$run),min, na.rm=TRUE, na.action=NULL)
  cat(sprintf("Avg solution: %10.5f +- %7.5f,  true minimum: %10.5f\n", 
              mean(z$x),sd(z$x),0.0))  #,"\n"
  
  z = stats::aggregate(distAll$distSolu,list(dfAll$ffc),mean,na.rm=T)
  names(z) <- c("ffc","distSolu")
  z = z[which(z$ffc<=ffcMax),]            # strip off ffc==41, if there
  graphics::lines(z$ffc,z$distSolu,col="red",lwd=2)
  cat("Final avg error (mean):  ", z$distSolu[nrow(z)],"\n")  
  
  z2 = stats::aggregate(distAll$distSolu,list(dfAll$ffc),stats::median,na.rm=T)
  names(z2) <- c("ffc","distSolu")
  z2 = z2[which(z2$ffc<=ffcMax),]            # strip off ffc==41, if there
  graphics::lines(z2$ffc,z2$distSolu,col="green",lwd=2)
  cat("Final avg error (median):", z2$distSolu[nrow(z)],"\n")  
  
  z3 = distAll$distSolu[dfAll$ffc==ffcMax]
  pTarget = length(which(z3 < target))/length(z3)
  cat("Target (",target,") reached with probability: ",pTarget,"\n")
  
  if (plotPDF) {
    grDevices::dev.off()
    cat("Plot saved to",pdfname,"\n")
  }
  return(dfXbest)
}
