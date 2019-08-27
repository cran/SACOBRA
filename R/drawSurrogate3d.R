##
## only for debugging purposes: visualize RBF surrogate (only d=2)
## 
## caller: defaultDebugRBF.R (from cobraPhaseII.R --> trainSurrogates.R, switch DEBUG_RBF)
##    and  trustRegion.R, switch DEBUG_TRU
##
drawSurrogate3d <- function(surrogate,fName,X,Y,Xr,Yr,lower,upper, 
                            TrueZ=NULL,overlayTrueZ=FALSE,
                            add=FALSE,newWindow=TRUE,device=NULL,title="",A){
  if (!is.null(surrogate)) {
    N = nrow(X)
    
    
    if(class(surrogate)!="function"){
      ZS = sapply(1:N,function(j)predict(surrogate,cbind(Xr[,j],Yr[,j])))
    }else{
      #browser()
      ZS = sapply(1:N,function(j)apply(cbind(Xr[,j],Yr[,j]),1,surrogate))
    }
   
   
    jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    graphics::filled.contour(X[1,],Y[,1],t(ZS),main=sprintf("%s, %s surrogate, N=%d",fName,title,nrow(A)),
                             xlab="x", ylab="y", color.palette=jet.colors)
    # if (newWindow) {
    #    rgl::open3d()      
    # } else {
    #   if (!add) rgl::rgl.clear()
    # }
    devNum<-device
    if(is.null(device) || newWindow){
      devNum<-rgl::open3d()
    }else{
      rgl::rgl.set(device)
    }
    if (!add) rgl::rgl.clear() 
    
    
    rngZS=max(ZS,TrueZ)-min(ZS,TrueZ)
    rngX = upper[1] - lower[1]
    rngY = upper[2] - lower[2]
    rngZ = max(ZS)-min(ZS)
    rgl::aspect3d(1/rngX,1/rngY,1/rngZ)
    drawSurface3d(X,Y,ZS)
    if (overlayTrueZ & !is.null(TrueZ))
      drawSurface3d(X,Y,TrueZ)
    rgl::axes3d(color="black") # draw axes and box with tickmarks
   # rgl::title3d(sprintf("%s, %s surrogate, N=%d",fName,title,surrogate$npts), '', 'x', 'y', 'z',color="black")
    rgl::rgl.bg(color="white")
    rgl::rgl.bringtotop(stay=F)
    return(list(ZS=ZS,devNum=devNum))
  }  
}

drawSurface3d <- function(X,Y,Z) {
  jet.colors <-   
    grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  colorzjet <- jet.colors(100)  # 100 separate colors 
  rgl::rgl.surface(x=X, y=Y, z=Z,
                   coords=c(1,2,3), 
                   color=colorzjet[ findInterval(Z, seq(min(Z), max(Z), length=100))]
  )        
}

snapShot3d <- function(fitnessSurrogate,every,prefix,A) {
  if (nrow(A) %% every == 0) {
    #if (fitnessSurrogate$npts==6) browser()  # optional pause to adjust 3d-orientation
    if (!file.exists("images.d")) dir.create("images.d")
    filename <- sprintf("images.d/%s-%03d.png",prefix, nrow(A))
    rgl::rgl.snapshot(filename)       
  }
}

