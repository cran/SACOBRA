

#' Constraint Optimization Problem Benchmark (G Function Suite)
#' 
#' COP is an object of class \code{R6ClassGenerator} which can be used to access G problems (aka G functions) implementations in R, by simply generating a new instance of 
#' COP for each G function \code{problem<-COP.new("problem")}. The COP instances have the following useful attributes:\cr
#' \itemize{
#' \item name : name of the problem given by the user
#' \item dimension: dimension of the problem. For the scalable problems \code{G02} and \code{G03}, the dimension should be given by users, otherwise it will be set automaticaly
#' \item lower: lower boundary of the problem
#' \item upper: upper boundary of the problem
#' \item fn: the COP function which can be passed to SACOBRA. (see \code{fn} description in \code{\link{cobraInit}})
#' \item nConstraints: number of constraints
#' \item xStart: The suggested optimization starting point
#' \item solu: the best known solution, (only for diagnostics purposes)
#' \item info: information about the problem}
#' G function suite is a set of 24 constrained optimization problems with various properties like 
#' dimension, number of equality/ inequality constraint, feasibilty ratio, etc. 
#' Although these problems were introduced as a suite in a technical report at CEC 2006, many of them have been used by different
#' autors earlier.
#' \cr For more details see:
#' Liang, J., Runarsson, T.P., Mezura-Montes, E., Clerc, M., Suganthan, P.,
#' Coello, C.C., Deb, K.: Problem definitions and evaluation criteria for the CEC
#' 2006 special session on constrained real-parameter optimization. Journal of
#' Applied Mechanics 41, 8 (2006), \url{http://www.lania.mx/~emezura/util/files/tr_cec06.pdf}
#' 
#'   
#' @examples 
#' ##creating an instance for G24 problem
#' G24<-COP$new("G24")
#' 
#' ##initializing SACOBRA
#' cobra <- cobraInit(xStart=G24$lower, fName=G24$name,
#'                    fn=G24$fn,  
#'                    lower=G24$lower, upper=G24$upper, feval=25)
#'                    
#' ## Run sacobra optimizer
#' cobra <- cobraPhaseII(cobra)
#' 
#' ## The true solution is at solu = G24$solu
#' ## The solution found by SACOBRA:
#' print(getXbest(cobra))
#' print(getFbest(cobra))                    
#' plot(abs(cobra$df$Best-G24$fn(G24$solu)[1]),log="y",type="l",
#' ylab="error",xlab="iteration",main=G24$name)
#' 
#' 
#' ## creating an instance for G03 in 2-dimensional space
#' G03<-COP$new("G03",2)
#' 
#' ## Initializing sacobra
#'cobra <- cobraInit(xStart=G03$lower, fn=G03$fn, 
#'fName=G03$name, lower=G03$lower, upper=G03$upper, feval=40)
#' 
#' @export
#' @author Samineh Bagheri, Wolfgang Konen
#' @keywords datasets
COP<-R6::R6Class("COP",
                 public=list(name=NA,
                             dimension=NA,
                             lower=NA,
                             upper=NA,
                             fn=NA,
                             nConstraints=NA,
                             xStart=NA,
                             solu=NA,
                             info=NA,
                             initialize=function(name,dimension){
                               allProbs=sprintf("G%02i",c(1:24))
                               if (!missing(name)){ 
                                 self$name <- name
                                 if(!any(name==allProbs))stop(paste("We cannot load the problem you chose, please choose one the following COPs:",capture.output(cat(allProbs, sep=", "))))
                               }else{
                                 stop(paste("Please choose one the following COPs:",capture.output(cat(allProbs, sep=", "))))
                               }
                               
                               if (!missing(dimension)){ 
                                 self$dimension <- dimension
                               }else{
                               }
                             switch(name,
                                    "G01"={private$callG01()},
                                    "G02"={private$callG02(self$dimension)},
                                    "G03"={private$callG03(self$dimension)},
                                    "G04"={private$callG04()},
                                    "G05"={private$callG05()},
                                    "G06"={private$callG06()},
                                    "G07"={private$callG07()},
                                    "G08"={private$callG08()},
                                    "G09"={private$callG09()},
                                    "G10"={private$callG10()},
                                    "G11"={private$callG11()},
                                    "G12"={private$callG12()},
                                    "G13"={private$callG13()},
                                    "G14"={private$callG14()},
                                    "G15"={private$callG15()},
                                    "G16"={private$callG16()},
                                    "G17"={private$callG17()},
                                    "G18"={private$callG18()},
                                    "G19"={private$callG19()},                                    
                                    "G20"={private$callG20()},
                                    "G21"={private$callG21()},
                                    "G22"={private$callG22()},
                                    "G23"={private$callG23()},
                                    "G24"={private$callG24()}) 
                             }
                             ),
                 private = list(                             
                   callG01=function(){
                     self$fn=function(x){
                       
                       obj<- sum(5*x[1:4])-(5*sum(x[1:4]*x[1:4]))-(sum(x[5:13]))
                       g1<- (2*x[1]+2*x[2]+x[10]+x[11] - 10)
                       g2<- (2*x[1]+2*x[3]+x[10]+x[12] - 10)
                       g3<- (2*x[2]+2*x[3]+x[11]+x[12] - 10)
                       
                       g4<- -8*x[1]+x[10]
                       g5<- -8*x[2]+x[11]
                       g6<- -8*x[3]+x[12]
                       
                       g7<- -2*x[4]-x[5]+x[10]
                       g8<- -2*x[6]-x[7]+x[11]
                       g9<- -2*x[8]-x[9]+x[12]
                       
                       res<-c(obj, g1 ,g2 , g3 
                              , g4 , g5 , g6 
                              , g7 , g8 , g9)
                       return(res)
                       
                     };
                     self$dimension=13   
                     self$lower=rep(0,self$dimension)
                     self$upper=c(rep(1,9),rep(100,3),1)
                     self$nConstraints=9
                     self$solu=c(rep(1,9),rep(3,3),1)  
                     self$xStart=rep(0,self$dimension) 
                   },
                   callG02=function(d){
                     d<-self$dimension
                     self$fn=function(x){
                       d<-self$dimension
                       den<-c()
                       for(i in 1:d){
                         den<-c(den,i*(x[i]^2))
                       }
                       obj<- -abs( (sum(cos(x)^4)-(2*prod(cos(x)^2)))/
                                     sqrt(sum(den))
                       )
                       g1<- 0.75-prod(x)
                       g2<- (sum(x)-7.5*d)
                       
                       res<-c(obj, g1 ,g2)
                       return(res)
                     };
                     self$lower=rep(1e-16,d);
                     self$upper=rep(10,d);
                     self$nConstraints=2;
                     if(d==2) self$solu=c(1.600859,0.4684985)
                     if(d==10) self$solu=c(3.1238477, 3.0690696, 3.0139085, 
                                           2.9572856, 1.4654789, 0.3684877, 
                                           0.3633289, 0.3592627, 0.3547453,
                                           0.3510025 )
                     
                   },
                   callG03=function(d){
                     d<-self$dimension
                     self$fn=function(x){
                       obj<--((sqrt(d))^d)*prod(x)
                       g1<-sum(x*x)-1
                       res<-c(obj=obj,equ=g1)
                       return(res)  
                     };
                     self$lower=rep(0,d);
                     self$upper=rep(1,d);
                     self$nConstraints=1; 
                     self$solu=rep(1/sqrt(d),d)
                   },
                   callG04=function(){
                     self$fn=function(x){
                       obj<-(5.3578547*(x[3]^2))+(0.8356891*x[1]*x[5])+(37.293239*x[1])-40792.141
                       g1<- -(85.334407+0.0056858*x[2]*x[5]+0.0006262*x[1]*x[4]-0.0022053*x[3]*x[5]) 
                       g2<- 85.334407+0.0056858*x[2]*x[5]+0.0006262*x[1]*x[4]-0.0022053*x[3]*x[5]-92
                       
                       
                       g3<-  90-(80.51249+(0.0071317*x[2]*x[5])+(0.0029955*x[1]*x[2])+(0.0021813*x[3]^2))
                       g4<- 80.51249+(0.0071317*x[2]*x[5])+(0.0029955*x[1]*x[2])+(0.0021813*x[3]^2)-110
                       
                       g5<-  20-(9.300961+(0.0047026*x[3]*x[5])+(0.0012547*x[1]*x[3])+(0.0019085*x[3]*x[4]))
                       g6<-  9.300961+(0.0047026*x[3]*x[5])+(0.0012547*x[1]*x[3])+(0.0019085*x[3]*x[4])-25
                       
                       
                       res<-c(objective=obj, g1=g1 , g2=g2 , g3=g3, g4=g4,g5=g5 , g6=g6)
                       return(res)
                       
                     };
                     self$dimension=5   
                     self$lower=c(78,33,rep(27,3))
                     self$upper=c(102,rep(45,4))
                     self$nConstraints=6
                     self$solu=c( 78.00000000000000000,
                                  33.00000000000000000,
                                  29.99525602568341398,
                                  44.99999999847009491,
                                  36.77581290640451783)
                     self$xStart=c(runif(1,min=78,max=102), #x1
                                   runif(1,min=33,max=45),  #x2
                                   runif(3,min=27,max=45)) 
                   },
                   callG05=function(){
                     self$fn=function(x){
                       
                       obj<-3*x[1]+(10^(-6))*(x[1]^3)+2*x[2]+((2*10^(-6))/3)*(x[2]^3)
                       
                       g1<- x[3]-x[4]-0.55
                       g2<- x[4]-x[3]-0.55
                       
                       g3<- 1000*sin(-x[3]-0.25)+1000*sin(-x[4]-0.25)+894.8-x[1] 
                       g4<- 1000*sin(x[3]-0.25)+1000*sin(x[3]-x[4]-0.25)+894.8-x[2] 
                       g5<- 1000*sin(x[4]-0.25)+1000*sin(x[4]-x[3]-0.25)+1294.8  
                       
                       res<-c(objective=obj , equ=g3, equ=g4,equ=g5 , g=g1 , g=g2)
                       
                     };
                     self$dimension=4   
                     self$lower=c(rep(0,2),rep(-0.55,2))
                     self$upper=c(rep(1200,2),rep(0.55,2))
                     self$nConstraints=5
                     self$solu=solu <- c(679.94531748791177961,
                                         1026.06713513571594376,
                                         0.11887636617838561,
                                         -0.39623355240329272)  
                     self$xStart=c(runif(2,min=0,max=1200), #x1 , x2
                                   runif(2,min=-0.55,max=0.55)) #x3,x4
                   },
                   callG06=function(){
                     self$fn=function(x){
                       
                       obj<- ((x[1]-10)^3 ) + ((x[2]-20)^3)
                       g1<- -((((x[1]-5)^2)+((x[2]-5)^2)-100)) 
                       g2<-  (((x[1]-6)^2)+((x[2]-5)^2)-82.81)
                       res<-c(objective=obj, g1=g1 , g2=g2 )
                       return(res)
                       
                     };
                     self$dimension=2  
                     self$lower=c(13,0)
                     self$upper=rep(100,self$dimension)
                     self$nConstraints=2
                     self$solu=c(14.095,(5 - sqrt(100-(14.095-5)^2)))
                     self$xStart=c(20.1,5.84) 
                   },
                   callG07=function(){
                     self$fn=function(x){
                         obj<- (x[1]^2)+(x[2]^2)+(x[1]*x[2])-(14*x[1])-(16*x[2])+((x[3]-10)^2)+(4*((x[4]-5)^2))+ ((x[5]-3)^2) +( 2*((x[6]-1)^2)) +(5*(x[7]^2))+( 7*((x[8]-11)^2))+(2*((x[9]-10)^2 ))+ ((x[10]-7)^2) + 45;
                         
                         g1<-((4.0*x[1])+(5.0*x[2])-(3.0*x[7])+(9.0*x[8]) - 105.0);
                         g2<-(10*x[1]-8*x[2]-17*x[7]+2*x[8]);
                         g3<-(-8*x[1]+2*x[2]+5*x[9]-2*x[10]-12);
                         g4<-(3*((x[1]-2)^2)+4*(x[2]-3)^2+2*x[3]^2-7*x[4]-120);
                         g5<-(5*x[1]^2+8*x[2]+(x[3]-6)^2-2*x[4]-40);
                         g6<-(x[1]^2+2*(x[2]-2)^2-(2*x[1]*x[2])+14*x[5]-6*x[6]);
                         g7<-(0.5*(x[1]-8)^2+2*(x[2]-4)^2+3*(x[5]^2)-x[6]-30);
                         g8<-(-3*x[1]+6*x[2]+12*(x[9]-8)^2-7*x[10]);
                         
                         res<-c(objective=obj, g1=g1 , g2=g2 , g3=g3,g4=g4,
                                g5=g5,g6=g6,g7=g7,g8=g8)
                         return(res)
                       
                     };
                     self$dimension=10  
                     self$lower=rep(-10,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=8
                     self$solu=c(  2.171997834812, 
                                   2.363679362798, 
                                   8.773925117415, 
                                   5.095984215855,
                                   0.990655966387, 
                                   1.430578427576, 
                                   1.321647038816, 
                                   9.828728107011,
                                   8.280094195305, 
                                   8.375923511901 );
                     self$xStart=runif(self$dimension,min=self$lower,max=self$upper)  
                   },
                   callG08=function(){
                     self$fn=function(x){
                       obj<- ((-sin(2*pi*x[1])^3)*(sin(2*pi*x[2])))/((x[1]^3)*(x[1]+x[2]));
                       
                       g1<- x[1]^2-x[2]+1
                       g2<- 1-x[1]+(x[2]-4)^2
                       
                       res<-c(objective=obj, g1=g1 , g2=g2 )
                       return(res)
                       
                     };
                     self$dimension=2  
                     self$lower=rep(0.00001,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=2
                     self$solu=c(1.2279713,4.2453733)  
                     self$xStart=NA 
                   },
                   callG09=function(){
                     self$fn=function(x){
                       
                       obj<-(x[1]-10)^2+5*(x[2]-12)^2+x[3]^4+3*(x[4]-11)^2+10*(x[5]^6)+7*x[6]^2 +x[7]^4-4*x[6]*x[7]-10*x[6]-8*x[7] ;
                       g1<- (2*x[1]^2+3*x[2]^4+x[3]+4*x[4]^2+5*x[5]-127)
                       g2<- (7*x[1]+3*x[2]+10*x[3]^2+x[4]-x[5]-282)
                       g3<- (23*x[1]+x[2]^2+6*x[6]^2-8*x[7]-196)
                       g4<- 4*x[1]^2+x[2]^2-3*x[1]*x[2]+2*x[3]^2+5*x[6]-11*x[7]
                       
                       res<-c(objective=obj, g1=g1 , g2=g2, g3=g3, g4=g4 )
                       return(res)
                       
                     };
                     self$dimension=7  
                     self$lower=rep(-10,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=4
                     self$solu=solu <- c(2.33049949323300210,
                                         1.95137239646596039,
                                         -0.47754041766198602,
                                         4.36572612852776931,
                                         -0.62448707583702823,
                                         1.03813092302119347,
                                         1.59422663221959926);
                     self$xStart=NA
                   },
                   callG10=function(){
                     self$fn=function(x){
                       
                       obj<- x[1]+x[2]+x[3];
                       
                       g1<- (-1+0.0025*(x[4]+x[6]))
                       g2<- (-1 + 0.0025*(-x[4]+x[5]+x[7]))
                       g3<- (-1+0.01*(-x[5]+x[8]))
                       g4<- (100*x[1]-(x[1]*x[6])+833.33252*x[4]-83333.333)
                       g5<- (x[2]*x[4]-x[2]*x[7]-1250*x[4]+1250*x[5])
                       g6<- (x[3]*x[5]-x[3]*x[8]-2500*x[5]+1250000)
                       
                       res<-c(objective=obj, g1=g1 , g2=g2, g3=g3, g4=g4, g5=g5, g6=g6)
                       return(res)
                       
                     };
                     self$dimension=8 
                     self$lower=c(100,1000,1000,rep(10,5))
                     self$upper=c(rep(10000,3),rep(1000,5))
                     self$nConstraints=6
                     self$solu=c(579.29340269759155,
                                 1359.97691009458777,
                                 5109.97770901501008,
                                 182.01659025342749,
                                 295.60089166064103,
                                 217.98340973906758,
                                 286.41569858295981,
                                 395.60089165381908); 
                     self$xStart=c(rep(1001,3),rep(100,self$dimension-3))  
                   },
                   callG11=function(){
                     self$fn=function(x){
                       
                       y<-x[1]*x[1]
                       y<-y+((x[2]-1)^2)
                       y<-as.numeric(y)
                       g1 <- as.numeric(+(x[2] - x[1]^2))
                       res<-c(objective=y, equ=g1)
                       return(res)
                       
                     };
                     self$dimension=2  
                     self$lower=c(-1,-1)
                     self$upper=c(1,1)
                     self$nConstraints=1
                     self$solu=c(-sqrt(0.5),0.5)
                     self$xStart=NA
                   },
                   callG12=function(){
                     self$fn=function(x){
                       obj <- -1+0.01*((x[1]-5)^2+(x[2]-5)^2+(x[3]-5)^2);
                       G<-c()
                       for( i in c(1:9)){
                         for(j in c(1:9)){
                           for(k in c(1:9)){
                             G<-c(G,g<- (x[1]-i)^2+(x[2]-j)^2+(x[3]-k)^2 - 0.0625)
                           }
                         }
                       }
                       # these are 729 disjoint spheres and a solution is feasible if it is within *one* of the 
                       # 729 spheres. Therefore we take the min over G:
                       #
                       
                       res<-c(obj,min(G))
                       names(res)<-c("objective","G")
                       return(res)
                       
                     };
                     self$dimension=3   
                     self$lower=rep(0,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=1
                     self$solu=c(5,5,5)
                     self$xStart= NA
                   },
                   callG13=function(){
                     self$fn=function(x){
                       obj<-exp(x[1]*x[2]*x[3]*x[4]*x[5])
                       g1<-x[1]^2+x[2]^2+x[3]^2+x[4]^2+x[5]^2-10
                       g2<-x[2]*x[3]-5*x[4]*x[5]
                       g3<-x[1]^3+x[2]^3+1
                       res<-c(objective=obj, equ=g1 , equ=g2, equ=g3)
                       
                       return(res)
                       
                     };
                     self$dimension=5  
                     self$lower=c(rep(-2.3,2),rep(-3.2,3))
                     self$upper=c(rep(2.3,2),rep(3.2,3))
                     self$nConstraints=3
                     solu0<-c(-1.7171435947203,1.5957097321519,1.8272456947885,-0.7636422812896,-0.7636439027742)
                     solu<-solu0
                     solu<-rbind(solu,c(solu0[1:2],-solu0[3],-solu0[4],+solu0[5]))
                     solu<-rbind(solu,c(solu0[1:2],-solu0[3],+solu0[4],-solu0[5]))
                     solu<-rbind(solu,c(solu0[1:2],+solu0[3],-solu0[4],-solu0[5]))
                     solu<-rbind(solu,c(solu0[1:2],-solu0[3],+solu0[4],-solu0[5]))
                     solu<-rbind(solu,c(solu0[1:2],-solu0[3],-solu0[4],+solu0[5]))
                     self$solu=solu
                     self$xStart=NA 
                     self$info="Please note that G13 has multiple global optima, all stored in solu"
                   },
                   callG14=function(){
                     self$fn=function(x){
                       cVector<-c(-6.089,-17.164,-34.054,-5.914,-24.721,-14.986,-24.1,-10.708,-26.662,-22.179)
                       sumX<-sum(x)
                       y<-log(x/sumX)
                       obj<-sum(x*(cVector+y))
                       g1<-x[1]+2*x[2]+2*x[3]+x[6]+x[10]-2
                       g2<-x[4]+2*x[5]+x[6]+x[7]-1
                       g3<-x[3]+x[7]+x[8]+2*x[9]+x[10]-1
                       
                       res<-c(objective=obj, equ=g1 , equ=g2, equ=g3)
                       return(res)
                       
                     };
                     self$dimension=10 
                     self$lower=rep(1e-06,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=3
                     self$solu=c(0.0406684113216282, 0.147721240492452, 0.783205732104114,
                                 0.00141433931889084, 0.485293636780388, 0.000693183051556082,
                                 0.0274052040687766,
                                 0.0179509660214818, 0.0373268186859717, 0.0968844604336845);
                     self$xStart=NA 
                   },
                   callG15=function(){
                     self$fn=function(x){
                       obj<-1000-(x[1]^2)-2*x[2]^2-x[3]^2-x[1]*x[2]-x[1]*x[3]
                       
                       h1<-x[1]^2+x[2]^2+x[3]^2-25
                       h2<-8*x[1]+14*x[2]+7*x[3]-56
                       
                       res<-c(objective=obj, equ=h1 , equ=h2)
                       return(res)
                       
                     };
                     self$dimension=3   
                     self$lower=rep(0,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=2
                     self$solu=c(3.51212812611795133,0.216987510429556135, 3.55217854929179921)  
                     self$xStart=NA  
                   },
                   callG16=function(){
                     self$fn=function(x){
                       y1<-x[2]+x[3]+41.6
                       c1<-0.024*x[4]-4.62
                       y2<-(12.5/c1)+12
                       c2<-0.0003535*(x[1]^2)+0.5311*x[1]+0.08705*y2*x[1]
                       c3<-0.052*x[1]+78+0.002377*y2*x[1]
                       y3<-c2/c3
                       y4<-19*y3
                       c4<-0.04782*(x[1]-y3)+(0.1956*(x[1]-y3)^2)/x[2]+0.6376*y4+1.594*y3
                       c5<-100*x[2]
                       c6<-x[1]-y3-y4
                       c7<-0.950-(c4/c5)
                       y5<-c6*c7
                       y6<-x[1]-y5-y4-y3
                       c8<-(y5+y4)*0.995
                       y7<-c8/y1
                       y8<-c8/3798
                       c9<-y7-(0.0663*(y7/y8))-0.3153
                       y9<-(96.82/c9)+0.321*y1
                       y10<-1.29*y5+1.258*y4+2.29*y3+1.71*y6
                       y11<-1.71*x[1]-0.452*y4+0.580*y3
                       c10<-12.3/752.3
                       c11<-(1.75*y2)*(0.995*x[1])
                       c12<-(0.995*y10)+1998
                       y12<-c10*x[1]+(c11/c12)
                       y13<-c12-1.75*y2
                       y14<-3623+64.4*x[2]+58.4*x[3]+146312/(y9+x[5])
                       c13<-0.995*y10+60.8*x[2]+48*x[4]-0.1121*y14-5095
                       y15<-y13/c13
                       y16<-148000-331000*y15+40*y13-61*y15*y13
                       c14<-2324*y10-28740000*y2
                       y17<-14130000-(1328*y10)-(531*y11)+(c14/c12)
                       c15<-(y13/y15)-(y13/0.52)
                       c16<-1.104-0.72*y15
                       c17<-y9+x[5]
                       
                       
                       obj<-(0.000117*y14)+0.1365+(0.00002358*y13)+(0.000001502*y16)+(0.0321*y12)+
                         (0.004324*y5)+(0.0001*c15/c16)+(37.48*(y2/c12))-(0.0000005843*y17)
                       
                       
                       g1<-(0.28/0.72)*y5-y4
                       g2<-x[3]-1.5*x[2]
                       g3<-3496*(y2/c12)-21
                       g4<-110.6+y1-(62212/c17)
                       g5<-213.1-y1
                       g6<-y1-405.23
                       g7<-17.505-y2
                       g8<-y2-1053.6667
                       g9<-11.275-y3
                       g10<-y3-35.03
                       g11<-214.228-y4
                       g12<-y4-665.585
                       g13<-7.458-y5
                       g14<-y5-584.463
                       g15<-0.961-y6
                       g16<-y6-265.916
                       g17<-1.612-y7
                       g18<-y7-7.046
                       g19<-0.146-y8
                       g20<-y8-0.222
                       g21<-107.99-y9
                       g22<-y9-273.366
                       g23<-922.693-y10
                       g24<-y10-1286.105
                       g25<-926.832-y11
                       g26<-y11-1444.046
                       g27<-18.766-y12
                       g28<-y12-537.141
                       g29<-1072.163-y13
                       g30<-y13-3247.039
                       g31<-8961.448-y14
                       g32<-y14-26844.086
                       g33<-0.063-y15
                       g34<-y15-0.386
                       g35<-71084.33-y16
                       g36<--140000+y16
                       g37<-2802713-y17
                       g38<-y17-12146108
                       
                       
                       res<-c(obj=obj
                              ,g1=g1,g2=g2,g3=g3,g4=g4
                              ,g5=g5,g6=g6,g7=g7,g8=g8
                              ,g9=g9,g10=g10,g11=g11,g12=g12
                              ,g13=g13,g14=g14,g15=g15,g16=g16
                              ,g17=g17,g18=g18,g19=g19,g20=g20
                              ,g21=g21,g22=g22,g23=g23,g24=g24
                              ,g25=g25,g26=g26,g27=g27,g28=g28
                              ,g29=g29,g30=g30,g31=g31,g32=g32
                              ,g33=g33,g34=g34,g35=g35,g36=g36,
                              g37=g37,g38=g38)
                       return(res)
                       
                     };
                     self$dimension=5  
                     self$lower=c(704.4148,68.6,0,193,25)
                     self$upper=c(906.3855,288.88,134.75,287.0966,84.1988)
                     self$nConstraints=38
                     self$solu=c(705.17454,  68.60000, 102.90000 ,282.32493,  37.58412); 
                     self$xStart=NA 
                   },
                   callG17=function(){
                     self$fn=function(x){
                       if(x[1]>=0 && x[1]<300){
                         f1<-30*x[1]
                       }else if(x[1]>=300 && x[1]<400){
                         f1<-31*x[1]
                       }
                       
                       
                       if(x[2]>=0 && x[2]<100){
                         f2<-28*x[2]
                       }else if(x[2]>=100 && x[2]<200){
                         f2<-29*x[2]
                       }else if(x[2]>=200 && x[2]<1000){
                         f2<-30*x[2]
                       }
                       
                       obj<-f1+f2
                       
                       h1<--x[1]+300-((x[3]*x[4])/131.078)*cos(1.48477-x[6])+((0.90798*x[3]^2)/131.078)*cos(1.47588)
                       h2<--x[2]    -((x[3]*x[4])/131.078)*cos(1.48477+x[6])+((0.90798*x[4]^2)/131.078)*cos(1.47588)
                       h3<--x[5]    -((x[3]*x[4])/131.078)*sin(1.48477+x[6])+((0.90798*x[4]^2)/131.078)*sin(1.47588)
                       h4<-200      -((x[3]*x[4])/131.078)*sin(1.48477-x[6])+((0.90798*x[3]^2)/131.078)*sin(1.47588)
                       
                       res<-c(objective=obj, equ=h1 , equ=h2,equ=h3,equ=h4)
                       return(res)
                       
                     };
                     self$dimension=6 
                     self$lower=c(0,0,340,340,-1000,0)
                     self$upper=c(400,1000,420,420,1000,0.5236)
                     self$nConstraints=4
                     self$solu=c(201.784467214523659, 99.9999999999999005,
                                 383.071034852773266, 420,-10.9076584514292652, 0.0731482312084287128) 
                     self$xStart=NA
                   },
                   callG18=function(){
                     self$fn=function(x){
                       obj<--0.5*(x[1]*x[4]-x[2]*x[3] + x[3]*x[9]-x[5]*x[9] + x[5]*x[8]- x[6]*x[7])
                       g1<-x[3]^2+x[4]^2-1
                       g2<-x[9]^2-1
                       g3<-x[5]^2+x[6]^2-1
                       g4<-x[1]^2+(x[2]-x[9])^2-1
                       g5<-(x[1]-x[5])^2+(x[2]-x[6])^2-1
                       g6<-(x[1]-x[7])^2+(x[2]-x[8])^2-1
                       g7<-(x[3]-x[5])^2+(x[4]-x[6])^2-1
                       g8<-(x[3]-x[7])^2+(x[4]-x[8])^2-1
                       g9<-x[7]^2+(x[8]-x[9])^2-1
                       g10<-x[2]*x[3]-x[1]*x[4]
                       g11<--x[3]*x[9]
                       g12<-x[5]*x[9]
                       g13<-x[6]*x[7]-x[5]*x[8]
                       res<-c(obj=obj,g1=g1,g2=g2,
                              g3=g3,g4=g4,
                              g5=g5,g6=g6,
                              g7=g7,g8=g8,
                              g9=g9,g10=g10,
                              g11=g11,g12=g12,
                              g13=g13)
                       
                       return(res)
                       
                     };
                     self$dimension=9 
                     self$lower=c(rep(-10,8),0)
                     self$upper=c(rep(10,8),20)
                     self$nConstraints=13
                     self$solu=c(-0.9890005492667746,0.1479118418638228,
                                 -0.6242897641574451,-0.7811841737429015,
                                 -0.9876159387318453,0.1504778305249072,
                                 -0.6225959783340022,-0.782543417629948,
                                 0.0)
                     self$xStart=NA 
                   },
                   callG19=function(){
                     # G19 fitness function
                     aMat19<-matrix(c(
                       -16 ,  2,  0,  1,   0,
                       +0  , -2,  0,0.4,   2,
                       -3.5,  0,  2,  0,   0,
                       +0  , -2,  0, -4,  -1,
                       +0  , -9, -2,  1,-2.8,
                       +2  ,  0, -4,  0,   0,
                       -1  , -1, -1, -1,  -1,
                       -1  , -2, -3, -2,  -1,
                       +1  ,  2,  3,  4,   5,
                       +1  ,  1,  1,  1,   1 ), byrow=T,nrow=10,ncol=5)
                     bVec19<-c(-40,-2,-0.25,-4,-4,-1,-40,-60,5,1)
                     cMat19<-matrix(c(
                       +30, -20, -10, 32, -10,
                       -20,  39,  -6,-31,  32,
                       -10,  -6,  10, -6, -10,
                       +32, -31,  -6, 39, -20,
                       -10,  32, -10,-20,  30), byrow=T,nrow=5,ncol=5)
                     dVec19 <- c(4, 8, 10, 6, 2)
                     eVec19 <- c(-15, -27, -36, -18, -12)
                     
                     
                     fitFunc<-function(x){
                       
                       obj <- -sum(bVec19*x[1:10]) + 2*sum(dVec19*x[11:15]*x[11:15]*x[11:15])
                       for (i in 1:5) obj <- obj + x[10+i]*sum(cMat19[i,]*x[11:15])
                       #obj <- obj 
                       return(obj)
                     }
                     conFunc <- function(x) {
                       
                       res <- rep(0,5)
                       for (j in 1:5) {
                         res[j] <- -2*sum(cMat19[,j]*x[11:15]) -3*dVec19[j]*x[10+j]*x[10+j] -eVec19[j] + sum(aMat19[,j]*x[1:10])
                         #res[j] <- res[j] 
                       }
                       return(res)
                     }
                     self$fn=function(x){
                       res<-c(objective=fitFunc(x), conFunc(x))
                       return(res)
                       
                     };
                     self$dimension=15  
                     self$lower=rep(0,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=5
                     self$solu=c( 
                       0, 0,  3.94600628013917,  0,    3.28318162727873,
                       10, 0, 0, 0,0,
                       0.370762125835098,  0.278454209512692,   0.523838440499861,     0.388621589976956,   0.29815843730292) 
                     self$xStart=NA
                   },
                   callG20=function(){
                     self$fn=function(x){
                       a<-c(0.0693,0.0577,0.05,0.2,
                            0.26,0.55,0.06,0.1,0.12,
                            0.18,0.1,0.09,0.0693,0.0577,
                            0.05,0.2,0.26,0.55,0.06,
                            0.1,0.12,0.18,0.1,0.09)
                       obj<-sum(a*x)
                       #g1,g2,g3
                       e<-c(0.1,0.3,0.4,0.3,0.6,0.3)
                       b<-c(44.094,58.12,58.12,137.4,120.9,170.9,62.501,84.94,
                            133.425,82.507,46.07,60.097,44.094,58.12,58.12,
                            137.4,120.9,170.9,62.501,84.94,133.425,82.507,46.07,60.097)
                       cVec<-c(123.7,31.7,45.7,14.7,84.7,27.7,49.7,7.1,2.1,17.7,0.85,0.64)
                       d<-c(31.244,
                            36.12,34.784,92.7,82.7,91.6,56.708,82.7,80.8,64.517,49.4,49.1)
                       k<-0.7302*530*( 14.7/40)
                       sumX<-sum(x)
                       g123<-c()
                       for(i in c(1:3)){
                         y<-(x[i]+x[i+12])/(sumX+e[i])
                         g123<-c(g123,y) 
                       }
                       
                       #g4,g5,g6
                       g456<-c()
                       for(i in c(4:6)){
                         y<-(x[i+3]+x[i+15])/(sumX+e[i])
                         g456<-c(g456,y) 
                         
                       }
                       #h1,..h12
                       h<-c()
                       for(i in c(1:12)){
                         h1<-x[i+12]/(b[i+12]*sum((x/b)[c(13:24)]))
                         h2<-cVec[i]*x[i]/(40*b[i]*sum((x/b)[c(1:12)]))
                         h<-c(h,(h1-h2))
                       }
                       names(h)<-rep("equ",length(h))
                       h13<-sumX-1
                       h14=sum((x/d)[c(1:12)])+k*sum((x/b)[c(13:24)])-1.671
                       res<-c(objective=obj, g123,g456, h,equ=h13 , equ=h14)
                       return(res)
                       
                     };
                     self$dimension=24   
                     self$lower=rep(0,self$dimension)
                     self$upper=rep(10,self$dimension)
                     self$nConstraints=20
                     self$solu=solu<-c(9.53E-7,
                                       0,
                                       4.21E-3,
                                       1.039E-4,
                                       0,
                                       0,
                                       2.072E-1,
                                       5.979E-1,
                                       1.298E-1,
                                       3.35E-2,
                                       1.711E-2,
                                       8.827E-3,
                                       4.657E-10,
                                       0,
                                       0,
                                       0,
                                       0,
                                       0,
                                       2.868E-4,
                                       1.193E-3,
                                       8.332E-5,
                                       1.239E-4,
                                       2.07E-5,
                                       1.829E-5) 
                     self$xStart=NA
                     self$info="As no feasible solution is known for G20, the provided solution is slightly infeasible, "
                   },
                   callG21=function(){
                     self$fn=function(x){
                       obj<-x[1]
                       g1<--x[1]+35*(x[2]^(0.6))+35*(x[3]^0.6)
                       h1<--300*x[3] + 7500*x[5]- 7500*x[6] - 25*x[4]*x[5] + 25*x[4]*x[6] + x[3]*x[4]
                       h2<-100*x[2] + 155.365*x[4] + 2500*x[7] - x[2]*x[4] - 25*x[4]*x[7] - 15536.5
                       h3<--x[5] + log(-x[4] + 900)
                       h4<--x[6] + log(x[4] + 300)
                       h5<--x[7] + log(-2*x[4] + 700)
                       res<-c(objective=obj,equ=h1 , equ=h2, equ=h3 , equ=h4, equ=h5, g1=g1)
                       return(res)
                       
                     };
                     self$dimension=7  
                     self$lower=c(0,0,0,100,6.3,5.9,4.5)
                     self$upper=c(1000,40,40,300,6.7,6.4,6.25)
                     self$nConstraints=6
                     self$solu=c(193.724510070034967,
                                 5.56944131553368433*(10^-27),
                                 17.3191887294084914,
                                 100.047897801386839,
                                 6.68445185362377892,
                                 5.99168428444264833,
                                 6.21451648886070451); 
                     self$xStart=NA 
                   },
                   callG22=function(){
                     self$fn=function(x){
                       obj<-x[1]
                       g1<--x[1]+x[2]^0.6+x[3]^0.6+x[4]^0.6
                       h1<-x[5] - 100000*x[8] + 10^7
                       h2<-x[6] + 100000*x[8] - 100000*x[9]
                       h3<-x[7] + 100000*x[9] - 5 * 10^7
                       h4<-x[5] + 100000*x[10] - 3.3 * 10^7
                       h5<-x[6] + 100000*x[11] - 4.4 * 10^7
                       h6<-x[7] + 100000*x[12] - 6.6 * 10^7
                       h7<-x[5] - 120*x[2]*x[13]
                       h8<-x[6] - 80*x[3]*x[14]
                       h9<-x[7] - 40*x[4]*x[15]
                       h10<-x[8] - x[11] + x[16]
                       h11<-x[9] - x[12] + x[17]
                       h12<--x[18] + log(x[10] - 100)
                       h13<--x[19] + log(-x[8] + 300)
                       h14<--x[20] + log(x[16])
                       h15<--x[21] + log(-x[9] + 400)
                       h16<--x[22] + log(x[17])
                       h17<--x[8] - x[10] + x[13]*x[18] - x[13]*x[19] + 400
                       h18<-x[8] - x[9] - x[11] + x[14]*x[20] - x[14]*x[21] + 400
                       h19<-x[9] - x[12] - 4.60517*x[15] + x[15]*x[22] + 100
                       res<-c(objective=obj,
                              g1=g1,
                              equ=h1 , equ=h2 , equ=h3 , equ=h4 , equ=h5,
                              equ=h6 , equ=h7 , equ=h8 , equ=h9 , equ=h10,
                              equ=h11, equ=h12, equ=h13, equ=h14, equ=h15,
                              equ=h16, equ=h17, equ=h18, equ=h19)
                       return(res)
                     };
                     self$dimension=22   
                     self$lower=c(rep(0,7),rep(100,2),100.01,rep(100,2),rep(0,3),rep(0.01,2),rep(-4.7,5))
                     self$upper=c(20000,rep(10^6,3),rep(4*10^7,3),299.99,399.99,300,400,600,rep(500,3),300,400,rep(6.25,5))
                     self$nConstraints=20
                     self$solu=c(236.430975504001054,
                                 135.82847151732463,
                                 204.818152544824585,
                                 6446.54654059436416,
                                 3007540.83940215595,
                                 4074188.65771341929,
                                 32918270.5028952882,
                                 130.075408394314167,
                                 170.817294970528621,
                                 299.924591605478554,
                                 399.258113423595205,
                                 330.817294971142758,
                                 184.51831230897065,
                                 248.64670239647424,
                                 127.658546694545862,
                                 269.182627528746707,
                                 160.000016724090955,
                                 5.29788288102680571,
                                 5.13529735903945728,
                                 5.59531526444068827,
                                 5.43444479314453499,
                                 5.07517453535834395)  
                     self$xStart=NA;
                     self$info="Please note that the provided solution is slightly infeasible"
                   },
                   callG23=function(){
                     self$fn=function(x){
                       obj<--9*x[5] - 15*x[8] + 6*x[1] + 16*x[2] + 10*(x[6] + x[7])
                       g1<-x[9]*x[3] + 0.02*x[6] - 0.025*x[5]
                       g2<-x[9]*x[4] + 0.02*x[7] - 0.015*x[8]
                       h1<-x[1] + x[2] - x[3] -x[4]
                       h2<-0.03*x[1] + 0.01*x[2] -x[9]*(x[3] + x[4])
                       h3<-x[3] + x[6] - x[5]
                       h4<-x[4] + x[7] - x[8]
                       
                       res<-c(objective=obj,equ=h1 , equ=h2,equ=h3 , equ=h4,g1=g1,g2=g2)
                       return(res)
                       
                     };
                     self$dimension=9 
                     self$lower=c(rep(0,8),0.01)
                     self$upper=c(300,300,100,200,100,300,100,200,0.03)
                     self$nConstraints=6
                     self$solu=c(0,  100, 0,  100,  0,  0,   100, 200, 0.01) 
                     self$xStart=NA  
                   },
                   callG24=function(){
                   self$fn=function(x){
                     obj<- -x[1] - x[2]
                     g1<- -2*x[1]^4 + 8*x[1]^3 - 8*x[1]^2 + x[2] - 2
                     g2<- -4*x[1]^4 +32*x[1]^3 -88*x[1]^2 + 96*x[1] + x[2] - 36
                     
                     res<-c(objective=obj,g1=g1,g2=g2)
                   };
                   self$dimension=2  
                   self$lower=c(0,0)
                   self$upper=c(3,4)
                   self$nConstraints=2
                   self$solu=c(2.329520197477607, 3.17849307411768)
                   self$xStart=NA 
                 })
               
                 )

# problems<-sprintf("G%02i",c(1:24))
# for (problem in problems){
#   print(sprintf("checking problem %s",problem))
#   if(problem=="G02"||problem=="G03"){
#     newProb<-COP$new(problem,10)
#   }else{
#     newProb<-COP$new(problem)
#   }
# 
#   if(problem=="G13"){
#     testit::assert("lower and upper should have the same length",length(newProb$lower)==length(newProb$upper))
#     testit::assert("lower and solu should have the same length",length(newProb$lower)==length(newProb$solu[1,]))
#     testit::assert("solu and upper should have the same length",length(newProb$solu[1,])==length(newProb$upper))
#     testit::assert("solu should be vector os size dimension",length(newProb$solu[1,])==newProb$dimension)
#     testit::assert("fn should return a vector of size nConstraints+1",length(newProb$fn(newProb$solu[1,]))==newProb$nConstraints+1)
#   }else{
#     testit::assert("lower and upper should have the same length",length(newProb$lower)==length(newProb$upper))
#     testit::assert("lower and solu should have the same length",length(newProb$lower)==length(newProb$solu))
#     testit::assert("solu and upper should have the same length",length(newProb$solu)==length(newProb$upper))
#     testit::assert("solu should be vector os size dimension",length(newProb$solu)==newProb$dimension)
#     testit::assert("fn should return a vector of size nConstraints+1",length(newProb$fn(newProb$solu))==newProb$nConstraints+1)
#   }
#   testit::assert("upper limits must be larger than the lower limits",all((newProb$upper-newProb$lower)>0))
#   testit::assert("solution must be within the lower and upper limit",all((newProb$upper-newProb$solu)>=0))
#   testit::assert("solution must be within the lower and upper limit",all((newProb$solu-newProb$lower)>=0))
# 
# }