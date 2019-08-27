## solve G11 problem nrun times  and plot the results of all nrun runs
nrun=4
feval=30

## Define the constraint problem: G11
fn <- function(x) {
  y<-x[1]*x[1]+((x[2]-1)^2)
  y<-as.numeric(y)

  g1 <- as.numeric(+(x[2] - x[1]^2))
  
  return(c(objective=y, g1=g1))
}
funcName="G11"
lower<-c(-1,-1) 
upper<-c(+1,+1) 

## Initializing and running cobra
cobra <- cobraInit(xStart=c(0,0), fn=fn, fName=funcName, lower=lower, upper=upper,
                   feval=feval, initDesPoints=3*2, DOSAC=1,cobraSeed=1)

mres <- multiCOBRA(fn,lower,upper,nrun=nrun,feval=feval,optim=0.75
                  ,cobra=cobra,funcName=funcName
                  ,ylim=c(1e-12,1e-0),plotPDF=FALSE,startSeed=42)
  

## There are two true solutions at 
## solu1 = c(-sqrt(0.5),0.5) and solu2 = c(+sqrt(0.5),0.5)
## where the true optimum is f(solu1) = f(solu2) = -0.75
## The solution from SACOBRA is close to one of those solutions:
print(getXbest(mres$cobra))
print(getFbest(mres$cobra))

print(mres$z2)
