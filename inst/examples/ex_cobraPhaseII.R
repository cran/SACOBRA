
## Initialize cobra. The problem to solve is the unconstrained sphere function sum(x^2).   
d=2
fName="sphere"
cobra <- cobraInit(xStart=rep(5,d), fn=function(x){c(obj=sum(x^2))}, fName=fName, 
                   lower=rep(-10,d), upper=rep(10,d), feval=40)

## Run cobra optimizer
cobra <- cobraPhaseII(cobra)

## The true solution is at solu = c(0,0)
## where the true optimum is fn(solu)[1] = optim = 0
## The solution found by SACOBRA:
print(getXbest(cobra))
print(getFbest(cobra))

## Plot the resulting error (best-so-far feasible optimizer result - true optimum)
## on a logarithmic scale:
optim = 0
plot(cobra$df$Best-optim,log="y",type="l",ylab="error",xlab="iteration",main=fName)

