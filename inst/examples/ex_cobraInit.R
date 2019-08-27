
## Initialize cobra. The problem to solve is the sphere function sum(x^2) with the equality   
## constraint that the solution is on a circle with radius 2 and center at c(1,0).
d=2
fName="onCircle"
cobra <- cobraInit(xStart=rep(5,d), fn=function(x){c(obj=sum(x^2),equ=(x[1]-1)^2+(x[2]-0)^2-4)},  
                   fName=fName,lower=rep(-10,d), upper=rep(10,d), feval=40)

## Run cobra optimizer
cobra <- cobraPhaseII(cobra)

## The true solution is at solu = c(-1,0) (the point on the circle closest to the origin)
## where the true optimum is fn(solu)[1] = optim = 1
## The solution found by SACOBRA:
print(getXbest(cobra))
print(getFbest(cobra))

## Plot the resulting error (best-so-far feasible optimizer result - true optimum)
## on a logarithmic scale:
optim = 1
plot(abs(cobra$df$Best-optim),log="y",type="l",main=fName,ylab="error",xlab="iteration")

