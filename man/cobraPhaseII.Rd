% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cobraPhaseII.R
\name{cobraPhaseII}
\alias{cobraPhaseII}
\title{Improve the feasible solution by searching new infill points}
\usage{
cobraPhaseII(cobra)
}
\arguments{
\item{cobra}{an object of class COBRA, this is a (long) list containing all settings
from \code{\link{cobraInit}}}
}
\value{
\code{cobra}, an object of class COBRA from \code{\link{cobraInit}}, 
   enhanced here by the following elements (among others):
     \item{\code{fn}}{ function accepting a \code{d}-dimensional vector \eqn{\vec{x}} and
           returning an \code{(1+m+r)}-vector \code{c(}\eqn{f,g_1,\ldots,g_m,h_1,\ldots,h_r}\code{)}. This
           function may be a rescaled and plog-transformed version of the original \code{fn} 
           passed into \code{\link{cobraInit}}. The original \code{fn} is stored in 
           \code{cobra$originalFn}. }
     \item{\code{df}}{  data frame with summary of the optimization run (see below)}
     \item{\code{df2}}{  data frame with additional summary information (see below)}
     \item{\code{dftr}}{  data frame with additional summary information for TR (see below)}
     \item{\code{A}}{ \code{(feval x d)}-matrix containing all evaluated points 
           in input space. If rescale==TRUE, all points are in \strong{rescaled} input space. }
     \item{\code{Fres}}{ a vector of the objective values of all evaluated points }
     \item{\code{Gres}}{ a \code{(feval x m)}-matrix of the constraint values of all evaluated points }
     \item{\code{predC}}{ a \code{(feval x m)}-matrix with the prediction of  
           \code{cobra$constraintSurrogates} at all evaluated points  }
     \item{\code{fbest}}{ the best feasible objective value found }
     \item{\code{xbest}}{ the point in input space yielding the best feasible objective value }
     \item{\code{ibest}}{ the corresponding iteration number (row of cobra$df, of cobra$A)}
     \item{\code{PLOG}}{ If TRUE, then the objective surrogate model is trained on the 
           \code{\link{plog}}-transformed objective function. }
  Note that \code{cobra$Fres}, \code{cobra$fbest}, \code{cobra$fbestArray} and similar contain 
  always the objective values of the orignial function \code{cobra$fn[1]}. (The surrogate models 
  may be trained on a \code{\link{plog}}-transformed version of this function.)
  
  \code{feval} = \code{cobra$feval} is the maximum number of function evaluations.\cr 
  
  The data frame \code{cobra$df} contains one row per iteration with columns 
  \describe{
     \item{iter}{ iteration index }
     \item{y}{   true objective value \code{Fres} }
     \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
           plog-transformed training data, but \code{predY} is transformed back to the original 
           objective range. NA for the initial design points.}
     \item{predSolu}{  surrogate objective value at best-known solution \code{cobra$solu}, if given. 
           If \code{cobra$solu} is NULL, take the current point instead. Note: The surrogate may be trained on  
           plog-transformed training data, but \code{predSolu} is transformed back to the original 
           objective range. NA for the initial design points.}
     \item{feasible}{ boolean indicating the feasibiltiy of infill point }
     \item{feasPred}{ boolean indicating if each infill point is feasible for \code{cobra$constraintSurrogates} }
     \item{nViolations}{ number of violated constraints }
     \item{maxViolation}{ maximum constraint violation. }
     \item{FEval}{  number of function evaluations in sequential optimizer. NA if it was a repair step }
     \item{Best}{  ever-best feasible objective value \code{fbest}. As long as there is 
           no feasible point, take among those with minimum number of violated constraints the
           one with minimum Fres. }
     \item{optimizer}{ e.g. "COBYLA"  }
     \item{optimizationTime}{  in sec}
     \item{conv}{ optimizer convergence code }
     \item{dist}{ distance of the current point (row of \code{cobra$A}) to the true solution 
           \code{cobra$solu} in rescaled space. If there is more than one solution, take the one
           which has the minimum distance element (since this is the solution to which the 
           current run converges). }
     \item{distOrig}{ same as \code{dist}, but in original space  }
     \item{XI}{  the DRC element used in the current iteration }
     \item{seed}{ the used seed in every run }
  }

  The data frame \code{cobra$df2} contains one row per phase-II-iteration with columns 
  \describe{
     \item{iter}{ iteration index}
     \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
           plog-transformed training data, but \code{predY} is transformed back to the original 
           objective range. NA for the initial design points.}
     \item{predVal}{   surrogate objective value + penalty }
     \item{predSolu}{   surrogate objective value at true solution (see \code{cobra$df$predSolu}) }
     \item{predSoluPenal}{   surrogate objective value + penalty at true solution (only diagnostics)}
     \item{sigmaD}{ the sigmaD element used in the current iteration (see \code{\link{cobraInit}})  }
     \item{penaF}{ penalty factor used in the current iteration (see \code{\link{cobraInit}}) }
     \item{XI}{  the DRC element used in the current iteration }
     \item{EPS}{ the current used margin for constraint function modeling (see \code{epsilonInit} in \code{\link{cobraInit}} ) }
  }
}
\description{
Improve the feasible solution using the SACOBRA optimizer phase II
by searching new infill points with the help of RBF surrogate models. 
May be even called if no feasible solution is found yet, then phase II will try to find
feasible solutions. \cr
The problem to solve iteratively is: \cr
\deqn{ \mbox{Minimize}\quad  f(\vec{x}) , \vec{x} \in [\vec{a},\vec{b}] \subset \mathbf{R}^d }
\deqn{ \mbox{subject to}\quad g_i(\vec{x}) \le 0, i=1,\ldots,m    }
\deqn{ \mbox{~~~~~~~~~~}\quad\quad h_j(\vec{x}) = 0, j=1,\ldots,r.    } \cr
In this phase the main optimization steps are repeated in a loop as long as the budget is not exhausted.
In every iteration the surrogate models are updated and an optimization on the surrogates is done in order 
to find a better feasible solution.
}
\examples{
## Initialize cobra. The problem to solve is the unconstrained sphere function sum(x^2).   
 
## In version 1.1 and higher there is no need for defining a dummy 
## constraint function for the unconstrained problems
d=2
fName="sphere"
cobra <- cobraInit(xStart=rep(5,d), fName=fName,
                   fn=function(x){c(obj=sum(x^2))},  
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

}
\seealso{
\code{\link{cobraPhaseI}}, \code{\link{cobraInit}}
}
\author{
Wolfgang Konen, Samineh Bagheri, Patrick Koch, Cologne University of Applied Sciences
}
