% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defaultRI.R
\name{defaultRI}
\alias{defaultRI}
\title{Default settings for \code{\link{repairInfeasRI2}} and \code{\link{repairChootinan}}.}
\usage{
defaultRI(repairMargin = 0.01)
}
\arguments{
\item{repairMargin}{[1e-2] repair only solutions whose infeasibility is less than this margin}
}
\value{
a list with the following elements:
     \item{RIMODE}{  [2] one out of \{0,1,2,3 \} with 0,1: deprecated older versions of RI2, 
         2: the recommended RI2-case, see \code{\link{repairInfeasRI2}}, 
         3: Chootinan's method, see \code{\link{repairChootinan}}  }
     \item{eps1}{  [1e-4] include all constraints not eps1-feasible into the repair mechanism  }
     \item{eps2}{  [1e-4] selects the solution with the shortest shift among all random 
         realizations which are eps2-feasible  }
     \item{q}{  [3.0] draw coefficients \eqn{\alpha_k} from uniform distribution \eqn{U[0,q]}  }
     \item{mmax}{  [1000] draw mmax random realizations  }
     \item{repairMargin}{  repair only solutions whose infeasibility is less than this margin. }
     \item{repairOnlyFresBetter}{  [FALSE] if TRUE, then repair only iterates with \cr
         \code{fitness < so-far-best-fitness + marFres}  }
     \item{marFres}{  [0.0] only relevant if \code{repairOnlyFresBetter==TRUE} }
     
     A solution \eqn{x} is said to be \strong{\eqn{\epsilon}-feasible} for constraint function \eqn{f}, if 
       \deqn{  f(x)+\epsilon \le 0 }
}
\description{
Sets suitable defaults for the repair-infeasible part of SACOBRA. \cr
With the call \code{\link{setOpts}(myRI,defaultRI())} it is possible to extend a partial list 
\code{myRI} to a list containing all \code{ri}-elements (the missing ones are taken from 
\code{defaultRI()})
}
\details{
The \strong{infeasibility} of a solution is its maximum constraint violation 
(0 for a feasible solution).
}
\seealso{
\code{\link{repairInfeasRI2}}, \code{\link{repairChootinan}}
}
\author{
Wolfgang Konen, Cologne University of Applied Sciences
}
