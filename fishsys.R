#Modified from 'http://www.r-bloggers.com/r-tools-for-dynamical-systems-r-pplane-to%C2%A0draw%C2%A0phase%C2%A0planes/'

library(deSolve)
require(rethinking)

fishsysDES <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = r*(x^a) - (r/K)*(x^(1+a)) - q*y*x
    dy = d*y*(p*q*x - (c + O))
    return(list(c(dx, dy)))
  })
}

fishsys <- function(c=2, O=1, p=4, r=1.1, d=0.125, K=20, q=0.25, a=2){
  function(x,y=NULL){
    if (is.null(y)) {
      y <- x[2]; x <- x[1];
    }
    dx = r*(x^a) - (r/K)*(x^(1+a)) - q*y*x
    dy = d*y*(p*q*x - (c + O))
    return(c(dx, dy))
  }
}

c <- 2
O <- 1
p <- 4
r <- 1.1
d <- 0.125
K <- 20
q <- 0.25
a <- 2

# please source() from http://www.macalester.edu/~kaplan/math135/pplane.r
nullclines(fishsys(),c(0,25),c(0,50),40)
phasearrows(fishsys(),c(0,25),c(0,50),20, add=TRUE)

# modification of phasetraj() in pplane.r

draw.traj <- function(func, Pars, tStart=0, tEnd=1, tCut=10, loc.num=1, color = "red") {
  traj <- list()
  print(paste("Click at initial values, press ESC when done"))
  x0 <- locator(1, "p", pch=16, col="red")
  i <- 1
  while(!is.null(x0)){
    out <- as.data.frame(ode(func=func, y=c(x=x0$x, y=x0$y), parms=Pars, times = seq(tStart, tEnd, length = tCut)))
    lines(out$x, out$y, col = color)
    points(out$x, out$y, pch=16, cex=0.5, col=col.alpha(acol=color, alpha=0.25))
    traj[[i]] <- out
    x0 <- locator(1, "p", pch=16, col="red")
    i <- i + 1
  }
  return(traj)
}


nullclines(fishsys(c=2, O=1, p=4, r=1.1, d=0.0005, K=20, q=0.25, a=2),c(0.001,25),c(0.001,30),250)
#showcontours(fishsys(c=2, O=1, p=4, r=1.1, d=0.0005, K=20, q=0.25, a=2),c(0.01,25),c(0.01,30),250)
phasearrows(fishsys(c=2, O=1, p=4, r=1.1, d=0.0005, K=20, q=0.25, a=2),c(0,25),c(0,30),20, add=TRUE)
traj <- draw.traj(func=fishsysDES, Pars=c(c=2, O=1, p=4, r=1.1, d=0.0005, K=20, q=0.25, a=2), tEnd=500, tCut=10000, loc.num=4)
phasetraj(fishsys(c=2, O=1, p=4, r=1.1, d=0.0005, K=20, q=0.25, a=2), tend=10)
