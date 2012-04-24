#Modified from 'http://www.r-bloggers.com/r-tools-for-dynamical-systems-r-pplane-to%C2%A0draw%C2%A0phase%C2%A0planes/'

library(deSolve)
require(rethinking)  #this isn't an easily available library, but it isn't neccessary.  It just has the function col.alpha to adjust colors.  You might need to change the color functions if you don't have it.

watersysDES <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = F - d*x - y
    dy = - (d + r)*y
    return(list(c(dx, dy)))
  })
}

watersys <- function(d=0.1, r=0.1, F=0.3){
  function(x,y=NULL){
    if (is.null(y)) {
      y <- x[2]; x <- x[1];
    }
    dx = F - d*x - y
    dy = - (d + r)*y
    return(c(dx, dy))
  }
}

d <- 0.1
r <- 0.1
F <- 0.3

# please source() from http://www.macalester.edu/~kaplan/math135/pplane.r


# modification of phasetraj() in pplane.r

draw.traj <- function(func, Pars, tStart=0, tEnd=1, tCut=10, color = "red") {
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

waterlevel <- function(t,d,r,F,W,T) {
  ((T*(W*d - F))/(exp(-T*(d + r)) - 1))*((exp(-t*(d + r)) - 1)/(1 + d*t)) + ((F*t + W)/(1 + d*t))
}
sales <- function(t,d,r,F,W,T) {
  ((T*(W*d - F)/(exp(-T*(d + r)) - 1)))* (exp(-t*(d + r))/(d+r))
}
plot.traj <- function(d=d,r=r,F=F,W=W,T=T) {
  t.seq <- seq(0,T, length.out=T*10)
  w <- rep(0, length(t.seq))
  s <- w
  for(i in 1:length(t.seq)) {
    w[i] = waterlevel(t=t.seq[i], d=d,r=r,F=F,W=W,T=T)
    s[i] = sales(t=t.seq[i], d=d,r=r,F=F,W=W,T=T)
  }
  points(w[1],s[1], col="red", pch=4, cex=2)
  text(w[1], s[1], bquote(list(omega==.(W), T==.(T))), pos=4)
  lines(w,s,col="red")
  points(w[seq(1,length(w),by=10)],s[seq(1,length(s),by=10)],col=col.alpha("red",0.75),cex=0.75, pch=16)
}

d <- 0.75
r <- 0.05
F <- 10
nullclines(watersys(d=d,r=r, F=F),c(0,20),c(-25,60),250,xlab="Revervoir Level (w) in megaliters", ylab="Sales rate (s) in megaliters/day ", family="serif")
Ws <- c(12, 12, 11, 11)
Ts <- c(10, 20, 10, 20)
for (j in 1:length(Ws)) {
  plot.traj(d=d,r=r,F=F,W=Ws[j],T=Ts[j])
}
phasearrows(watersys(d=d,r=r, F=F),c(0,20),c(-25,60),20)
