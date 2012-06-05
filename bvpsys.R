#Modified from 'http://www.r-bloggers.com/r-tools-for-dynamical-systems-r-pplane-to%C2%A0draw%C2%A0phase%C2%A0planes/'
library(deSolve)
require(bvpSolve)
require(rethinking)  #this isn't an easily available library, but it isn't neccessary.  It just has the function col.alpha to adjust colors.  You might need to change the color functions if you don't have it.

bvpDES <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = y-p
    dy = ry - c/(x^2)
    return(list(c(dx, dy)))
  })
}

bvp <- function(alpha=0.4675, beta=1, a=1.1, b=1.1, r=0.1){
  function(x,y=NULL){
    if (is.null(y)) {
      y <- x[2]; x <- x[1];
    }
    dx = y-p
    dy = ry - c/(x^2)
    return(c(dx, dy))
  }
}



# please source() from http://www.macalester.edu/~kaplan/math135/pplane.r


# modification of phasetraj() in pplane.r

draw.traj <- function(func, Pars, tStart=0, tEnd=1, tCut=10, color = "red") {
  traj <- list()
  print(paste("Click at initial values, press ESC when done"))
  x0 <- locator(1, "p", pch=16, col="red")
  i <- 1
  while(!is.null(x0)){
    out <- as.data.frame(ode(func=func, y=c(x=x0$x, y=x0$y), parms=Pars, times = seq(tStart, tEnd, length.out = tCut), method="bdf"))
    lines(out$x, out$y, col = color)
    points(out$x, out$y, pch=16, cex=0.5, col=col.alpha(acol=color, alpha=0.25))
    traj[[i]] <- out
    x0 <- locator(1, "p", pch=16, col="red")
    i <- i + 1
  }
  return(traj)
}

bvp <- function(a=a, B=B, K=K,n=n){
  function(x,y=NULL){
    if (is.null(y)) {
      y <- x[2]; x <- x[1];
    }
    dx = y-p
    dy = ry - c/(x^2)
    return(c(dx, dy))
  }
}

r <- 1
p <- 1
c <- 0.5

parms <- c(                   # Set our parameters
  r <- 0.1,
  p <- 1,
  c <- 0.5
  )

par(col="#45462f", col.axis="#45462f", col.lab="#45462f", col.main="#45462f", col.sub="#45462f", family="serif", mfrow=c(1,1))
nullclines(bvp(a=a, B=B, K=K,n=n),c(0,2),c(0,2),250,xlab="Fish stock (x)", ylab="Harvest rate (h)", colors=c("#45462f", "#45462f"), xaxs="i", yaxs="i")
phasearrows(bvp(a=a, B=B, K=K,n=n),c(0,2),c(0,2),30, add=TRUE, col="grey60")
points(x0,y0)

traj <- draw.traj(bvpDES, Pars=parms, tStart=0, tEnd=10, tCut=50)




x0 <- 0.2
xT = 0.6
T <- 10
yini <- c(x=x0, y=NA)         # Set initial conditions, leaving unknowns as NA
yend <- c(x=xT, y=NA)         # Set final conidtions, leaving unknowns as NA
times <- seq(0,T, length.out=(T+1)) # Give the the desired time outputs
#out1s <- bvpshoot(yini=yini, x=times, func=bvpDES, yend=yend, parms=parms, guess = 0.1, method="bdf")
out1 <- bvpcol(yini=yini, x=times, func=bvpDES, yend=yend, parms=parms)
points(out1[1,2], out1[1,3], col="springgreen4", pch=16)
points(out1[,2], out1[,3], col=col.alpha("springgreen4", alpha=0.5), pch=16, cex=0.75)
lines(out1[,2], out1[,3], col="springgreen4")

T = 20
times <- seq(0,T, length.out=T+1) # Set teh desired time output
#out2s <- bvpshoot(yini=yini, x=times, func=bvpDES, yend=yend, parms=parms, guess=0.1279, method="bdf")
out2 <- bvpcol(yini=yini, x=times, func=bvpDES, yend=yend, parms=parms)
points(out2[1,2], out2[1,3], col="tomato4", pch=16)
points(out2[,2], out2[,3], col=col.alpha("tomato4", alpha=0.5), pch=16, cex=0.75)
lines(out2[,2], out2[,3], col="tomato4", lwd=2)

T = 30
times <- seq(0,T, length.out=T+1) # Set teh desired time output
out3 <- bvpcol(yini=yini, x=times, func=bvpDES, yend=yend, parms=parms)
points(out3[1,2], out3[1,3], col="slateblue", pch=16)
points(out3[,2], out3[,3], col=col.alpha("slateblue", alpha=0.5), pch=16, cex=0.75)
lines(out3[,2], out3[,3], col="slateblue")

T = 10
yini <- c(x=0.1, y=NA)         # Set initial conditions, leaving unknowns as NA
yend <- c(x=0.3, y=NA)         # Set final conidtions, leaving unknowns as NA
times <- seq(0,T, length.out=900) # Set teh desired time output
out4 <- bvpcol(yini=yini, x=times, func=bvpDES, yend=yend, parms=parms)
par(col="#45462f", col.axis="#45462f", col.lab="#45462f", col.main="#45462f", col.sub="#45462f", family="serif", cex.main=1.3, cex.axis=1.3, cex.lab=1.3)
nullclines(bvp(alpha=alpha, beta=beta, a=a,b=b,r=r),c(0.285,0.315),c(0.229,0.235),250,xlab="Stock (x)", ylab="Harvest rate (h)", colors=c("#45462f", "#45462f"), xaxs="i", yaxs="i")
phasearrows(bvp(alpha=alpha, beta=beta, a=a,b=b,r=r),c(0.285,0.315),c(0.229,0.235),30, add=TRUE, col="grey60")
points(out4[1,2], out4[1,3], col="blue", pch=16)
points(out4[,2], out4[,3], col=col.alpha("blue", alpha=0.5), pch=16, cex=0.75)
lines(out4[,2], out4[,3], col="blue")
points(out4[1,2], out4[1,3], col="blue", pch=16)
points(out4[,2], out4[,3], col=col.alpha("blue", alpha=0.5), pch=16, cex=0.75)
lines(out4[,2], out4[,3], col="blue")