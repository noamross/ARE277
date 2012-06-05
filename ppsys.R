library(deSolve)
require(bvpSolve)
ppDES <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = r*x*(1 - x/K) - a*x*y - q1*E1*x
    dy = (B + a*x)*y - (n/2)*(y^2) - q2*E2*y
    dE1 = (1/(2*c1))*(p1*q1*x*(r*(1 - 2*x/K) - a*y - d) + (c0 + 2*c1*E1)*(r*x/K + d) + q1*a*x*(p2*y - (c0 + 2*c1*E2)/q2))
    dE2 = (1/(2*c1))*(p2*q2*y*(B + a*x - n*y - d) + (c0 + 2*c1*E2)*(0.5*n*y + d) - q2*a*y*(p1*x - (c0 + 2*c1*E1)/q1))
    return(list(c(dx, dy, dE1, dE2)))
  })
}

r <- 0.85
n <- 1
K <- 1
a <- 0.5
B <- 0.25
q1 <- 1
q2 <- 1
d <- 0.05
p1 <- 2
p2 <- 2
c0 <- 0.2
c1 <- 0.8
T <- 25
x0 <- (n*r - 2*a*B)/(2*a^2 + (n*r/K))
y0 <- (2*r*(a*K+B))/(n*r + 2*(a^2)*K)
x0; y0
E1T <- 0
E2T <- 0

parms <- c(r=r, n=n, K=K, a=a, B=B, q1=q1, q2=q2, d=d, p1=p1, p2=p2, c0=c0, c1=c1)

yini <- c(x=x0, y=y0, E1=NA, E2=NA)         
yend <- c(x=0, y=0, E1=NA, E2=NA)         
times <- seq(from=0, to=T, length.out=T+1)
out1 <- bvpshoot(yini=yini, x=times, func=ppDES, yend=yend, parms=parms, guess=c(0,0))
plot(out1, xlab="Time", ylab="Value")

a <- 0; parms <- c(r=r, n=n, K=K, a=a, B=B, q1=q1, q2=q2, d=d, p1=p1, p2=p2, c0=c0, c1=c1)
x0 <- (n*r - 2*a*B)/(2*(a^2) + (n*r/K))
y0 <- 2*r*(a*K+B)/(n*r + 2*(a^2)*K)
yini <- c(x=x0, y=y0, E1=NA, E2=NA)         
yend <- c(x=0, y=0, E1=NA, E2=NA)
out2 <- bvpshoot(yini=yini, x=times, func=ppDES, yend=yend, parms=parms, guess=c(0,0))
plot(out2, xlab="Time", ylab="Value")

ppDES <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = r*x*(1 - x/K) - a*x*y - q1*E1*x
    dy = (B + a*x)*y - (n/2)*(y^2) - q2*E2*y
    dE1 = (1/(2*c1))*(p1*q1*x*(r*(1 - 2*x/K) - a*y - d) + (c0 + 2*c1*E1)*(r*x/K + d) + q1*a*x*(p2*y - (c0 + 2*c1*E2)/q2))
    dE2 = 0
    return(list(c(dx, dy, dE1, dE2)))
  })
}
a <- 0.5; parms <- c(r=r, n=n, K=K, a=a, B=B, q1=q1, q2=q2, d=d, p1=p1, p2=p2, c0=c0, c1=c1)
out3 <- bvpshoot(yini=yini, x=times, func=ppDES, yend=yend, parms=parms, guess=c(0,0))
plot(out3)



bvp <- function(a=a, B=B, K=K,n=n){
  function(x,y=NULL){
    if (is.null(y)) {
      y <- x[2]; x <- x[1];
    }
    dx = r*x*(1-x/K) - a*x*y
    dy = (B + a*x)*y - (n/2)*(y^2)
    return(c(dx, dy))
  }
}

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

par(col="#45462f", col.axis="#45462f", col.lab="#45462f", col.main="#45462f", col.sub="#45462f", family="serif")
nullclines(bvp,c(0,1),c(0,0.4),250,xlab="Fish stock (x)", ylab="Harvest rate (h)", colors=c("#45462f", "#45462f"), xaxs="i", yaxs="i")
phasearrows(bvp,c(0,1),c(0,0.4),30, add=TRUE, col="grey60")
