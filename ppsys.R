library(deSolve)
require(bvpSolve)
require(rethinking)
ppDES <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = r*x*(1 - x/K) - a*x*y - q1*E1*x
    dy = (B + a*x)*y - (n/2)*(y^2) - q2*E2*y
    dE1 = (1/(2*c1))*(p1*q1*x*(r*(1 - 2*x/K) - a*y - d) + (c0 + 2*c1*E1)*(r*x/K + d) - q1*a*x*(p2*y - (c0 + 2*c1*E2)/q2))
    dE2 = (1/(2*c1))*(p2*q2*y*(B + a*x - n*y - d) + (c0 + 2*c1*E2)*(0.5*n*y + d) + q2*a*y*(p1*x - (c0 + 2*c1*E1)/q1))
    return(list(c(dx, dy, dE1, dE2)))
  })
}

r <- 0.85
n <- 1
K <- 1
a <- 0
B <- 0.25
q1 <- 1
q2 <- 1
d <- 0.05
p1 <- 2
p2 <- 2
c0 <- 0.2
c1 <- 0.8
T <- 25
x0 <- (n*r - 2*a*B)/(2*(a^2) + (n*r/K))
y0 <- 2*a*r*(K+B)/(n*r + 2*(a^2)*K)
x0; y0
E1T <- 0
E2T <- 0

parms <- c(r=r, n=n, K=K, a=a, B=B, q1=q1, q2=q2, d=d, p1=p1, p2=p2, c0=c0, c1=c1)

yini <- c(x=1, y=1, E1=NA, E2=NA)         # Set initial conditions, leaving unknowns as NA
yend <- c(x=NA, y=NA, E1=0, E2=0)         # Set final conidtions, leaving unknowns as NA
times <- seq(from=0, to=T, length.out=T+1)
out1 <- bvpcol(yini=yini, x=times, func=ppDES, yend=yend, parms=parms)

par(col="#45462f", col.axis="#45462f", col.lab="#45462f", col.main="#45462f", col.sub="#45462f", family="serif", mfrow=c(2,1))

points(out2[1,2], out2[1,3], col="tomato4", pch=16)
points(out2[,2], out2[,3], col=col.alpha("tomato4", alpha=0.5), pch=16, cex=0.75)
lines(out2[,2], out2[,3], col="tomato4", lwd=2)

some <- ode(y=c(x=x0,y=y0,E1=1,E2=1), times, ppDES, parms, method="bdf"); plot(some);tail(some)

