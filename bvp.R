## Solving Boundary Point Value Problems in R
require(bvpSolve)   # This calls the package

fishsysDES <- function (Time, State, Pars) {   #Define the system of ODEs
  with(as.list(c(State, Pars)), {
    dx = r*(x^a) - (r/K)*(x^(1+a)) - q*y*x
    dy = d*y*(p*q*x - (c + O))
    return(list(c(dx, dy)))
  })
}

yini <- c(x=NA, y=20)         # Set initial conditions, leaving unknowns as NA
yend <- c(x=NA, y=20)         # Set final conidtions, leaving unknowns as NA
times <- seq(0,20, by=0.05)   # Set teh desired time output
parms <- c(                   # Set our parametsrs
  c=2,
  O=1,
  p=4,
  r=1.1,
  d=0.125,
  K=20,
  q=0.25,
  a=2
)
  
time <- system.time(     #Do this inside system.time to clock the program
  out <- bvpshoot(yini=yini, x=times, func=fishsysDES, yend=yend, parms=parms) 
)
plot(out)