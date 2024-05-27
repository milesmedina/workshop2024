# Ecosystem Dynamics and Causality
# CCM with Lotka-Volterra predator-prey model
# Miles Medina, ECCO Scientific, 2024
#

rm(list=ls(all=TRUE))
if(!require(deSolve)){install.packages('deSolve')}; library(deSolve)
if(!require(rEDM)) { install.packages('rEDM') }; library(rEDM)

# Lotka-Volterra ODE
model <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dx <- a*x - b*x*y  # prey
    dy <- c*x*y - d*y   # predator
    list(c(dx,dy))
  })
}
state <- c( x=80, y=20  )  # initial conditions
parameters <- c( a = 0.1, b = 0.002, c = 0.0025, d = 0.2 )
times <- seq(0,500,1)

# Run LV model
out <- ode( y = state, times = times, func = model, parms = parameters, method = 'lsoda' )
x <- out[,2] ;  y <- out[,3]


# Plot time series
par(mfrow=c(2,1))
plot( x, ylim = range(c(0,x,y)),
      type = 'l', lwd = 3, col = rgb(0,0.4,0.8,0.8),
      main = 'Lotka-Volterra solutions', xlab = 'time', ylab = 'population' )
  lines( y, lwd = 3, col = rgb(1,0.1,0.2,0.9) )
  legend( 'bottomleft', bty='n',
          legend=c('prey','predator'), text.font = 2,
          text.col = c(rgb(0,0.4,0.8,0.8),rgb(1,0.1,0.2,0.9)) )
# Plot phase space
plot( y ~ x, type = 'l', main = 'phase space', col = rgb(1,0,0.3,0.7),
      xlab = 'prey', ylab = 'predator' )


# No correlation detected
cor.test( x, y )
lm(y~x) |> summary()


# CCM detects bidirectional causal relationship
library(rEDM)
df1 <- data.frame( time = 1:length(x),
                   prey = x,
                   predator = y )
par(mfrow=c(1,1))
ccm <- CCM( dataFrame = df1,
            E = 2,   # embedding dimension
            tau = -9,   # embedding delay
            exclusionRadius = 11,   # Theiler window
            target = "prey",   # prediction target (cause)
            columns = "predator",   # library (effect) 
            libSizes = "5 485 40",  # string for sequence 'from, to, by'
            sample = 20,   # number of replicate tests at each libSize
            showPlot = TRUE,
            parameterList = TRUE,
            includeData = TRUE
            )
ccm$LibMeans

