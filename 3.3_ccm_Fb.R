# Ecosystem Dynamics and Causality
# CCM for Jon's FL Bay data
# Miles Medina, ECCO Scientific, 2024
#

rm(list=ls(all=TRUE)) 

# Load libraries
if(!require(rEDM)) { install.packages('rEDM') }; library(rEDM)


# Load data
dat <- read.csv("./dat/ssa_FB_ac.csv")
colnames( dat )

# Time series plots
par(mfrow=c(4,2))
for( i in 2:ncol(dat) ){
  plot( dat[,i], type = 'l', xlab='time',
        main=colnames(dat)[i] )
}


#   Phase space reconstruction
##

colnames( dat )
x <- dat[,8]

# Find embedding parameters
  par(mfrow=c(3,1))
  # Set embedding delay (d) using AMI function
  ami <- mutual( x, lag.max = 20 )  # average mutual information function
  local.d.min <- rollapply( ami, 3, function (x) which.min(x)==2 )  # local minima
  d <- as.numeric( which(local.d.min==TRUE)[1] )  # first local min
  points( x = d, y = ami[d+1], cex = 2, col = rgb(1,0,0,0.7) )
  d

  # Set Theiler window (tw) using ACF
  ac <- acf( x, lag.max = 20 )  # autocorrelation function
  ac$acf_abs <- ac$acf |> abs()  # absolute values of acf
  local.tw.min <- rollapply(ac$acf_abs, 3, function (x) which.min(x)==2 )  # abs local minima
  tw <- as.numeric( which(local.tw.min==TRUE)[1] )  # first abs local min
  points( x = tw, y = ac$acf[tw+1], cex = 2, col = rgb(1,0,0,0.7) )
  tw

  # Set embedding dimension (m) using false nearest neighbors test
  fnn <- false.nearest( series=x, m=6, d=d, t=tw, eps=sd(x), rt=10 )
  threshold <- 0.15
  plot( fnn[1,], type = 'b', main = "FNN", pch = 16, cex = 2,
        xlab = 'dim', ylab = 'proportion of false neighbors' )
    abline( h = threshold, lty = 2, col = rgb(0,0,0,0.6) )
  m <- as.numeric( which( fnn[1,] <= threshold )[1] )
  points( x = m, y = fnn[1,m], cex = 3, col = rgb(1,0,0,0.7) )
  
  # space-time separation plot
  par(mfrow=c(1,1))
  stp <- stplot( series = x, m = m, d = d, mdt = length(x) )
  
# Time-delay embedding
  Mx <- embedd( x, m = m, d = d ) |> as.data.frame()

# Plot phase space reconstruction
  # Rename columns
  for(i in 1:m){
    if(i==1){ names(Mx)[i]<-'x(t)'
    } else {
      names(Mx)[i] <- paste0('x(t+',d*(i-1),')')
    }
  }  # // end i
  # Plotly
  plot_ly( Mx, x = ~Mx[,1], y = ~Mx[,2], z = ~Mx[,3],
           type = 'scatter3d', mode = 'lines',
           opacity = 0.75, line = list(width = 6, reverscale = FALSE) ) |> 
  layout( title = 'Reconstructed attractor',
          scene = list( xaxis = list(title=names(Mx)[1]),
                        yaxis = list(title=names(Mx)[2]),
                        zaxis = list(title=names(Mx)[3])
          ) )
  



#   CCM
##

# Select variables for CCM test
y <- "gbchl"  # effect
x <- "rbphos"  # cause
df1 <- dat
df1[,c(2:ncol(df1))] <- apply( df1[,c(2:ncol(df1))], 2, scale )  # scale signals to mean=0, sd=1
dim( df1 )


# Run CCM and plot results
par(mfrow=c(1,1))
ccm <- CCM( dataFrame = df1,
            E = 3,   # embedding dimension
            tau = -3,   # embedding delay
            exclusionRadius = 7,   # Theiler window
            target = x,   # prediction target (cause)
            columns = y,   # library (effect) 
            libSizes = "6 96 5",  # string for sequence 'from, to, by'
            sample = 100,   # number of replicate tests at each libSize
            showPlot = TRUE,
            parameterList = TRUE,
            includeData = TRUE
)

# Output: CCM summary table
ccm$LibMeans

# Output: Results for each 'y xmap x' test
ccm$CCM1_PredictStat |> tail()


# A nicer plot for 'y xmap x' tests
plot( x = ccm$LibMeans$LibSize,
      y = ccm$LibMeans[,2],
      main = paste( y, 'xmap', x),
      ylab = "Prediction skill", xlab = "Library size",
      ylim = range( 0, range(ccm$LibMeans[,2]), 1 ),
      type = 'l', col = 1, lwd = 1 )
# grid lines
abline( h = axTicks(2), col = rgb(0,0,0,0.2) )
abline( v = axTicks(1), col = rgb(0,0,0,0.2) )
abline( h = 0 )
# results of individual tests
points( x = ccm$CCM1_PredictStat$LibSize,
        y = ccm$CCM1_PredictStat$rho, 
        pch = 16, col = rgb(1,0,0,0.1)
)
# Redraw mean prediction skill curve
lines( x = ccm$LibMeans$LibSize,
       y = ccm$LibMeans[,2],
       lwd = 3, col = 1 )
