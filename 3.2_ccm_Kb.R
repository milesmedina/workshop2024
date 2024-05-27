# Ecosystem Dynamics and Causality
# CCM for Karenia brevis in southwest Florida (Medina et al., 2022)
# Miles Medina, ECCO Scientific, 2024
#

rm(list=ls(all=TRUE)) 

# Load libraries
if(!require(rEDM)) { install.packages('rEDM') }; library(rEDM)


# Load data
dat <- read.csv("./dat/CCM_Kb.csv")
colnames( dat )

# Select variables for CCM test
y <- "Kbrevis"  # effect
x <- "S79.TN"  # cause
df1 <- dat[,c("Date",x,y)] |> na.omit()
df1$Date <- df1$Date |> as.Date(format="%m/%d/%Y") # format dates
df1[,c(2,3)] <- apply( df1[,c(2,3)], 2, scale )  # scale signals to mean=0, sd=1
dim( df1 )

# Run CCM and plot results
ccm <- CCM( dataFrame = df1,
            E = 3,   # embedding dimension
            tau = -6,   # embedding delay
            exclusionRadius = 19,   # Theiler window
            target = x,   # prediction target (cause)
            columns = y,   # library (effect) 
            libSizes = "15 235 20",  # string for sequence 'from, to, by'
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
