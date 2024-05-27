# Ecosystem Dynamics and Causality
# Signals and Fourier transformation
# Miles Medina, ECCO Scientific, 2024



# cartoon: time series = signal + noise
rm(list=ls(all=TRUE)) 
  # generate data
    t <- 0:200
    x <- sin( 12*t ) * exp(-0.001*t)
    x <- x + 0.1*rev(x) + x * 1.4*cos( 12*t + 1.5 )
    x <- x - mean(x)
    xn <- x + 0.2*rnorm(length(t))
    xn <- xn - mean(xn)
  
  # plots
    par(mfrow=c(3,1), mar=c(5,5,2,1))
    ylims <- range(c(x,xn))
    
    # time series
    plot( x = t, y = xn, type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'Time series', xlab = "", ylab = "", las = 1, bty = "L",
          cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, ylim = ylims )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
    
    # signal
    plot( x = t, y = xn, type = 'l', lwd = 4, col = rgb(0,0,0,0.2),
          main = 'Signal', xlab = "", ylab = "", las = 1, bty = "L",
          cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, ylim = ylims )
    lines( x = t, y = x, type = 'l', lwd = 3, col = rgb(1,0,0.2,0.6) )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )
    
    # noise
    plot( x = t, y = (x-xn), type = 'l', lwd = 3, col = rgb(0,0.4,0.8,0.5),
          main = 'Noise', xlab = "time", ylab = "", las = 1, bty = "L",
          cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, ylim = ylims )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1) )
    abline( h = axTicks(2), col = rgb(0,0,0,0.1) )



# example 1: Fourier transform with a clean sine wave
rm(list=ls(all=TRUE)) 
  # generate data
    t <- 0:200
    x <- sin( 2*pi*t/12 )
    
  # Fourier transform
    spec <- spectrum( x, method = 'ar', plot = FALSE )
    df <- data.frame( power = spec$spec, period = 1/spec$freq )
    df <- df[ order( df$period ), ]
  
  # plots
    par(mfrow=c(2,1), mar=c(5,4,2,1))
    
    # signal
    plot( x = t, y = x, type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'Time series', xlab = "t", ylab = "x(t)", las = 1, bty = "L",
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    
    # periodogram
    plot( power ~ period, data = df[1:480,],
          type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'Periodogram', las = 1, bty = "L",
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3 )
    points( power ~ period, data = df, pch = 16, cex = 0.8 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    
  # power spectrum
    df[ order(df$power, decreasing = TRUE ), ] |> head()

    
    
# example 2: Fourier transform with a more complex signal
rm(list=ls(all=TRUE)) 
  # generate data
    t <- 1:200
    x_1 <- sin( 2*pi*t/24 )
    x_2 <- sin( 2*pi*t/12 ) * 1.5
    x_3 <- sin( 2*pi*t/6 - pi/8 )
    x <- x_1 + x_2 + x_3
    
  # Fourier transform
    spec <- spectrum( x, method = 'ar', plot = FALSE )
    df <- data.frame( power = spec$spec, period = 1/spec$freq )
    df <- df[ order( df$period ), ]
    
  # ts and pgram plots
    par(mfrow=c(2,1), mar=c(5,4,2,1))
    
    # signal
    plot( x = t, y = x, type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'Time series', xlab = "t", ylab = "x(t)", las = 1, bty = "L",
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    
    # periodogram
    plot( power ~ period, data = df[1:490,],
          type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'Periodogram', las = 1, bty = "L",
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3 )
    points( power ~ period, data = df, pch = 16, cex = 0.8 )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    
  # power spectrum
    df[ order(df$power, decreasing = TRUE ), ] |> head(12)
    
  # signal component plots
    par(mfrow=c(3,1), mar=c(5,4,2,1))
    ylims <- range( c(x_1,x_2,x_3) )
    
    # x_1
    plot( x = t, y = x_1, type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'x_1 series', xlab = "t", ylab = "x_1(t)", las = 1, bty = "L",
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, ylim = ylims )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    
    # x_2
    plot( x = t, y = x_2, type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'x_2 series', xlab = "t", ylab = "x_2(t)", las = 1, bty = "L",
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, ylim = ylims )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    
    # x_3
    plot( x = t, y = x_3, type = 'l', lwd = 3, col = rgb(0,0,0,0.5),
          main = 'x_3 series', xlab = "t", ylab = "x_3(t)", las = 1, bty = "L",
          cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, ylim = ylims )
    abline( v = axTicks(1), col = rgb(0,0,0,0.1))
    abline( h = axTicks(2), col = rgb(0,0,0,0.1))
    