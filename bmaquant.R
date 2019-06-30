# Written: Sept 12, 2003 (Adrian Raftery), 
# Modified: Nov 15, 2007, June 2019 (Hana Sevcikova)

bmaquant <- function (alpha, w, means, vars, niter=14) {
    # Find the alpha quantile of the n-component mixture
    #  using the bisection method.
    # Inputs:
    #  alpha        required quantile
    #  w            vector of weights
    #  means        means of the mixture components
    #  vars         variances of the mixture components
    #  niter        number of iterations in the bisection method
    #
    sig <- sqrt (vars)
    
    # Initialize: Find lower and upper bounds
    lower <- min (means-3*sig)
    upper <- max (means+3*sig)
    Flower <- bmacdf (lower, w, means, vars)
    Fupper <- bmacdf (upper, w, means, vars)
    if (Flower>alpha || Fupper<alpha) return(NA)
    
    # Bisection method
    for (iter in 1:niter) {
        mid <- (lower+upper)/2
        Fmid <- bmacdf (mid,w, means, vars)
        if (Fmid>alpha) upper <- mid
        if (Fmid<alpha) lower <- mid
    }
    mid
}

get.mixture.comp <- function(x, mean, sd, transform=FALSE) {
  p <- dnorm (x, mean=mean, sd=sd)
  if (transform) {
    return (1/(2*x)*p)
  }
  return (p)
}
  
plot.mixtures <- function(w, means, vars, transform=FALSE,
                          simulated.posterior=NULL, mini=NULL, 
                          maxi=NULL, ci=90, trunc.at.0=TRUE, ...) {
  sig <- sqrt(vars)
  x.min <- min (means - 3*sig)
  x.max <- max (means + 3*sig)
  if(is.null(mini)) mini <- x.min
  if(is.null(maxi)) maxi <- x.max
  l <- 1000
  x <- seq(from=x.min, to=x.max, length=l)
  nmix <- length(vars)
  sum.mean <- sum(w * means)
  sum.sd <- sqrt(sum(w^2*vars))
  if (trunc.at.0)
  	y1 <- dnorm.trunc(x, sum.mean, sum.sd, low=0)
  else y1 <- dnorm(x, sum.mean, sum.sd)
  if (transform) {
  	xdenom <- if (trunc.at.0) pmax(x, 0.0001) else x
  	y1 <- 1/(2*xdenom)*y1
  	}
  # Find and plot CIs
  q <- (1-ci/100)/2
  lower <- bmaquant (q, w, means, vars)
  upper <- bmaquant (1-q, w, means, vars)
  plot.x.min <-  max(0, min (means - 3*sig, lower, mini))
  plot.x.max <- max (means + 3*sig, upper, maxi)
  x.coord <- x
  if (transform) {
    x.coord <- xdenom^2
    plot.x.min <- plot.x.min^2
    plot.x.max <- plot.x.max^2
  }
  
  if (!is.null(simulated.posterior)) {
    hist(simulated.posterior, prob=T, xlim=c(plot.x.min,plot.x.max),
         #ylim=c(0, max(y1)), 
         ...)
    lines(density(simulated.posterior), lty=2)
    lines (x.coord, y1, lwd=3, col='red')
  }       
  else {plot (x.coord, y1, type="l", lwd=3, xlim=c(plot.x.min,plot.x.max), ...)}
  for (j in 1:nmix) {
    y1 <- w[j]*get.mixture.comp(x, means[j], sig[j], transform)
    lines (x.coord,y1, col='black')
  }

  m <- bmaquant (.5, w, means, vars)
  if (transform) {
    lower <- lower^2
    upper <- upper^2
    m <- m^2
  }
  abline (v=lower, lty=2)
  abline (v=upper, lty=2)
  abline (v=m, lty=1)
  return (list(xlim=c(plot.x.min, plot.x.max), quant=c(m, lower, upper)))
}

#--------------------------------------
# Here's the R code for finding the mixture CDF at a point
#--------------------------------------

bmacdf <- function (x, w, means, vars) {
    # Sept 12, 2003 (Adrian Raftery)
    # Find the BMA cdf at x
    # Inputs:
    #  x    value at which cdf is required
    
    sig <- sqrt (vars)
    sum (w * pnorm (rep(x,length(means)), means, sig) )
}

dnorm.trunc<-function(x, mean, sd, low = NULL, high  =NULL){
    # normal density truncated between low and high
    plow <- if (!is.null(low)) pnorm(low,mean=mean,sd=sd) else 0
    phigh <- if (!is.null(high)) pnorm(high,mean=mean,sd=sd) else 1
    out<-dnorm(x, mean = mean,sd = sd)/(phigh - plow)
    if (!is.null(low)) out[x<low]<-0
    if (!is.null(high)) out[x>high]<-0
    out
}
