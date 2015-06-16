## Function for plotting priors, likelihoods, and posteriors for binomial data
## Output consists of a plot and various statistics
## PS and PF determine the shape of the prior distribution.
## PS = prior success, PF = prior failure for beta dist. 
## PS = 1, PF = 1 corresponds to uniform(0,1) and is default. If left at default, posterior will be equivalent to likelihood
## k = number of observed successes in the data, n = total trials. If left at 0 only plots the prior dist.
## null = Is there a point-null hypothesis? null = NULL leaves it out of plots and calcs
## CI = Is there a relevant X% credibility interval? .95 is recommended and standard

plot.beta <- function(PS = 1, PF = 1, k = 0, n = 0, null = NULL, CI = NULL) {
      
                x = seq(.001, .999, .001) ##Set up for creating the distributions
                y1 = dbeta(x, PS, PF) # data for prior curve
                y3 = dbeta(x, PS + k, PF + n - k) # data for posterior curve
                y2 = dbeta(x, 1 + k, 1 + n - k) # data for likelihood curve, plotted as the posterior from a beta(1,1)
                
                
        plot(x, y1, xlim=c(0,1), ylim=c(0, 1.25 * max(y1,y3,1.6)), type = "l", ylab= "Density", lty = 2,
             xlab= "Probability of success", las=1, main="Prior-to-Posterior Transformation with Binomial Data",lwd=3,
             cex.lab=1.5, cex.main=1.5, col = "skyblue", axes=FALSE)
        
        axis(1, at = seq(0,1,.2)) #adds custom x axis
        axis(2, las=1) # custom y axis
        
        
        
        if(n != 0){
                #if there is new data, plot likelihood and posterior
                lines(x, y2, type = "l", col = "darkorange", lwd = 2, lty = 3)
                lines(x, y3, type = "l", col = "darkorchid1", lwd = 5)
                #legend("topright", c("Prior", "Posterior", "Likelihood"), col = c("skyblue", "darkorchid1", "darkorange"),
                      # lty = c(2,1,3), lwd = c(3,5,2), bty = "n", inset = c(-.3, 0), y.intersp = .55,
                      # x.intersp = .40)
                legend("topleft", c("Prior", "Posterior", "Likelihood"), col = c("skyblue", "darkorchid1", "darkorange"), 
                       lty = c(2,1,3), lwd = c(3,5,2), bty = "n", y.intersp = .55, x.intersp = .1, seg.len=.7)
        }
        
        ##Specified CI% but no null? Calc and report only CI
        if(is.numeric(CI) == T && is.numeric(null) == F){
                CI.low <- qbeta((1-CI)/2, PS + k, PF + n - k)
                CI.high <- qbeta(1-(1-CI)/2, PS + k, PF + n - k)
                return( list( "Posterior CI lower" = CI.low, "Posterior CI upper" = CI.high))
        }
        
        ##Specified null but not CI%? Calculate and report BF only 
        if(is.numeric(null) == T && is.numeric(CI) == F){
                null.H0 <- dbeta(null, PS, PF)
                null.H1 <- dbeta(null, PS + k, PF + n - k)
                CI.low <- qbeta((1-CI)/2, PS + k, PF + n - k)
                CI.high <- qbeta(1-(1-CI)/2, PS + k, PF + n - k)
                return( list("BF01 (in favor of H0)" = null.H1/null.H0, "BF10 (in favor of H1)" = null.H0/null.H1
                             ))
        }
        
        ##Specified both null and CI%? Calculate and report both
        if(is.numeric(null) == T && is.numeric(CI) == T){
                null.H0 <- dbeta(null, PS, PF)
                null.H1 <- dbeta(null, PS + k, PF + n - k)
                CI.low <- qbeta((1-CI)/2, PS + k, PF + n - k)
                CI.high <- qbeta(1-(1-CI)/2, PS + k, PF + n - k)
                return( list("BF01 (in favor of H0)" = null.H1/null.H0, "BF10 (in favor of H1)" = null.H0/null.H1,
                             "Posterior CI lower" = CI.low, "Posterior CI upper" = CI.high))
        }
        
        ## adds null points on prior and posterior curve if null is specified and there is new data
        if(is.numeric(null) == T && n != 0){
                ## Adds points on the distributions at the null value if there is one and if there is new data
                points(null, dbeta(null, PS, PF), pch = 21, bg = "blue", cex = 1.5)
                points(null, dbeta(null, PS + k, PF + n - k), pch = 21, bg = "darkorchid", cex = 1.5)
                abline(v=null, lty = 5, lwd = 1, col = "grey73")
        }
        
}