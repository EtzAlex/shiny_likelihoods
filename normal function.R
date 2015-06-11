## This is a simple function for calculating and visualizing bayesian inference
## The inference method is bayes factors by savage-dickey density ratio
## The function works by inputting the various distribution hyperparameters and data summaries 
## The variables are all pretty straightforward, but I'll define them anyway
## priorMean is the location parameter (mu) for the normal (gaussian) prior distribution
## priorSD is the scale parameter (sigma) for the prior distribution
## sampleMean is the mean value calculated from sample data
## sampleSE is the standard **ERROR** calculated from sample data
## Remember, t is the sampleMean divided by sampleSE, so if you know any two you know the other
## CI is the credible interval limits you want reported and shaded. .95 is recommended and standard
## null is the null value of interest (if left NULL it won't calculate a bayes factor)
## this function returns plots, bayes factors, posterior distribution hyperparameters, and 95% credible interval limits


## The function
plot.norm <- function(priorMean = 0, priorSD = 1, sampleMean = NULL, sampleSE = NULL, CI = NULL, null = NULL){
        
        w <- seq(from=priorMean-15*priorSD, to=priorMean+15*priorSD, by=priorSD/1000) ##Set up for distributions
        y1 <- dnorm(w, priorMean, priorSD) ## data for plotting prior curve
        
        if(is.numeric(sampleMean) == T){
                y2 <- dnorm(w, sampleMean, sampleSE) ## data for plotting likelihood curve
                
                ##Calculate and prepare to plot posterior distribution
                priorPrec <- (1/priorSD)^2 ## prior precision, 1 / sigma^2
                dataPrec <- (1/sampleSE)^2 ## data (likelihood) precision, 1 / SE^2
                postPrec <- priorPrec + dataPrec 
                postSD <- sqrt(1/postPrec) ## 1/precision = var
                postMean <- postSD^2 * (priorMean/priorSD^2 + sampleMean/sampleSE^2) ##calculate mean of posterior
                y3 <- dnorm(w, postMean, postSD) #data for posterior line
                
        }
        
        
        ##Defining y axis limits and ticks
        ##If there is sample data, base graphic dimensions off of posterior scale
        if(is.numeric(sampleMean) == T){
                yupper <- 1.25 * max(y1,y2,y3)
                ylower <- 0
                xlower <- qnorm(.005, postMean, postSD)
                xupper <- qnorm(.995, postMean, postSD)
        }
        
        ##If there is no sample data, base dimensions off of prior scale
        else{
                yupper <- 1.25 * max(y1)
                ylower <- 0
                xlower <- qnorm(.001,priorMean,priorSD)
                xupper <- qnorm(.999,priorMean,priorSD)
        }
        
        xlims <- pretty(c(xlower,xupper)) ##To be used for x axis ticks
        ylims <- pretty(c(ylower,yupper)) ##To be used for y axis ticks
        xrange <- max(xlims)-min(xlims)
        
        ##Plot only the prior distribution if there is no sample data. Notice the titles are different
        if(is.numeric(sampleMean) == F){
        plot(w, y1, xlim= c(min(xlims),max(xlims)), ylim=c(ylower,max(ylims)), type = "l", ylab= "Density",
             xlab= "Mu", las=1, main="Prior distribution with normal shape",lwd=5,
             cex.lab=1.5, cex.main=1.5, col = "skyblue",axes=F)
        
        #Format the axes so they look good
        axis(1, at = xlims, pos= -.02 * (max(ylims))) #adds x axis 
        axis(2, at = ylims, las=1, pos= min(xlims)- abs(.02*xrange ))#adds y axis 
               
        }
        
        ##Plot all three distributions if there is sample data. Notice the titles and x axis caption are different
        if(is.numeric(sampleMean) == T){
                plot(w[w > (min(xlims) - .02 * xrange)], y1[w > (min(xlims) - .02 * xrange)], xlim= c(min(xlims),max(xlims)), 
                ylim=c(0,max(ylims)), type = "l", ylab= "Density", xlab= "Mu", las=1,lwd=3, lty = 2,
                cex.lab=1.5, cex.main=1.5, col = "skyblue",axes=F)
                title(main="Prior-to-Posterior Transformation", sub = "(normal prior with normal data)")
                
                #Format the axes so they look good
                axis(1, at = xlims, pos= -.02 * (max(ylims))) #adds x axis 
                axis(2, at = ylims, las=1, pos= min(xlims)- abs(.02* xrange ))#adds y axis
                
                ##Add lines for likelihood and posterior
                lines(w[w> (min(xlims) - .02 * xrange)], y2[w> (min(xlims) - .02 * xrange)], lwd =2, lty =3, col = "darkorange") #if there is sample data, plot the likelihood
                lines(w[w> (min(xlims) - .02 * xrange)], y3[w> (min(xlims) - .02 * xrange)], lwd = 4, col = "darkorchid1") # if there is sample data, plot the posterior
                
                ##Add legend in the top left
                legend("topleft", c("Prior", "Posterior", "Likelihood"), col = c("skyblue", "darkorchid1", "darkorange"), 
                       lty = c(2,1,3), lwd = c(3,5,2), bty = "n", y.intersp = .55, x.intersp = .1, seg.len=.7)
        }
        
        ##If there is sample data and a null value specified, but no CI cutoff then add points and return BFs
        if(is.numeric(null) == T && is.numeric(sampleMean) == T && is.numeric(CI) == F){
                ## Adds points on the distributions at the null value 
                points(null, dnorm(null, priorMean, priorSD), 
                       pch = 21, bg = "blue", cex = 1.5) #adds point to prior
                points(null, dnorm(null, postMean, postSD), 
                       pch = 21, bg = "darkorchid", cex = 1.5) #adds point to posterior
                lines(c(0,0),c(0,1.15*max(y1,y2,y3)), lwd=1, lty = 5, col ="grey73")##adds vertical line at null value
                ##Calculate BF using Savage-Dickey density ratio
                null.H0 <- dnorm(null, priorMean, priorSD)
                null.H1 <- dnorm(null, postMean, postSD)
                return( list("BF01 (in favor of H0)" = round(null.H1/null.H0,2), "BF10 (in favor of H1)" = round(null.H0/null.H1,2), 
                             "Posterior Mean" = round(postMean,2), "Posterior SD" = round(postSD,2)))
        }
                
        ##If there is sample data, a CI cutoff, but no null value, only return posterior descriptives and shade graph
        if(is.numeric(sampleMean) == T && is.numeric(null) == F && is.numeric(CI) == T){
                ##Calculate CI limits
                CI.low <- qnorm((1-CI)/2, postMean, postSD) 
                CI.high <- qnorm(1-(1-CI)/2, postMean, postSD)
                
                ##Adds shading for the area beyond the curve
                CI.lowtail <- mean(xlims) - abs(.55* xrange) ##Set up the max range of the shading
                CI.hightail <- mean(xlims) + abs(.55 * xrange)
                SEQlow<-seq(CI.lowtail, CI.low, (CI.low - CI.lowtail)/1000)
                SEQhigh <- seq(CI.hightail,CI.high, (CI.high - CI.hightail)/1000)
                ##Adds shaded area for x% Posterior CIs
                cord.x <- c(CI.lowtail, SEQlow, CI.low) ##set up for shading
                cord.y <- c(0,dnorm(SEQlow,postMean, postSD),0) ##set up for shading
                polygon(cord.x,cord.y,col='orchid', lty= 3) ##shade left tail
                cord.xx <- c(CI.hightail, SEQhigh,CI.high) 
                cord.yy <- c(0,dnorm(SEQhigh,postMean, postSD), 0)
                polygon(cord.xx,cord.yy,col='orchid', lty=3) ##shade right tail
                
                return( list("Posterior Mean" = round(postMean,2), "Posterior SD" = round(postSD,2), 
                             "Posterior 95% CI lower" = round(CI.low,2), "Posterior 95% CI upper" = round(CI.high,2)))
        }
        
        ##If there is sample data and a null value, add null line, null points,
        ##and return posterior descriptives and Savage-Dickey bayes factors
        if(is.numeric(null) == T && is.numeric(sampleMean) == T && is.numeric(CI) == T){
                ##Calculate CI limits
                CI.low <- qnorm((1-CI)/2, postMean, postSD) 
                CI.high <- qnorm(1-(1-CI)/2, postMean, postSD)
                
                ##Adds shading for the area beyond the curve
                CI.lowtail <- mean(xlims) - abs(.55* xrange) ##Set up the max range of the shading
                CI.hightail <- mean(xlims) + abs(.55 * xrange)
                SEQlow<-seq(CI.lowtail, CI.low, (CI.low - CI.lowtail)/1000)
                SEQhigh <- seq(CI.hightail,CI.high, (CI.high - CI.hightail)/1000)
                ##Adds shaded area for x% Posterior CIs
                cord.x <- c(CI.lowtail, SEQlow, CI.low) ##set up for shading
                cord.y <- c(0,dnorm(SEQlow,postMean, postSD),0) ##set up for shading
                polygon(cord.x,cord.y,col='orchid', lty= 3) ##shade left tail
                cord.xx <- c(CI.hightail, SEQhigh,CI.high) 
                cord.yy <- c(0,dnorm(SEQhigh,postMean, postSD), 0)
                polygon(cord.xx,cord.yy,col='orchid', lty=3) ##shade right tail
                
                ## Adds points on the distributions at the null value 
                points(null, dnorm(null, priorMean, priorSD), 
                       pch = 21, bg = "blue", cex = 1.5) #adds point to prior
                points(null, dnorm(null, postMean, postSD), 
                       pch = 21, bg = "darkorchid", cex = 1.5) #adds point to posterior
                lines(c(null,null),c(-.02 * (max(ylims)),1.11*max(y1,y2,y3)), lwd=1, lty = 5, col ="grey73")##adds vertical line at null value
                ##Calculate BF using Savage-Dickey density ratio
                null.H0 <- dnorm(null, priorMean, priorSD)
                null.H1 <- dnorm(null, postMean, postSD)
                
                return( list("BF01 (in favor of H0)" = round(null.H1/null.H0,2), "BF10 (in favor of H1)" = round(null.H0/null.H1,2), 
                             "Posterior Mean" = round(postMean,2), "Posterior SD" = round(postSD,2), 
                             "Posterior 95% CI lower" = round(CI.low,2), "Posterior 95% CI upper" = round(CI.high,2)))
        }
       
}