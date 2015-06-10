## This is a simple function for calculating and visualizing bayesian inference
## The inference method is bayes factors by savage-dickey density ratio
## The function works by inputting the various distribution hyperparameters and data summaries 
## The variables are all pretty straightforward, but I'll define them anyway
## priorMean is the location parameter (mu) for the normal (gaussian) prior distribution
## priorSD is the scale parameter (sigma) for the prior distribution
## sampleMean is the mean value calculated from sample data
## sampleSE is the standard **ERROR** calculated from sample data
## Remember, t is the sampleMean divided by sampleSE, so if you know any two you know the other
## null is the null value of interest (if left NULL it won't calculate a bayes factor)
## This function returns plots, bayes factors, posterior distribution hyperparameters, and 95% credible interval limits

## The function
plot.norm <- function(priorMean, priorSD, sampleMean = NULL, sampleSE = NULL, null = NULL){
        
        w <- seq(from=priorMean-10*priorSD, to=priorMean+10*priorSD, by=priorSD/1000) ##Set up for distributions
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
        
        ##Defining x axis limits and ticks
        xlower <- qnorm(.01,priorMean,priorSD)
        xupper <- qnorm(.99,priorMean,priorSD)
        xlims <- pretty(c(xlower,xupper)) ##To be used for x axis ticks
        
        ##Defining y axis limits and ticks
        ##If there is sample data
        if(is.numeric(sampleMean) == T){
                yupper <- 1.25 * max(y1,y2,y3)
                ylower <- 0
        }
        ##If there is no sample data
        else{
                yupper <- 1.25 * max(y1)
                ylower <- 0
        }
        ylims <- pretty(c(ylower,yupper)) ##To be used for y axis ticks
        
        ##Plot only the prior distribution if there is no sample data. Notice the titles are different
        if(is.numeric(sampleMean) == F){
        plot(w, y1, xlim= c(min(xlims),max(xlims)), ylim=c(ylower,max(ylims)), type = "l", ylab= "Density",
             xlab= "Mu", las=1, main="Prior distribution with normal shape",lwd=5,
             cex.lab=1.5, cex.main=1.5, col = "skyblue",axes=F)
        
        #Format the axes so they look good
        axis(1, at = xlims, pos= -.02 * (max(ylims))) #adds x axis 
                if(min(xlims) < 0){ 
                        axis(2, at = ylims, las=1, pos= min(xlims)+.06*min(xlims) )#adds y axis 
                }
                if(min(xlims) > 0){
                        axis(2, at = ylims, las=1, pos= min(xlims)-.06*min(xlims) )#adds y axis 
                }
        }
        
        ##Plot all three distributions if there is sample data. Notice the titles and x axis caption are different
        if(is.numeric(sampleMean) == T){
                plot(w, y1, xlim= c(min(xlims),max(xlims)), ylim=c(0,max(ylims)), type = "l", ylab= "Density",
                xlab= "Mu", las=1,lwd=5,
                cex.lab=1.5, cex.main=1.5, col = "skyblue",axes=F)
                title(main="Prior-to-Posterior Transformation", sub = "(normal prior with normal data)")
                
                #Format the axes so they look good
                axis(1, at = xlims, pos= -.02 * (max(ylims))) #adds x axis 
                if(min(xlims) < 0){
                        axis(2, at = ylims, las=1, pos= min(xlims)+.06*min(xlims) )#adds y axis 
                }
                if(min(xlims) > 0){
                        axis(2, at = ylims, las=1, pos= min(xlims)-.06*min(xlims) )#adds y axis 
                }
                lines(w, y2, lwd =3, lty =3, col = "darkorange") #if there is sample data, plot the likelihood
                lines(w, y3, lwd = 5, col = "darkorchid1") # if there is sample data, plot the posterior
                legend("topleft", c("Prior", "Posterior", "Likelihood"), col = c("skyblue", "darkorchid1", "darkorange"), 
                       lty = c(1,1,3), lwd = c(5,5,3), bty = "n", y.intersp = .55, x.intersp = .1, seg.len=.7)
        }
        
        ##If there is sample data but no null value, only return posterior descriptives
        if(mode(sampleMean) == "numeric" && is.numeric(null) == F){
                ##Calculate CI limits
                CI.low <- qnorm(.025, postMean, postSD) 
                CI.high <- qnorm(.975, postMean, postSD)
                
                return( list("Posterior Mean" = round(postMean,2), "Posterior SD" = round(postSD,2), 
                             "95% cred lower" = round(CI.low,2), "95% cred upper" = round(CI.high,2)))
        }
        ##If there is sample data and a null, add null line, null points,
        ##and return posterior descriptives and Savage-Dickey bayes factors
        if(mode(null) == "numeric" && mode(sampleMean) == "numeric"){
                ## Adds points on the distributions at the null value if there is one and if there is sample data
                points(null, dnorm(null, priorMean, priorSD), 
                       pch = 21, bg = "blue", cex = 1.5) #adds point to prior
                points(null, dnorm(null, postMean, postSD), 
                       pch = 21, bg = "darkorchid", cex = 1.5) #adds point to posterior
                abline(v=null, lty = 5, lwd = 1, col = "grey73") ##adds vertical line at null value
                
                ##Calculate BF using Savage-Dickey density ratio
                null.H0 <- dnorm(null, priorMean, priorSD)
                null.H1 <- dnorm(null, postMean, postSD)
                
                ##Calculate CI limits
                CI.low <- qnorm(.025, postMean, postSD) 
                CI.high <- qnorm(.975, postMean, postSD)
                
                return( list("BF01 (in favor of H0)" = round(null.H1/null.H0,2), "BF10 (in favor of H1)" = round(null.H0/null.H1,2), 
                             "Posterior Mean" = round(postMean,2), "Posterior SD" = round(postSD,2), 
                             "95% cred lower" = round(CI.low,2), "95% cred upper" = round(CI.high,2)))
        }
       
}