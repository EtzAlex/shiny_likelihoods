##

plot.norm <- function(priorMean, priorSD, sampleMean = NULL, sampleSE = NULL, null = NULL){
        
        w <- seq(from=priorMean-10*priorSD, to=priorMean+10*priorSD, by=priorSD/1000) ##Set up for distributions
        y1 <- dnorm(w, priorMean, priorSD) ## data for prior curve
        y2 <- dnorm(w, sampleMean, sampleSE) ## data for likelihood curve
        if(mode(sampleMean) == "numeric"){
                ##Calculate and plot posterior distribution
                priorPrec <- (1/priorSD)^2 ## prior precision, 1 / sigma^2
                dataPrec <- (1/sampleSE)^2 ## data (likelihood) precision, 1 / SE^2
                postPrec <- priorPrec + dataPrec 
                postSD <- sqrt(1/postPrec) ## 1/precision = var
                postMean <- postSD^2 * (priorMean/priorSD^2 + sampleMean/sampleSE^2) ##calculate mean of posterior
                y3 <- dnorm(w, postMean, postSD) #data for posterior line
                
        }
        
        plot(w, y1, xlim= c(priorMean-6*priorSD,priorMean+6*priorSD), ylim=c(0,1.25*max(y1,y2,y3)), type = "l", ylab= "Density",
             xlab= "Mu", las=1, main="Prior-to-Posterior Transformation with Normal Data",lwd=5,
             cex.lab=1.5, cex.main=1.5, col = "skyblue", axes=FALSE)
        axis(1, at = seq(from=priorMean-6*priorSD, to = priorMean+6*priorSD, by = 2*priorSD)) #adds x axis spanning 6 prior SDs
        axis(2, las=1) #adds y axis 
        
        if(mode(sampleMean) == "numeric"){
        lines(w, y2, lwd =3, lty =3, col = "darkorange") #if there is new data, plot the likelihood
        lines(w, y3, lwd = 5, col = "darkorchid1") # if there is new data, plot the posterior
        legend("topright", c("Prior", "Posterior", "Likelihood"), col = c("skyblue", "darkorchid1", "darkorange"),
               lty = c(1,1,3), lwd = c(5,5,3), bty = "n", inset = c(-.3, 0), y.intersp = .55,
               x.intersp = .40)
        }
        
        if(mode(null) == "numeric" && mode(sampleMean) == "numeric"){
                ## Adds points on the distributions at the null value if there is one and if there is new data
                points(null, dnorm(null, priorMean, priorSD), pch = 21, bg = "blue", cex = 1.5) #adds point to prior
                points(null, dnorm(null, postMean, postSD), pch = 21, bg = "darkorchid", cex = 1.5) #adds point to posterior
                abline(v=null, lty = 5, lwd = 1, col = "grey73") ##adds vertical line at null value
                
                ##Calculate BF using Savage-Dickey density ratio
                null.H0 <- dnorm(null, priorMean, priorSD)
                null.H1 <- dnorm(null, postMean, postSD)
                CI.low <- qnorm(.025, postMean, postSD)
                CI.high <- qnorm(.975, postMean, postSD)
                
                return( list("BF01 (in favor of H0)" = null.H1/null.H0, "BF10 (in favor of H1)" = null.H0/null.H1, 
                             "Posterior Mean" = postMean, "Posterior SD" = postSD, 
                             "95% cred lower" = CI.low, "95% cred upper" = CI.high))
        }
}