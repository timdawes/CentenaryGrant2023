# Understanding the determinants of Low Cardiac Output Syndrome (LCOS) after paediatric cardiac surgery
# Tim Dawes
# May 2023

# Code to simulate sample sizes:

library(colorRamps)
library(readxl)
library(Hmisc)
library(RColorBrewer)

cols<- c(brewer.pal(3,"Reds")[-1], "black", rev(brewer.pal(3, "Blues")[-1]))


# Plot the data

      Rs<- seq(0.8,1.2,length.out=5)

# Function to generate data for groups of size N
      corr_func <- function(N, R) {
        age <- rbeta(N,2,20)*50
        age.beta <- rnorm(N, 0.103, 0.217 / R)  # Age
        
        sats <- rnorm(N, 95, 10) # Age
        sats[sats>100]<- 100
        sats.beta <- rnorm(N, 0.127, 0.220 / R)  # Sats
        
        shunt <- rbinom(N, 1, 0.33) 
        shunt.beta <- rnorm(N, 0.131, 0.127 / R)  # Shunt
        
        CPBduration <- rbeta(N,4,8)*200
        CPBduration.beta <- rnorm(N, 0.146, 0.170 / R)  # CPB duration
        
        residualshunt <- rbinom(N, 1, 0.1)
        residualshunt.beta <- rnorm(N, 0.109, 0.113 / R) 
        
        LCOS<- age*age.beta + sats*sats.beta + shunt*shunt.beta + CPBduration*CPBduration.beta + residualshunt*residualshunt.beta
        
        d<- data.frame(age=age, sats=sats, shunt=shunt, CPBduration=CPBduration, residualshunt=residualshunt, LCOS=LCOS)
        
        fit<- lm(LCOS ~ age + sats + shunt + CPBduration + residualshunt, d, family=binomial)
        f<- summary(fit)$fstatistic
        p<- pf(f[1], f[2], f[3], lower.tail=F)
        
        return(c(f=f[1], p=p, sig=(p < .05)))
        # return a named vector with the results we want to keep
      }


# Loops to check whether effect is detected in successive trials
      df<- list()
      
      par(mfrow=c(1,1), mar=c(5,6,2,3), mgp=c(1,1,1))
      calc<- TRUE
      maxN<- 1500
      maxTrials<- 300
      sample.size<- rep(0, length(Rs))
      options(warn=-1)
      s<- seq(2,maxN, length.out=50)
      
      for (k in 1:length(Rs))
      {
        cat(k)
        
      sig<- rep(0, length(s))
      
      if (calc==TRUE){
            for (i in s)
            {
              cat(".")
              for (j in 1:maxTrials)
              {
              corr.results<- corr_func(i,Rs[k])
              sig[match(i,s)]<- sig[match(i,s)] + corr.results[3] 
              }
            }
        df[[k]]<- data.frame(X=s, Y= 100 * sig / maxTrials)
      }
        
        if (k==1) {plot(df[[k]]$X, df[[k]]$Y, type='p', lwd=1, col=cols[k], pch=19, cex=1, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', xlim=c(0,maxN), ylim=c(0,100))}
        if (k>1) {for (l in 2:k) {points(df[[l]]$X, df[[l]]$Y, type='p', lwd=1, col=cols[l], pch=19, cex=1, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', xlim=c(0,100), ylim=c(0,100))}}
        
        axis(side=1, lwd=5, at=seq(0,maxN,length.out=7), line=0, cex.axis=2, font=2)
        mtext("Number of patients in study", side=1, line=3.5, cex=3, font=2)
        axis(side=2, lwd=5, at=seq(0,100,length.out=5), line=-1, cex.axis=2, font=2, las=2)
        mtext("Power (%)", side=2, line=3, cex=3, font=2)
        lo<- loess(Y~X, span=0.75, data=df[[k]])
        sample.size[k]<- with(df[[1]], which.min(abs(predict(lo,X)-90)))
        
        new.y<- with(df[[k]], predict(lo, X))
        new.y[new.y>100]<- 100
        new.x<- df[[k]]$X
        points(new.x, new.y, type='l', lwd=7, col=cols[k])
        
      } 
      
      
# Sample size needed in each group for 90% power
      
      segments(40,90,1450,90,lwd=5,col="black", lty=2)
      segments(s[sample.size[3]],0,s[sample.size[3]],100,lwd=5,col="black", lty=2)
      
      text(230, 95, "90% Power", cex=2.2, font=2)
      
      step<- 12
      text(1100, (6*step)-5, "Standard deviation", cex=2.4, font=2)
      
      for (k in 1:length(Rs))
          {
          segments(750,(k*step)-5,1450,(k*step)-5, lwd=40, col=cols[k])
          t<- c("Standard deviation 120%", "Standard deviation 110%", "Standard deviation 100%", "Standard deviation 90%", "Standard deviation 80%")[k]
          text(1100,(k*step)-5, labels=t, col="white", cex=1.8)
      }
      
      # Export as 5" x 10"
      