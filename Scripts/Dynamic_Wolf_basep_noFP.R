#########################################################################################

#  Southwest Alberta Grey Wolf Patch Occupancy Model

#  Make sure EH and survey data are set up correctly before running this
#  This is the baseline template to for future models
#  This model does not account for false positive detections
#  This is the base model with SURVEY on detection
#  Make sure cov[] is correct for each model
#  Make sure to change sink, out, and write.table text files for each model

################################################################################
#  Data and packages
#  Pull in encounter histories and covariate data from source script
#  Data has been organized into appropriate arrays (ABwolf; cov[], survey)

source("C:/Sarah B/Thesis POM/Scripts/Final_Models/No_FP_Final_Models/Input_Data.R")

setwd("C:/Sarah B/Thesis POM/Scripts/Final_Models/No_FP_Final_Models")

library(R2jags)
library(mcmcplots)

################################################################################
#  Specify model in BUGS language

sink("Dynamic_Wolf_basep_noFP.txt")
cat("
model {

################################################################################
#  Specify priors
#  Everything is on the logit scale- defining priors based on my linar models
#  Need to transform to probability scale later

  # Prior for detection probability (varies by year and survey method)

  for(k in 1:K){
    B0.p.rend[k] ~ dnorm(0,0.001)
    B0.p.hunt[k] ~ dnorm(0,0.001)
  }#k

  # Priors for occupancy parameters
  B0.psi1 ~ dnorm(0,0.001)T(-5,5)       # Prior for psi1 coefficient

  # B0.phi.const ~ dnorm(0, 0.001)T(-5,5)
  # B0.gamma.const ~ dnorm(0, 0.001)T(-5,5)
  # 
  # for(k in 1:K-1){
  #   B0.phi[k] <- B0.phi.const      # Prior for phi coefficient
  #   B0.gamma[k] <- B0.gamma.const    # Prior for gamma coefficient
  # }#k

  for(k in 1:K-1){
    B0.phi[k] ~ dnorm(0,0.001)T(-5,5)      # Prior for phi coefficient
    B0.gamma[k] ~ dnorm(0,0.001)T(-5,5)    # Prior for gamma coefficient
  }#k

  # Priors for covariates
  # B7.p.SVY ~ dnorm(0,0.001)       # Survey effort (rnd & hunter survey) on detection


################################################################################
#  Ecological process/submodel
#  Define State conditional on parameters- No. sites occupied

#  logit.psi1 is on the logit scale (psi is the probability of B0.psi + B1.psi.FC...)
#  logit() transforms the logit scale to the probability scale
#  gamma/phi are the transition probabilities FOR year 1 to the next
  
  for(i in 1:nSite){
	  logit.psi1[i] <- B0.psi1       
    logit(psi1[i]) <- logit.psi1[i]                                     
    z[i,1] ~ dbern(psi1[i])                                             
	  for(k in 1:K-1){                                                    
      logit.phi[i,k] <- B0.phi[k]                                       
      logit.gamma[i,k] <- B0.gamma[k]
      logit(phi[i,k]) <- logit.phi[i,k]                         
      logit(gamma[i,k]) <- logit.gamma[i,k]
  	 }#k
      for(k in 2:K){
        muZ[i,k] <- z[i,k-1]*phi[i,k-1] + (1-z[i,k-1])*gamma[i,k-1]
    	  z[i,k] ~ dbern(muZ[i,k])
  	  }#k
	}#i


################################################################################
#  Observation process/submodel: indexed under one for-loop

#  z is either 0 or 1 (not occupied or occupied) but JAGS needs it to be 1 or 2
#  y (observation dependent on the state z) can be 0,1 (no obs or obs)
#  but needs to be 1,2 for JAGS

#  logit() transfroms from logit scale to probability scale

  for(i in 1:nSite){
      for(k in 1:K){
        logit.p[i,1,k] <- B0.p.rend[k]
        for(j in 2:nOcc){
          logit.p[i,j,k] <- B0.p.hunt[k]
        }#j
        for(j in 1:nOcc){
          logit(p[i,j,k]) <- logit.p[i,j,k]  
          muy[i,j,k] <- z[i,k]*p[i,j,k]
          y[i,j,k] ~ dbern(muy[i,j,k])
        }#j
      }#k
    }#i


################################################################################
#  Derived Parameters   

#  Must define growthr in yr1, otherwise growthr[k] in yr2 has nothing to work with
	
    for(i in 1:nSite){
      psi[i,1] <- psi1[i]
      growthr[i,1] <- 1                                        
      for (k in 2:K){                                          
        psi[i,k] <- psi[i,k-1]*phi[i,k-1] + (1-psi[i,k-1])*gamma[i,k-1]
        growthr[i,k] <- psi[i,k]/psi[i,k-1]
        turnover[i,k-1] <- (1 - psi[i,k-1]) * gamma[i,k-1]/psi[i,k]
      }#k
    }#i


################################################################################
# Annual parameter averages       
  
  for(k in 1:K){
    psik[k] <- mean(psi[,k])
    n.occ[k] <- sum(z[1:nSite, k])
  }#k

  for(k in 1:K-1){				
    gammak[k] <- mean(gamma[,k])
    phik[k] <- mean(phi[,k])
    turnoverk[k] <- mean(turnover[,k])
  }#k

  growthrk[1] <- 1   
  for(k in 2:K){
    growthrk[k] <- mean(growthr[,k])
  }#k


#  Mean detection across sites for each occasion and year
  
  for(j in 1:nOcc){	                                        
    for(k in 1:K){
     ik.p[j,k] <- mean(p[,j,k])
    }#j
  }#k
  

#  Mean detection across all sites and occasions by survey method per year
#  These will be the same without survey covariates
  
  for(k in 1:K){		                                       
    k.p.rnd[k] <- ik.p[1,k]                       
    k.p.hs[k] <- mean(ik.p[2:9,k])
  }#k
 
}
", fill=TRUE)
sink()


############# Bundle data, specify parameters & MCMC settings, and run JAGS ###############

# Define & Bundle Data
Site <- length(unique(dynamic.ABwolf$cell_ID))
win.data <- list("y"=ABwolf, "nSite"=Site, "nOcc"=9, "K"=3)


# Initial Values	### PAY ATTENTION to this if I change the format of my encounter histories!!!
## Use naive occupancy estimate as initial value ## 

# Calculating naive occupancy based on raw encouter histories
# When encouter histories are in 1,2,3 format for JAGS
zst <- apply(ABwolf,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Apply function turns NAs into -Inf so have to change it back
zst[zst==0] <- 0  			# Makes 1's => 0 for counting naive occ (use for initial z value)
zst[zst==1] <- 1  			# Makes 2's => 1 for counting naive occ (use for initial z value)


inits <- function(){list(z=zst)}

# Parameters to keep track of and report
params <- c("psik", "gammak", "phik", "growthrk", "turnoverk", "n.occ", "k.p.rnd",
            "k.p.hs", "B0.psi1", "B0.phi", "B0.gamma", "B0.p.rend", "B0.p.hunt", "psi")

# MCMC Settings 
ni <- 300
nt <- 4
nb <- 150
nc <- 3

# Call JAGS 
out <- jags(win.data, inits, params, "Dynamic_Wolf_basep_noFP.txt", n.chains=nc, 
            n.thin=nt, n.iter=ni, n.burnin=nb, jags.module = c("glm","dic"))

jag.sum <- out$BUGSoutput$summary
#write.table(x=jag.sum, file="C:/Sarah B/Thesis POM/Model_Outputs/Dynamic_Wolf_basep_noFP.txt", sep="\t")

print(out, dig=2)

out2 <- out
#mcmcplot(out2)




