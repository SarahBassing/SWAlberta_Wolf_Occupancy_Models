################################################################################
#  Southwest Alberta Grey Wolf Patch Occupancy Model
#  Sarah B. Bassing
#  Montana Cooperative Wildlife Research Unit
#  Aug 2016

#  This is the baseline template for future models
#  This is a Dynamic MULTI-SEASON occupancy model
#  This model accounts for FALSE-POSITIVE DETECTIONS in the data
#  There are 9 sampling occasions per season:
#     1st is based on genetic samples collected during rendezvous site surveys
#     2nd is based on hunter observations of live wolves during big game season
#  Encounter histories include:
#     1 = no detection 
#     2 = uncertain detection (two hunter observations of 2+ live wolves)
#     3 = certain detection (three or more hunter observations of 2+ live wolves;
#         genetic detections of two or more wolves)
#  Make sure cov[] is correct for each model
#  Make sure to change sink, out, and write.table text files for each model

################################################################################
# Pull in encounter histories and covariate data from source script
# Data has been organized into appropriate arrays (ABwolf; cov[], survey[])

source("C:/Sarah B/Thesis POM/Scripts/Final_Models/Input_Data.R")

setwd("C:/Sarah B/Thesis POM/Scripts/Final_Models")

library(R2jags)
library(mcmcplots)
################################################################################
#  Specify model in BUGS language

sink("Dynamic_FP_Wolf_base.txt")
cat("
model {

################################################################################
#  Specify priors

#  Everything is on the logit scale- defining priors based on my linar models ##
#  Need to transform them to probability scale later ##

  #  Priors for occupancy parameters
  #  Prior for psi1 coefficient (occupancy in year 1)
  #  Prior truncated

  B0.psi1 ~ dnorm(0,0.001)T(-5,5)	        

  #  Priors for transition probabilities (survival and colonization)
  #  Priors truncated

  for(k in 1:K-1){
    B0.phi[k] ~ dnorm(0,0.001)T(-5,5)
    B0.gamma[k] ~ dnorm(0,0.001)T(-5,5)
  }#k


  # #  Holding occupancy, colonization, and survial parameters constant over time
  # B0.phi.const ~ dnorm(0,0.001)T(-5,5)
  # B0.gamma.const ~ dnorm(0, 0.001)T(-5,5)
  # 
  # for(k in 1:K-1){
  #   B0.phi[k] <- B0.phi.const      # Prior for phi coefficient
  #   B0.gamma[k] <- B0.gamma.const    # Prior for gamma coefficient
  # }#k



  # Priors for detection probabilities (only varies by year)

  for(k in 1:K){
  B0.p11[k] ~ dnorm(0,0.001)
  B0.p10[k] ~ dnorm(0,0.001)
  B0.b[k] ~ dnorm(0,0.001)
  }#k


  # #  Holding detection parameters constant over time
  # B0.p11.const ~ dnorm(0, 0.001)
  # B0.p10.const ~ dnorm(0, 0.001)
  # B0.b.const ~ dnorm(0, 0.001)
  # 
  # for(k in 1:K){
  #   B0.p11[k] <- B0.p11.const
  #   B0.p10[k] <- B0.p10.const
  #   B0.b[k] <- B0.b.const
  # }#k



  #  Priors for covariates
  #  Survey effort (SVY)
  #
  # B7.p11.SVY ~ dnorm(0,0.001)
  # B7.p10.SVY ~ dnorm(0,0.001)
  # B7.b.SVY ~ dnorm(0,0.001)

################################################################################
# Ecological process/submodel
# Define State (z) conditional on parameters- Nmbr sites occupied

	#  logit.psi1 is on the logit scale (psi is the probability of B0.psi + B1.cov...)
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
#  Observation process/submodel
#  z is either 0 or 1 (unoccupied or occupied)
#  y (observation dependent on the state z) can be 0,1,2 (no obs, uncertain obs,
#  certain obs) but JAGS's multinomial link function, dcat(), needs y to be 1,2,3

	#  Observation process: define observations [y,i,j,k,z]
  #  y|z has a probability of...
	#  Detection probabilities are site, occasion, and year specific
 
	for(i in 1:nSite){
    for (j in 1:nOcc){
      for(k in 1:K){			                                  
	  	  p[1,i,j,k,1] <- (1-p10[i,j,k])
 	  	  p[1,i,j,k,2] <- (1-p11[i,j,k])
  	    p[2,i,j,k,1] <- p10[i,j,k]
  	    p[2,i,j,k,2] <- (1-b[i,j,k])*p11[i,j,k]
  	    p[3,i,j,k,1] <- 0
        p[3,i,j,k,2] <- b[i,j,k]*p11[i,j,k]
		  }#k
	  }#j
  }#i

  
# # Need mulitnomial link function for false positive detections: dcat function in JAGS
# 
# 	#  Observation model
#                                          
# 	for(i in 1:nSite){
#   	for(j in 1:nOcc){
#       for(k in 1:K){
#         multilogit.p11[i,j,k] <- B0.p11[k]+B7.p11.SVY*SVY[i,j,k]
#         logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
#         multilogit.p10[i,j,k] <- B0.p10[k]+B7.p10.SVY*SVY[i,j,k]
#         logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
#         multilogit.b[i,j,k] <- B0.b[k]+B7.b.SVY*SVY[i,j,k]
#         logit(b[i,j,k]) <- multilogit.b[i,j,k]
#     	  y[i,j,k] ~ dcat(p[,i,j,k,z[i,k]+1])
# 		  }#k
#   	}#j
# 	}#i

  #  Observation model without covariates on detection probabilities
	for(i in 1:nSite){
  	for(j in 1:nOcc){
      for(k in 1:K){
        multilogit.p11[i,j,k] <- B0.p11[k]
        logit(p11[i,j,k]) <- multilogit.p11[i,j,k]
        multilogit.p10[i,j,k] <- B0.p10[k]
        logit(p10[i,j,k]) <- multilogit.p10[i,j,k]
        multilogit.b[i,j,k] <- B0.b[k]
        logit(b[i,j,k]) <- multilogit.b[i,j,k]
    	  y[i,j,k] ~ dcat(p[,i,j,k,z[i,k]+1])
		  }#k
  	}#j
	}#i

################################################################################
#  Derived Parameters    

	#  Derived parameters for occupancy and transitions
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

  #  Mean Annual Parameter           
  
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
     ik.p11[j,k] <- mean(p11[,j,k])
     ik.p10[j,k] <- mean(p10[,j,k])
     ik.b[j,k] <- mean(b[,j,k])
    }#j
  }#k
  
  #  Mean detection across all sites and occasions by survey method per year
  #  Don't need to estimate k.p10.rnd[k] b/c no FP detections in RND surveys
  
  for(k in 1:K){		                                       
    k.p11.rnd[k] <- ik.p11[1,k]                       
    k.p11.hs[k] <- mean(ik.p11[2:9,k])
    k.p10.hs[k] <- mean(ik.p10[2:9,k])                     
    k.b.rnd[k] <- ik.b[1,k]                           
    k.b.hs[k] <- mean(ik.b[2:9,k])
  }#k
 
}
", fill=TRUE)
sink()

################################################################################
#  Bundle data, specify parameters & MCMC settings, and run JAGS

# Define & Bundle Data
Site <- length(unique(dynamic.ABwolf$cell_ID))
win.data <- list("y"=ABwolf, "nSite"=Site, "nOcc"=9, "K"=3, "SVY"=survey)


#  Initial Values	
#  PAY ATTENTION to this if format of encounter histories change
#  Use naive occupancy estimate as initial value 

#  Calculating naive occupancy based on raw encouter histories
#  Necessary when encouter histories are in 1,2,3 format for JAGS
zst <- apply(ABwolf,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
zst[zst=="-Inf"] <- NA 	# Change -Inf backs back to NAs (weird apply fucntion issue)
zst[zst==1] <- 0  			# Makes 1's => 1
zst[zst==3] <- 1  			# Makes 3's => 1
zst[zst==2] <- 1  			# Makes 2's => 1

inits <- function(){list(z=zst)}
# inits <- function(){list(z=zst, B0.psi1=0.5, B0.gamma.const=-2, B0.phi.const=0,
#                          B0.p11.const=0, B0.p10.const=-0.5, B0.b.const=0.5)}


# Parameters to keep track of and report
params <- c("B0.psi1", "B0.phi", "B0.gamma", "B0.p11", "B0.p10", "B0.b",
            "psik", "gammak", "phik", "growthrk", "turnoverk", "n.occ",
            "k.p11.rnd", "k.b.rnd","k.p11.hs","k.p10.hs", "k.b.hs") 
            # include "psi" if want site specific probabilities of occupancy
            # include "B7.p11.SVY", "B7.p10.SVY", "B7.b.SVY" if includig survey-specific covariates


# MCMC Settings 
ni <- 300000
nt <- 4
nb <- 150000
nc <- 3


# Call JAGS 
out <- jags(win.data, inits, params, "Dynamic_FP_Wolf_base.txt", n.chains=nc, 
            n.thin=nt, n.iter=ni, n.burnin=nb, jags.module = c("glm","dic"))

jag.sum <- out$BUGSoutput$summary
write.table(x=jag.sum, file="C:/Sarah B/Thesis POM/Model_Outputs/Dynamic_FP_Wolf_base.txt", sep="\t")

print(out, dig=2)

out2 <- out
#mcmcplot(out2)
