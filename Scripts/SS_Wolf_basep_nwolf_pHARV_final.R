################################################################################
  #  Southwest Alberta Grey Wolf Patch Occupancy Model
  #  Sarah B. Bassing
  #  Montana Cooperative Wildlife Research Unit
  #  Aug 2016
  
  #  This is a SINGLE SEASON occupancy model
  #  This model does NOT account for false positive detections in the data
  #  There are 9 sampling occasions per season:
  #     1st is based on genetic samples collected during rendezvous site surveys
  #     2nd is based on hunter observations of live wolves during big game season
  #  Encounter history includes ALL detections: 
  #     uncertain and certain detections are considered certain detections
  #     (i.e., two or more observations of 2+ live wolves = detection of a pack)
  #  Detection prob is parsed out by method (rendezvous site vs hunter surveys)
  #  Covariates on occupancy: percent forest cover (FC)
  #  Make sure cov[] is correct for each model
  #  Make sure to change sink, out, and write.table text files for each model
  
################################################################################
  #  Data and packages
  #  Pull in encounter histories and covariate data from source script
  #  Data has been organized into appropriate arrays (ABwolf; cov[], survey)
  
  source("C:/Sarah B/Thesis POM/Scripts/Final_Models/SS_Models/Input_Data.R")
  
  setwd("C:/Sarah B/Thesis POM/Final_Models_4_Thesis")
  
  library(R2jags)
  library(mcmcplots)
  
################################################################################
  #  Specify model in BUGS language
  
  sink("SS_Wolf_basep_nwolf_pHARV_final.txt")
  cat("
  model {
  

  #  Specify priors
  #  Everything is on the logit scale- defining priors based on my linar models
  #  Need to transform to probability scale later
  
    #  Prior for occupancy and detection probability
    #  Prior on psi is truncated
    #  Detection probability varies by year and is split out by survey method
    for(k in 1:K){
      B0.psi[k] ~ dnorm(0,0.001)T(-10,10)
      B0.p.rend[k] ~ dnorm(0,0.001)
      B0.p.hunt[k] ~ dnorm(0,0.001)
    }#k
  
    #  Prior for covariates
    #B1.FC ~ dnorm(0,0.001)T(-5,5)
    B1.pHvMed ~ dnorm(0,0.001)
    B1.pHvHi ~ dnorm(0,0.001)
  

  #  Ecological process/submodel
  #  Define True State (z) conditional on parameters- Nmbr. sites occupied
  
  #  logit.psi is on the logit scale (psi is the probability of B0.psi + B1.cov...)
  #  logit() transforms the logit scale to the probability scale
  
    for(i in 1:nSite){
       for (k in 1:K){
        logit.psi[i,k] <- B0.psi[k] # + B1.FC*FC[i,k]
        logit(psi[i,k]) <- logit.psi[i,k]
        z[i,k] ~ dbern(psi[i,k])
      }#k
    }#i
  

  #  Observation process/submodel: indexed under one for-loop
  
  #  z is either 0 or 1 (not occupied or occupied)
  #  y (observation dependent on the state z) can be 0,1 (no obs or obs)
  #  logit() transfroms from logit scale to probability scale
  
    for(i in 1:nSite){
      for(k in 1:K){
        logit.p[i,1,k] <- B0.p.rend[k] + B1.pHvMed*HARV[i,k] + B1.pHvHi*HARV2[i,k]
        for(j in 2:nOcc){
          logit.p[i,j,k] <- B0.p.hunt[k] + B1.pHvMed*HARV[i,k] + B1.pHvHi*HARV2[i,k]
        }#j
        for(j in 1:nOcc){
          logit(p[i,j,k]) <- logit.p[i,j,k]  
          muy[i,j,k] <- z[i,k]*p[i,j,k]
          y[i,j,k] ~ dbern(muy[i,j,k])
        }#j
      }#k
    }#i
  

  #  Derived parameters
  #  Annual parameter averages  
  #  n.occ sums occupancy across sites for a total nmbr of occupied sites per year
    
    for(k in 1:K){
      n.occ[k] <- sum(z[1:nSite, k])
      k.psi[k] <- mean(psi[,k])
    }#k
  
  #  Mean detection across sites for each occasion and year
  
    for(j in 1:nOcc){	                                        
      for(k in 1:K){
       ik.p[j,k] <- mean(p[,j,k])
      }#j
    }#k
  
  #  Mean detection across all sites and occasions per year by survey method 
  
    for(k in 1:K){
      k.p.rnd[k] <- ik.p[1,k]
      k.p.hs[k] <- mean(ik.p[2:9,k])
    }#k
  
  #  Annual pack (n.packs) and wolf (n.wolves) abundance, 
  #  area occupied (occ.area, ann.occ, and mean.occ), and mean # packs across 
  #  years (mean.pack)

  for(i in 1:nSite){
    for(k in 1:K){
      occ.area[i,k] <- psi[i,k]*area[i,k]
    }#j
  }#k

  for(k in 1:K){
    ann.occ[k] <- sum(occ.area[,k])
  }#k

  mean.occ <- mean(ann.occ[])

  for(k in 1:K){
      n.packs[k] <- (sum(occ.area[,k]))/mu.terr
    }#k

  mean.pack <- mean(n.packs[])

  for(k in 1:K){
    n.wolves[k] <- ((sum(occ.area[,k]))/mu.terr)*mu.pack
  }#k

  
  }
  ", fill=TRUE)
  sink()
  
################################################################################
  #  Bundle data, specify parameters & MCMC settings, and run JAGS
  
  # Define & Bundle Data
  Site <- length(unique(dynamic.ABwolf$cell_ID))
  win.data <- list("y"=ABwolf, "nSite"=Site, "nOcc"=9, "K"=3, "FC"=cov[,3,], 
                   "HARV"=cov[,5,], "HARV2"=cov[,6,], "area"=cov[,1,], 
                   "mu.terr"=1000, "mu.pack"=6.76)
  
  #  Initial Values
  #  Use naive occupancy estimate as initial value 
  
  #  Calculating naive occupancy based on raw encouter histories
  zst <- apply(ABwolf,c(1,3), max,na.rm=TRUE) # Take max value for each row (1) for each year (3)
  zst[zst=="-Inf"] <- NA 	# Apply function turns NAs into -Inf so have to change it back
  zst[zst==0] <- 0        # Probably not necessary- relic of false-positive encounter history
  zst[zst==1] <- 1 
  
  inits <- function(){list(z=zst)}
  
  #  Parameters to keep track of and report
  #  DON'T FORGET to ask for "psi" on the models that have covariates on them!
  params <- c("B0.psi","B1.pHvMed", "B1.pHvHi", "B0.p.rend", "B0.p.hunt", 
              "n.occ", "k.psi", "k.p.rnd", "k.p.hs", "n.packs", "mean.pack", 
              "n.wolves", "ann.occ", "mean.occ", "psi") # "B1.FC", 
  
  #  MCMC Settings 
  ni <- 300000
  nt <- 4
  nb <- 150000
  nc <- 3
  
  # Call JAGS 
  out <- jags(win.data, inits, params, "SS_Wolf_basep_nwolf_pHARV_final.txt", n.chains=nc, 
              n.thin=nt, n.iter=ni, n.burnin=nb, jags.module = c("glm","dic"))
  
  jag.sum <- out$BUGSoutput$summary
  #write.table(x=jag.sum, file="./Output/SS_Wolf_basep_nwolf_pHARV_final.txt", sep="\t")
  SS_Wolf_basep_nwolf_pHARV_final <- out
  #save(SS_Wolf_basep_nwolf_FC_pHARV_output, file ="./Output/SS_Wolf_basep_nwolf_pHARV_final.RData")

  print(out, dig=2)
  
  out2 <- out
  #mcmcplot(out2)
