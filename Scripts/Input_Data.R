    #  Input data for wolf occ model
    #  May 2016
    #  Updated with correctly scaled covariates January 2017
    
    #  Reads in encounter history data and covariate data
    #  Creates arrays for EH (ABwolf), covariates on occ. parameters (cov), 
    #  and covariates on detection probabilities (survey)
################################################################################
    #  Load packages
    
    library(dplyr)
    library(tidyr)
################################################################################
    #  Set up encounter history array
    
    #  Read in data
    load("C:/Sarah B/Thesis POM/Final_Models_4_Thesis/Input/adj_POM_EH.RData")
    dynamic.ABwolf <- adj_POM_EH 
    
    #  adj_POM_EH includes false-positive detections- remove and correct EH for JAGS
    dynamic.ABwolf[,3:11] <- lapply(dynamic.ABwolf[,3:11], function(x){replace(x, x == 1, 0)})
    dynamic.ABwolf[,3:11] <- lapply(dynamic.ABwolf[,3:11], function(x){replace(x, x == 2, 1)})
    dynamic.ABwolf[,3:11] <- lapply(dynamic.ABwolf[,3:11], function(x){replace(x, x == 3, 1)})
    #  str(dynamic.ABwolf)
    
    
    #  Create 3D array from 2D data (dim = c(52,9,3) when not adjusted for small grid cells)
    ABwolf <- array(NA, dim=c(50,9,3)) # place holder for 41 sites, 9 occasions, 3 years
    
    for(i in 1:3){                           # making site x occasion matrix by year
      sel.rows <- dynamic.ABwolf$year==i
      ABwolf[,,i] <- as.matrix(dynamic.ABwolf)[sel.rows,3:11]
    }#i
    
    
################################################################################
    # Set up covariate arrays
    
    #  Read in data
    load("C:/Sarah B/Thesis POM/Final_Models_4_Thesis/Input/adj_POM_zCovs2.RData")
    dyn.cov <- adj_zCovs2
    #  str(dyn.cov)
    
    #  Create covariate array for occupancy covariates
    cov <- array(NA, dim=c(50,7,3), dimnames = list((unique(as.character(dyn.cov$cell_ID))), 
                                                    (colnames(dyn.cov[,3:9])), 
                                                    (unique(as.character(dyn.cov$year)))))
    
    #  Fill array with occupancy covariates
    for(k in 1:3){
      sel.rows <- dyn.cov$year==k
      cov[,,k] <- (as.numeric(as.matrix(dyn.cov)[sel.rows,3:9]))   # note:  need to change STOCK covariate to as.character for model to recognize it as categorical (or do I- there is a rank to them- L, M, H)
    }#i
    
    #  Create covariate array for detection probability covariates
    survey <- array(NA, dim=c(50,9,3), dimnames=list((unique(as.character(dyn.cov$cell_ID))), 
                                                     (colnames(dyn.cov[,10:18])), 
                                                     (unique(as.character(dyn.cov$year)))))
    #  Fill array with detection probability covariates
    for(k in 1:3){
      sel.rows <- dyn.cov$year==k
      survey[,,k] <- (as.numeric(as.matrix(dyn.cov)[sel.rows,10:18]))
    }