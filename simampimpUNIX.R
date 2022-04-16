##### SIMULATION #####

## Looped conditions
# Number of individuals: ns = 50, 100, 200
# Number of repeated measures: nis = 10, 30

## Manual conditions
# a0 = .2 
# v2 = 1.1
# c21 = 0.3
# Proportion missing (%): p = 10

## preliminaries
#install.packages("mice")
#install.packages("data.table")
#install.packages("miceadds")
#install.packages("dplyr")
#install.packages("modeest")

library("mice")
library("data.table")
library("miceadds")
library("dplyr")
library("statip")

## fixed population values (from values in Table 4)
beta0 = 0 # intercept
beta1 = .2 # X1 coefficient 
beta2 = .1 # X2.b coefficient
gamma1 = .1 # W1 coefficient 
gamma2 = .1 # W2.b coefficient
a1 = .2 # W1 location coefficient
a2 = .1 # W2.b location coefficient
ta0 = .2 # WS variance (WS intercept)
ta1 = .1 # X1 scale coefficient 
ta2 = .2 # X2.b scale coefficientd
ta3 = .2 # W1 scale coefficient
ta4 = .1 # W2.b scale coefficient

# varying population values
a0 = .2 # location random effect variance (BS intercept) 
v2 = 1.1 # scale random effect variance (BS intercept)
c21 = 0.3 # correlation between the random effects

## simulation
rm_levels = c(10, 30)
repeatedmeasures = list()
for (rm in rm_levels){
  
  subject_levels = c(50, 200, 500)
  subjects = list()
  for (sl in subject_levels){
    
    ns =  sl # number of subjects
    nis = rm # number of rm per subject
    
    reps = 1000
    replications = list()
    for(k in 1:reps){
      
      y <- array(NA, dim = c(ns*nis))
      X1 <- array(NA, dim = c(ns*nis)) # LEVEL 1 continuous covariate corresponding to beta1 and ta1
      X2.b <- array(NA, dim = c(ns*nis)) # LEVEL 1 binary covariate corresponding to beta2 and ta2
      W1 <- array(NA, dim = c(ns*nis)) # LEVEL 2 continuous covariate corresponding to gamma1 and ta3 and a1
      W2.b <- array(NA, dim = c(ns*nis)) # LEVEL 2 binary covariate corresponding to gamma2 and ta4 and a2
      X1.c <- array(NA, dim = c(ns*nis)) # centered version of X1
      
      subj <- array(NA, dim = c(ns*nis))
      y.ni <- array(NA, dim = c(nis))
      x1.ni <- array(NA, dim = c(nis))
      x2.b.ni <- array(NA, dim = c(nis))
      sub.ni <- array(NA, dim = c(nis))
      
      for(n in 1:ns){ # for each individual:				
        sub1 = rnorm(1) # single random value from a standard normal distribution	(standardized random location effect)					
        sub2 = rnorm(1) # single random value from a standard normal distribution	(standardized random scale effect)
        
        # LEVEL 2
        w1 = rnorm(1) # single random value from a standard normal distribution (standardized level 2 continuous covariate)	
        w2.b = rbinom(1, size = 1, prob = .1) # single random value from a binomial distribution, p = .1 (standardized level 2 binary covariate)
        
        # form Cholesky elements
        s1  = sqrt(exp(a0 + w1*a1 + w2.b*a2)) # BS location random effect sd as a function of covariates (top left Cholesky element) (Hedeker et al., 2008)
        s21 = c21 / s1 # correlation between s1 and s2 (bottom left Cholesky element) (Hedeker et al., 2008)							
        s2  = sqrt(v2 - (c21*c21 / s1*s1)) # BS scale random effect sd (bottom right Cholesky element) (Hedeker et al., 2008)
        
        # LEVEL 1
        for(ni in 1:nis){	# for each repeated measure:						
          err =  rnorm(1)							
          x1 = rnorm(1) # single random value from a standard normal distribution (standardized level 1 continuous covariate)	
          x2.b = rbinom(1, size = 1, prob = .1) # single random value from a binomial distribution, p = .1 (standardized level 1 binary covariate)
          g   = sqrt(exp(ta0 + x1*ta1 + x2.b*ta2 + w1*ta3 + w2.b*ta4 + s21*sub1 + s2*sub2)) # WS scale random effect sd as a function of covariates (Hedeker et al., 2008)
          y.ni[ni]  =  beta0 + x1*beta1 + x2.b*beta2 + w1*gamma1 + w2.b*gamma2 + s1*sub1 + g*err # mean of y as a function of covariates, the BS location random effect sd, and the WS scale random effect sd (Hedeker et al., 2008)
          x1.ni[ni] <-  x1
          x2.b.ni[ni] <- x2.b
          sub.ni[ni] <- n
        }
        y[((n-1)*nis+1):(n*nis)] <- y.ni
        x1.ni.mean <- mean(x1.ni) # for each individual, calculate the mean of the repeated measures (calculate the within individual mean)
        X1[((n-1)*nis+1):(n*nis)] <- x1.ni
        X1.c[((n-1)*nis+1):(n*nis)] <- x1.ni - x1.ni.mean # X1.c is each repeated measure minus the within individual mean (center X1 by the within individual mean)
        W1[((n-1)*nis+1):(n*nis)] <- rep(w1, nis) 
        X2.b[((n-1)*nis+1):(n*nis)] <- x2.b.ni
        W2.b[((n-1)*nis+1):(n*nis)] <- rep(w2.b, nis)
        subj[((n-1)*nis+1):(n*nis)] <- sub.ni
      }
      W1_mean <- mean(W1) # calculate the grand mean of W1
      W1.c <- W1 - W1_mean # center W1 by the grand mean
      time = rep(c(0:(nis-1)), ns)
      data <- data.frame(y, X1, W1, X2.b, W2.b, subj, time)  
      
      ## amputation
      
      # ampute the time-invariant variables W2.b and W1.c (missingness is a function of the time-varying variables y, X2.b, and X1.c)
      # create person-specific averages for time-varying variables
      y.avg <- aggregate(data$y, list(data$subj), FUN = mean)
      names(y.avg) <- c("subj", "y.avg")
      X2.b.mode <- aggregate(data$X2.b, list(data$subj), FUN = mfv1)
      names(X2.b.mode) <- c("subj", "X2.b.mode")
      X1.avg <- aggregate(data$X1, list(data$subj), FUN = mean)
      names(X1.avg) <- c("subj", "X1.avg")
      merge1 <- merge(y.avg, X2.b.mode, by = "subj")
      merge2 <- merge(merge1, X1.avg, by = "subj") 
      
      # order by subj for merging purposes
      merge2$subj <- as.numeric(merge2$subj)
      merge2 <- merge2[order(merge2$subj), ]
      data$subj <- as.numeric(data$subj)
      data <- data[order(data$subj), ]
      
      # merge person-specific averages (merge2) with data
      data_avg <- merge(merge2, data, by = "subj")
      
      # transpose long to wide (spreading y, X3.b, X1.c, and time)
      setDT(data_avg)
      widedata_avg <- dcast(data_avg, subj + y.avg + X2.b.mode + X1.avg + W2.b + W1 ~ time, value.var = c("y", "X2.b", "X1"))
      
      # create patterns (ampute W2.b and W1) and weights (for MAR, weights of variables that will be made incomplete should be zero)
      mypatterns <- rbind(c(1, 1, 1, 1, 0, 0, rep(1, nis), rep(1, nis), rep(1, nis)),
                          c(1, 1, 1, 1, 1, 0, rep(1, nis), rep(1, nis), rep(1, nis)),
                          c(1, 1, 1, 1, 0, 1, rep(1, nis), rep(1, nis), rep(1, nis))) # 0: missing, 1: complete
      # I don't want subject, y_0:y_9, X2.b_0:X2.b_9, or X1_0:X1_9 to have an influence on the probability of missingness
      myweights <-  rbind(c(0, 1, 1, 1, 0, 0, rep(0, nis), rep(0, nis), rep(0, nis)), 
                          c(0, 1, 1, 1, 1, 0, rep(0, nis), rep(0, nis), rep(0, nis)), 
                          c(0, 1, 1, 1, 0, 1, rep(0, nis), rep(0, nis), rep(0, nis))) # 0: missing, 1: complete
      
      # ampute data
      amputed_data <- ampute(widedata_avg, prop = .1, patterns = mypatterns, mech = "MAR", weights = myweights, type = "RIGHT")
      # mids objects as data frames
      widedata_amputed <- amputed_data$amp
      
      # drop the person-specific averages for time-varying variables
      drop = c("y.avg", "X2.b.mode", "X1.avg")
      widedata_amputed = widedata_amputed[,!(names(widedata_amputed) %in% drop)]
      
      # transpose wide to long
      longdata_amputed <- reshape (widedata_amputed, direction = 'long', 
                                   varying=list(grep("y", colnames(widedata_amputed), value=T), 
                                                grep("X2.b", colnames(widedata_amputed), value=T), 
                                                grep("X1", colnames(widedata_amputed), value=T)),
                                   timevar = "time", 
                                   times = as.character(c(0:(nis-1))),
                                   v.names = c("y", "X2.b", "X1"),
                                   idvar = "subj")
      
      # order by subj
      longdata_amputed$subj <- as.numeric(longdata_amputed$subj)
      longdata_amputed <- longdata_amputed[order(longdata_amputed$subj), ]
      
      ## imputation 
      
      # settings
      imp0 <- mice(longdata_amputed, maxit = 0) # setup run of mice()
      impmethod <- imp0$method # default method
      impmethod[c("W2.b", "W1")] <- c("2lonly.pmm", "2lonly.norm") # 2lonly.norm imputes level two continuous variables allowing for residual variance to vary across individuals (2lonly.pan assumes homogeneous residual variance)
      
      #pm <- imp0$predictorMatrix # default predictor matrix
      pm <- matrix(1, ncol(longdata_amputed), ncol(longdata_amputed))
      rownames(pm) <- colnames(pm) <- colnames(longdata_amputed)
      pm[, "subj"] <- -2 # identify cluster variable
      pm[, "time"] <- 0
      pm["W2.b", c("subj", "W2.b", "W1", "time", "y", "X2.b", "X1")] <- c(-2, 0, 1, 0, 1, 1, 1)
      pm["W1",   c("subj", "W2.b", "W1", "time", "y", "X2.b", "X1")] <- c(-2, 1, 0, 0, 1, 1, 1)
      
      # run multiple imputations
      imp <- mice(longdata_amputed, m=20, predictorMatrix=pm, 
                  method=impmethod, maxit=10)
      
      # store simulated amputed data in a list called replications
      replications[[k]] <- complete(imp, action = "long")
    }
    subjects[[sl]] <- replications
  }
  repeatedmeasures[[rm]] <- subjects
}

# store each combination of conditions of sample size and number of repeated measures inside of a data frame
data10_50 <- bind_rows(repeatedmeasures[[10]][[50]], .id = "replication")
 
data10_200 <- bind_rows(repeatedmeasures[[10]][[200]], .id = "replication")
 
data10_500 <- bind_rows(repeatedmeasures[[10]][[500]], .id = "replication")
 
data30_50 <- bind_rows(repeatedmeasures[[30]][[50]], .id = "replication")

data30_200 <- bind_rows(repeatedmeasures[[30]][[200]], .id = "replication")
 
data30_500 <- bind_rows(repeatedmeasures[[30]][[500]], .id = "replication")

# output dataframes to .txt files
setwd("/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data10_50")
write.table(x = data10_50, file='data10_50.txt', row.names = F, col.names = T)

setwd("/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data10_200")
write.table(x = data10_200, file='data10_200.txt', row.names = F, col.names = T)

setwd("/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data10_500")
write.table(x = data10_500, file='data10_500.txt', row.names = F, col.names = T)

setwd("/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data30_50")
write.table(x = data30_50, file='data30_50.txt', row.names = F, col.names = T)

setwd("/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data30_200")
write.table(x = data30_200, file='data30_200.txt', row.names = F, col.names = T)

setwd("/scratch/mcraft/LSMsimulation/p_10/sim0.2_1.1_0.3/data30_500")
write.table(x = data30_500, file='data30_500.txt', row.names = F, col.names = T)

# References:
# Hedeker, D., Mermelstein, R. J., & Demirtas, H. (2008). 
# An application of a mixed-effects location scale model for analysis of ecological momentary assessment (EMA) data. 
# Biometrics, 64, 627–634.
