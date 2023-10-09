# Install packages and library
#install.packages("tsintermittent")
#install.packages("xlsx")
#install.packages("lightgbm")
#install.packages("readxl")
#install.packages("boot")
#install.packages("expm")
#install.packages("igraph")
#install.packages("matlab")
#install.packages("markovchain")
#install.packages("Rcpp")
#install.packages("TraMineR")
#install.packages("dplyr")
#install.packages("RSNNS")
#install.packages("tidyr")
#install.packages("collapse")
#install.packages("runner")
#install.packages("ggplot2")
#install.packages("data.table")
#install.packages("greybox")
#install.packages("beepr")
#install.packages("knitr")
#install.packages("greybox")
#install.packages("sjmisc")
#install.packages("weights")
#install.packages("cNORM")
#install.packages("keras")
#install.packages("randomForest")
#install.packages("inventorize")
#install.packages("keras")
#install.packages("caret")
#install.packages("cvTools")
#install.packages("h2o")
#install.packages("doParallel")
library(xlsx)
library(tsintermittent)
library(lightgbm)
library(readxl)
library(boot)
library(expm)
library(igraph)
library(matlab)
library(markovchain)
library(Rcpp)
library(TraMineR)
library(dplyr)
library(RSNNS)
library(tidyr)
library(collapse)
library(runner)
library(ggplot2)
library(data.table)
library(greybox)
library(beepr)
library(knitr)
library(sjmisc)
library(weights)
library(cNORM)
library(keras)
library(tensorflow)
library(reticulate)
library(randomForest)
library(inventorize)
library(Rtools)
library(caret)
library(cvTools)
library(h2o)
library(tuneRanger)
library(ranger)
library(doParallel)
############################# PERFORMANCE MEASURES #############################

# Forecasting Accuracy Measures (FAM)

##GMAE function
GMAE_fun <- function(predicted, actual) {
  absolute_errors <- abs(predicted - actual)
  gmae <- exp(mean(log(absolute_errors), na.rm = TRUE))
  return(gmae)
}

##FAM
FAM <- function(train_data, test_data, predictions) {
  MSE_score <- NULL
  MASE_score <- NULL
  RMSSE_score <- NULL
  GMAE_score <- NULL
  for (i in 1:ncol(test_data)){
    MSE_score <- cbind(MSE_score,  MSE(t(test_data[i]), t(predictions[,i])))
    MASE_score <- cbind(MASE_score,  MASE(t(test_data[i]), t(predictions[,i]), mean(abs(t(train_data[i])))))
    RMSSE_score <- cbind(RMSSE_score,  RMSSE(t(test_data[i]), t(predictions[,i]), mean(abs(t(train_data[i])))))
    GMAE_score <- cbind(GMAE_score, GMAE_fun(t(predictions[,i]), t(test_data[i])))
  }
  MSE_score <- mean(MSE_score)
  MASE_score <- mean(MASE_score)
  RMSSE_score <- mean(RMSSE_score)
  GMAE_score <- mean(GMAE_score)
  
 return(c(MSE_score, MASE_score, RMSSE_score, GMAE_score))
}

# Inventory performance measures (IPM)
### ESC GAMMA ###

ESC_GAMMA <- function(k, alpha, R) {
  ESC = k/alpha - R - k/alpha * pgamma(alpha*R, k+1, 1) + R * pgamma(alpha*R, k, 1)
  return(ESC)
}

IPM_GAMMA <- function(data, train_data, test_data, predictions, prices, data_name, Method) {
  # Containing fill rate for each item for each TFR (col ~ item, row ~ TFR)
  # Each entry = total supply / total demand of an item for a TFR
  FillRates <- as.data.frame(matrix(ncol = ncol(test_data), nrow = nTFR))
  colnames(FillRates) <- colnames(data)
  rownames(FillRates) <- targetFillRates
  
  # Containing total supply for each item for each TFR (col ~ item, row ~ TFR)
  TotalSupply <- as.data.frame(matrix(ncol = ncol(test_data), nrow = nTFR))
  colnames(TotalSupply) <- colnames(data)
  rownames(TotalSupply) <- targetFillRates
  
  # Containing total demand for each item
  TotalDemand <- colSums(test_data)
  
  # Containing average holding cost for each item for each TFR (col ~ item, row ~ TFR)
  HoldingCosts <- as.data.frame(matrix(ncol = ncol(test_data), nrow = nTFR))
  colnames(HoldingCosts) <- colnames(data)
  rownames(HoldingCosts) <- targetFillRates
  
  t = nrow(train_data)
  h = nrow(test_data)
  
  for (i in 1:ncol(test_data)) {
    for (j in 1:nTFR) {
      r <- targetFillRates[j]
      inventory <- data.frame(TFR = rep(r, length = h), Sigma = NA, Mu = NA, K = NA, Lambda = NA,
                              R = NA, Q = NA, ILpreD = NA, D = NA, S = NA, ILpostD = NA, Cost = NA)
      inventory$D <- rep(test_data[,i])
      for (k in 1:h) {
        inventory$Sigma[k] <- sd(data[1:(t+k-1),i])
      }
      inventory$Mu <- predictions[,i]
      inventory$K      <- inventory$Mu^2 / inventory$Sigma^2   # shape
      inventory$Lambda <- inventory$Sigma^2 / inventory$Mu     # scale
      alpha  <- 1 / inventory$Lambda                           # rate
      
      for (g in 1:h) {
        if (inventory$Mu[g] <= 0) {
          inventory$R[g] <- 0
          next
        }
        maxR <- round(max(data[1:(t+g-1),i]))
        R <- 0
        lossTarget <- (1 - r) * inventory$Mu[g]
        tempQ <- inventory$Mu[g]
        shortagePerCycle <- ESC_GAMMA(inventory$K[g], alpha[g], R) #- ESC_GAMMA(inventory$K[g], alpha[g], R + tempQ)
        while (shortagePerCycle > lossTarget) {
          R <- R + 1
          if (R == maxR) {
            break
          }
          shortagePerCycle <- ESC_GAMMA(inventory$K[g], alpha[g], R) #- ESC_GAMMA(inventory$K[g], alpha[g], R + tempQ)
        }
        inventory$R[g] <- R
      }
      
      inventory$Q[1] <- 0                              # Initial order quantity
      inventory$ILpreD[1] <- inventory$R[1]
      inventory$S[1] <- min(inventory$ILpreD[1], inventory$D[1])
      inventory$ILpostD[1] <- inventory$ILpreD[1] - inventory$D[1]
      
      for (g in 2:h) {
        inventory$Q[g] <- max(0, inventory$R[g] - inventory$ILpostD[g-1])
        inventory$ILpreD[g] <- max(inventory$ILpostD[g-1], inventory$R[g])
        inventory$S[g] <- min(inventory$ILpreD[g], inventory$D[g])
        inventory$ILpostD[g] <- inventory$ILpreD[g] - inventory$D[g]
      }
      
      Cost <- 0.25 * inventory$ILpreD * prices[i]
      Cost[Cost < 0] <- 0
      inventory$Cost <- Cost
      
      if (TotalDemand[i] == 0) {
        FillRates[j,i] <- 0
      } else {
        FillRates[j,i] <- sum(inventory$S) / sum(inventory$D)
      }
      TotalSupply[j,i] <- sum(inventory$S)
      HoldingCosts[j,i] <- mean(inventory$Cost)
    } # end loop r
    print(paste("Done", i))
  } # end loop i
  
  totalAFR <- rowSums(TotalSupply) / sum(TotalDemand)
  avgAFR <- rowMeans(FillRates)
  
  # Holding Cost (for each target fill rate) is the sum of holding costs of all items
  holdingCost <- rowSums(HoldingCosts)
  
  # Service Level data frame
  ServiceLevel <- data.frame(avgAFR, totalAFR, holdingCost, targetFillRates, Method)
  colnames(ServiceLevel) <- c("AchievedFillRates_Avg", "AchievedFillRates_Total", "HoldingCosts", "TargetFillRates", "Method")
  
  # Saving
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\IPM\\GAMMA_IPM_", Method, "_", data_name, ".Rda", sep = "")
  save(ServiceLevel, FillRates, TotalDemand, TotalSupply, file = file_name)
  
  return(ServiceLevel)
}


#############################STATISTICAL METHODS################################
############################# CROSTON, SBA & Pennings DLP #########################

# Croston
Croston_Perf <- function(train_data, test_data, data_name) {
  h <- nrow(test_data)                               # Number of periods to be predicted (for each item)
  numItems <- ncol(train_data)                       # Number of items
  predictions <- matrix(nrow = h, ncol = numItems)   # Matrix storing the predictions
  
  # For each item
  for (j in 1:numItems) {
    x <- rbind(train_data[j], test_data[j])          # All demands (train + test) of item j
    t <- nrow(train_data)                            # Number of periods used for predicting (updated every iteration below)
    
    for (i in 1:h) {
      predictions[i,j] <- crost(x[1:t,], h = 1, w = NULL, nop = 2, type = "croston", cost = "mar", init = "naive", init.op = FALSE, na.rm = TRUE)$frc.out
      t = t + 1
    }
    
    print(paste("Done item", j))
  }
  
  # Store predictions as data frame
  predictions <- as.data.frame(predictions)
  
  # Saving
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\Stat_Methods\\Croston_", data_name, ".Rda", sep = "")
  save(predictions, file = file_name)
  return(predictions)
}


# SBA
SBA_Perf <- function(train_data, test_data, data_name) {
  h <- nrow(test_data)                               # Number of periods to be predicted (for each item)
  numItems <- ncol(train_data)                       # Number of items
  predictions <- matrix(nrow = h, ncol = numItems)   # Matrix storing the predictions
  
  # For each item
  for (j in 1:numItems) {
    x <- rbind(train_data[j], test_data[j])          # All demands (train + test) of item j
    t <- nrow(train_data)                            # Number of periods used for predicting (updated every iteration below)
    
    for (i in 1:h) {
      predictions[i,j] <- crost(x[1:t,], h=1, w=NULL, nop=2, type="sba", cost="mar", init="naive", init.opt=FALSE, na.rm=TRUE)$frc.out
      t = t + 1
    }
    
    print(paste("Done item", j))
  }
  
  # Store predictions as data frame
  predictions <- as.data.frame(predictions)
  
  # Saving
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\Stat_Methods\\SBA", data_name, ".Rda", sep = "")
  save(predictions, file = file_name)
  return(predictions)
}

#DLP 
#Code from Van Dalen 
# Evaluate DLP for h-period ahead forecast, h = nrow(test_data)
DLP_Perf <- function(train_data, test_data, data_name) {
  h <- nrow(test_data)                               # Number of periods to be predicted (for each item)
  numItems <- ncol(train_data)                       # Number of items
  predictions <- matrix(nrow = h, ncol = numItems)   # Matrix storing the predictions
  
  # For each item
  for (j in 1:numItems) {
    x <- rbind(train_data[j], test_data[j])          # All demands (train + test) of item j
    t <- nrow(train_data)                            # Number of periods used for predicting (updated every iteration below)
    
    # Initialize variables for DLP
    fc.y <- rep(NA, t + h)                           # Forecast demand
    fc.tau <- rep(NA, t + h)                         # Forecast lead-time
    since.Last <- rep(NA, t + h)                     # Periods since last demand
    p.t <- 1                                        # Initial value of p.t (intermittent factor)
    
    for (i in 1:h) {
      # Croston's method to forecast demand and lead-time
      croston_result <- crost(x[1:t,], h = 1, w = NULL, nop = 2, type = "croston", cost = "mar", init = "naive", init.opt = FALSE, na.rm = TRUE)
      
      # Update DLP variables
      fc.y[t+i] <- croston_result$frc.out
      fc.tau[t+i] <- t
      since.Last[t+i] <- 1
      
      # Calculate DLP forecast
      p <- 1 / fc.tau[t+i]
      predictions[i,j] <- (fc.y[t+i] / fc.tau[t+i]) * (h + (since.Last[t+i] - (1 - p) / p) * (1 - (1 - p) ^ h))
      
      t = t + 1
    }
    
    print(paste("Done item", j))
  }
  
  # Store predictions as data frame
  predictions <- as.data.frame(predictions)
  
  # Saving
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\DLP.R", data_name, ".Rda", sep = "")
  save(predictions, file = file_name)
  
  return(predictions)
}

DLP_Perf(trainSIM1, testSIM1, "SIM1")

############################# RUNNING METHODS ON DATA SETS #############################


# Set seed for reproducibility.
set.seed(8377)


########## CROSTON ##########
ptm <- proc.time()

##### SIM1 #####
Croston_SIM1 <- Croston_Perf(trainSIM1,testSIM1, "SIM1")
Croston_FAM_SIM1 <- FAM(trainSIM1, testSIM1, Croston_SIM1)
Croston_IPM_GAMMA_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1, Croston_SIM1, pricesSIM1, "SIM1" , "Croston")
Croston_FAM_SIM1
Croston_IPM_GAMMA_SIM1

##### SIM2 #####
Croston_SIM2 <- Croston_Perf(trainSIM2, testSIM2, "SIM2")
Croston_FAM_SIM2 <- FAM(trainSIM2, testSIM2, Croston_SIM2)
Croston_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2, Croston_SIM2, pricesSIM2, "SIM2", "Croston")
Croston_FAM_SIM2
Croston_IPM_SIM2
##### SIM3 #####
Croston_SIM3 <- Croston_Perf(trainSIM3, testSIM3, "SIM3")
Croston_FAM_SIM3 <- FAM(trainSIM3, testSIM3, Croston_SIM3)
Croston_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3, Croston_SIM3, pricesSIM3, "SIM3", "Croston")
Croston_FAM_SIM3
##### SIM4 #####
Croston_SIM4 <- Croston_Perf(trainSIM4, testSIM4, "SIM4")
Croston_FAM_SIM4 <- FAM(trainSIM4, testSIM4, Croston_SIM4)
Croston_IPM_SIM4 <- IPM_GAMMA(SIM4 ,trainSIM4, testSIM4, Croston_SIM4, pricesSIM4, "SIM4", "Croston")
Croston_FAM_SIM4

##### MAN #####
Croston_MAN <- Croston_Perf(trainMAN, testMAN, "MAN")
Croston_FAM_MAN <- FAM(trainMAN, testMAN, Croston_MAN)
Croston_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN, Croston_MAN, pricesMAN, "MAN","Croston")
Croston_FAM_MAN
##### BRAF #####
Croston_BRAF <- Croston_Perf(trainBRAF, testBRAF, "BRAF")
Croston_FAM_BRAF <- FAM(trainBRAF, testBRAF, Croston_BRAF)
Croston_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF, Croston_BRAF, pricesBRAF, "BRAF", "Croston")
Croston_FAM_BRAF
##### AUTO #####
Croston_AUTO <- Croston_Perf(trainAUTO, testAUTO, "AUTO")
Croston_FAM_AUTO <- FAM(trainAUTO, testAUTO, Croston_AUTO)
Croston_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO, Croston_AUTO, pricesAUTO, "AUTO", "Croston")
Croston_IPM_AUTO
#### OIL ####
Croston_OIL <- Croston_Perf(trainOIL, testOIL, "OIL")
Croston_FAM_OIL <- FAM(trainOIL, testOIL, Croston_OIL)
Croston_IPM_OIL <- IPM_GAMMA(OIL, trainOIL, testOIL, Croston_OIL, pricesOIL, "OIL", "Croston")



##### Print running time #####
print("time Croston total")
proc.time() - ptm




########## SBA ##########
ptm <- proc.time()

##### SIM1 #####
SBA_SIM1 <- SBA_Perf(trainSIM1, testSIM1, "SIM1")
SBA_FAM_SIM1 <- FAM(trainSIM1, testSIM1, SBA_SIM1)
SBA_IPM_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1, SBA_SIM1, pricesSIM1, "SIM1", "SBA")
SBA_FAM_SIM1
##### SIM2 #####
SBA_SIM2 <- SBA_Perf(trainSIM2, testSIM2, "SIM2")
SBA_FAM_SIM2 <- FAM(trainSIM2, testSIM2, SBA_SIM2)
SBA_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2, SBA_SIM2, pricesSIM2, "SIM2","SBA")
SBA_FAM_SIM2
##### SIM3 #####
SBA_SIM3 <- SBA_Perf(trainSIM3, testSIM3, "SIM3")
SBA_FAM_SIM3 <- FAM(trainSIM3, testSIM3, SBA_SIM3)
SBA_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3, SBA_SIM3, pricesSIM3, "SIM3", "SBA")
SBA_FAM_SIM3
##### SIM4 #####
SBA_SIM4 <- SBA_Perf(trainSIM4, testSIM4, "SIM4")
SBA_FAM_SIM4 <- FAM(trainSIM4, testSIM4, SBA_SIM4)
SBA_IPM_SIM4 <- IPM_GAMMA(SIM4, trainSIM4, testSIM4, SBA_SIM4, pricesSIM4, "SIM4", "SBA")
SBA_FAM_SIM4
##### MAN #####
SBA_MAN <- SBA_Perf(trainMAN, testMAN, "MAN")
SBA_FAM_MAN <- FAM(trainMAN, testMAN, SBA_MAN)
SBA_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN, SBA_MAN, pricesMAN, "MAN", "SBA")
SBA_FAM_MAN
SBA_IPM_MAN
##### BRAF #####
SBA_BRAF <- SBA_Perf(trainBRAF, testBRAF, "BRAF")
SBA_FAM_BRAF <- FAM(trainBRAF, testBRAF, SBA_BRAF)
SBA_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF, SBA_BRAF, pricesBRAF, "BRAF", "SBA")
SBA_FAM_BRAF
##### AUTO #####
SBA_AUTO <- SBA_Perf(trainAUTO, testAUTO, "AUTO")
SBA_FAM_AUTO <- FAM(trainAUTO, testAUTO, SBA_AUTO)
SBA_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO, SBA_AUTO, pricesAUTO, "AUTO", "SBA")
SBA_FAM_AUTO
####OIL####
SBA_OIL <- SBA_Perf(trainOIL, testOIL, "OIL")
SBA_FAM_OIL <- FAM(trainOIL, testOIL, SBA_OIL)
SBA_IPM_OIL <- IPM_GAMMA(OIL ,trainOIL ,testOIL, SBA_OIL, pricesOIL, "OIL", "SBA")
SBA_FAM_OIL
##### Print running time #####
print("time SBA total")
proc.time() - ptm

####### DLP ####### 
#Code from Van Dalen 
# Evaluate DLP for h-period ahead forecast, h = nrow(test_data)
# DLP 

DLP_Perf <- function(train_data, test_data, data_name) {
  h <- nrow(test_data)                               # Number of periods to be predicted (for each item)
  numItems <- ncol(train_data)                       # Number of items
  predictions <- matrix(nrow = h, ncol = numItems)   # Matrix storing the predictions
  
  # For each item
  for (j in 1:numItems) {
    x <- rbind(train_data[j], test_data[j])          # All demands (train + test) of item j
    t <- nrow(train_data)                            # Number of periods used for predicting (updated every iteration below)
    
    # Initialize variables for DLP
    fc.y <- rep(NA, t)         # Forecast demand
    fc.tau <- rep(NA, t)       # Forecast lead-time
    since.Last <- rep(NA, t)   # Periods since the last demand
    p.t <- 1                   # Initial value of p.t (intermittent factor)
    
    for (i in 1:h) {
      # Croston's method to build DLP on it
      croston_result <- crost(x[1:t,], h = 1, w = NULL, nop = 2, type = "croston", cost = "mar", init = "naive", init.opt = FALSE, na.rm = TRUE)
      
      # Update DLP variables
      fc.y[t+i] <- croston_result$frc.out
      fc.tau[t+i] <- t
      since.Last[t+i] <- 1
      
      # Calculate DLP forecast
      p.t <- 1 / fc.tau[t+i]
      
      predictions[i,j] <- (fc.y[t+i] / fc.tau[t+i]) * (h + (since.Last[t+i] - (1 - p.t) / p.t) * (1 - (1 - p.t) ^ h))
      
      t = t + 1
    }
    
    print(paste("Done item", j))
  }
  
  # Store predictions as data frame
  predictions <- as.data.frame(predictions)
  
  # Saving
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\Stat_Methods\\DLP_", data_name, ".Rda", sep = "")
  save(predictions, file = file_name)
  
  return(predictions)
}

#SIM1
DLP_SIM1 <- DLP_Perf(trainSIM1, testSIM1, "SIM1")
DLP_FAM_SIM1 <- FAM(trainSIM1, testSIM1, DLP_SIM1)
DLP_IPM_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1, DLP_SIM1, pricesSIM1, "SIM1", "DLP")
DLP_FAM_SIM1

#SIM2
DLP_SIM2 <- DLP_Perf(trainSIM2, testSIM2, "SIM2")
DLP_FAM_SIM2 <- FAM(trainSIM2, testSIM2, DLP_SIM2)
DLP_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2, DLP_SIM2, pricesSIM2, "SIM2", "DLP")
DLP_FAM_SIM2

#SIM3
DLP_SIM3 <- DLP_Perf(trainSIM3, testSIM3, "SIM3")
DLP_FAM_SIM3 <- FAM(trainSIM3, testSIM3, DLP_SIM3)
DLP_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3, DLP_SIM3, pricesSIM3, "SIM3", "DLP")
DLP_FAM_SIM3

#SIM4
DLP_SIM4 <- DLP_Perf(trainSIM4, testSIM4, "SIM4")
DLP_FAM_SIM4 <- FAM(trainSIM4, testSIM4, DLP_SIM4)
DLP_IPM_SIM4 <- IPM_GAMMA(SIM4, trainSIM4, testSIM4, DLP_SIM4, pricesSIM4, "SIM4", "DLP")
DLP_FAM_SIM4

#MAN
DLP_MAN <- DLP_Perf(trainMAN, testMAN, "MAN")
DLP_FAM_MAN <- FAM(trainMAN, testMAN, DLP_MAN)
DLP_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN, DLP_MAN, pricesMAN, "MAN", "DLP")
DLP_FAM_MAN

#BRAF
DLP_BRAF <- DLP_Perf(trainBRAF, testBRAF, "BRAF")
DLP_FAM_BRAF <- FAM(trainBRAF, testBRAF, DLP_BRAF)
DLP_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF, DLP_BRAF, pricesBRAF, "BRAF", "DLP")
DLP_FAM_BRAF

#AUTO
DLP_AUTO <- DLP_Perf(trainAUTO, testAUTO, "AUTO")
DLP_FAM_AUTO <- FAM(trainAUTO, testAUTO, DLP_AUTO)
DLP_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO, DLP_AUTO, pricesAUTO, "AUTO", "DLP")
DLP_FAM_AUTO

#OIL
DLP_OIL <- DLP_Perf(trainOIL, testOIL, "OIL")
DLP_FAM_OIL <- FAM(trainOIL, testOIL, DLP_OIL)
DLP_IPM_OIL <- IPM_GAMMA(OIL, trainOIL, testOIL, DLP_OIL, pricesOIL, "OIL", "DLP")
DLP_FAM_OIL
