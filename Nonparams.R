############################NON-PARAMETRIC METHODS#########################################################

################## Willemain ###############
targetFillRates <- c(75:99)
targetFillRates <- targetFillRates / 100
nTFR <- length(targetFillRates)


dataName <- "MAN"
data <- MAN
train_data <- trainMAN
test_data <- testMAN
prices <- pricesMAN
sample <- sampleMAN

dataName <- "AUTO"
data <- AUTO
train_data <- trainAUTO
test_data <- testAUTO
prices <- pricesAUTO
sample <- sampleAUTO

dataName <- "BRAF"
data <- BRAF
train_data <- trainBRAF
test_data <- testBRAF
prices <- pricesBRAF
sample <- sampleBRAF

dataName <- "OIL"
data <- OIL
train_data <- trainOIL
test_data <- testOIL
prices <- pricesOIL
sample <- sampleOIL

dataName <- "SIM1"
data <- SIM1
train_data <- trainSIM1
test_data <- testSIM1
prices <- pricesSIM1
sample <- sampleSIM1

dataName <- "SIM2"
data <- SIM2
train_data <- trainSIM2
test_data <- testSIM2
prices <- pricesSIM2
sample <- sampleSIM2

dataName <- "SIM3"
data <- SIM3
train_data <- trainSIM3
test_data <- testSIM3
prices <- pricesSIM3
sample <- sampleSIM3

dataName <- "SIM4"
data <- SIM4
train_data <- trainSIM4
test_data <- testSIM4
prices <- pricesSIM4
sample <- sampleSIM4

####################################

ESC_WILL <- function(probLTD, mu, R) {
  if (R >= max(probLTD$Val)) {
    return(-1)
  } else {
    iStart <- which(probLTD$Val > R)[1]
    iEnd <- nrow(probLTD)
    ESC <- sum((probLTD$Val[iStart:iEnd] - R) * probLTD$Prob[iStart:iEnd])
    return(ESC)
  }
}

Method <- "Willemain"

####################################

seqFull <- data
seqFull[data != 0] <- 1
seqFull <- seqdef(seqFull)


predictions <- matrix(nrow = nrow(test_data), ncol = ncol(test_data))   # One-step ahead forecasts
# predictions2 <- matrix(nrow = nrow(test_data), ncol = ncol(test_data))   # Two-step ahead forecasts
h <- nrow(test_data)
t <- nrow(train_data) 

listLTD <- list()
FillRates <- as.data.frame(matrix(ncol = ncol(test_data), nrow = nTFR))
colnames(FillRates) <- colnames(data)
rownames(FillRates) <- targetFillRates

TotalSupply <- as.data.frame(matrix(ncol = ncol(test_data), nrow = nTFR))
colnames(TotalSupply) <- colnames(data)
rownames(TotalSupply) <- targetFillRates

TotalDemand <- colSums(test_data)

HoldingCosts <- as.data.frame(matrix(ncol = ncol(test_data), nrow = nTFR))
colnames(HoldingCosts) <- colnames(data)
rownames(HoldingCosts) <- targetFillRates

set.seed(8377)

b <- 1000       # Number of bootstrapping iterations
for (i in 1:ncol(seqFull)) { 
  leadtime <- 1
  LTD <- matrix(NA, nrow = b, ncol = h)
  for (k in 1:h) {
    # Define new train set / train sequence
    train <- data[1:(t+k-1),i]
    trainSeq <- seqFull[1:(t+k-1),i]
    nNonZeros <- length(which(trainSeq == 1))
    
    statesNames = c(0,1)
    # STEP 1: Estimate the transition probabilities for the two-state Markov model from historical demands (zero or non-zero)
    transprob <- createSequenceMatrix(trainSeq, toRowProbs = TRUE, sanitize = FALSE)
    
    # # The following piece of code only hold for AUTO, for the other type we remove
    # # Check type 2 (only AUTO & SIM3 has this type), and change probabilities of row 1 to the probability 0/1 appears in the train data
    # if (sum(transprob[1,]) == 0) {
    #   transprob[1,2] <- nNonZeros / length(trainSeq)
    #   transprob[1,1] <- 1 - transprob[1,2]
    # } else if (i == 874 || i == 1196 || i == 2106) { # Check type 6 (specific items), if yes then do the same as for type 2
    #   transprob[1,2] <- nNonZeros / length(trainSeq)
    #   transprob[1,1] <- 1 - transprob[1,2]
    # }
    # 
    # # The following piece of code only hold for SIM3, for the other type we remove
    # # Check type 2 (only AUTO & SIM3 has this type), and change probabilities of row 1 to the probability 0/1 appears in the train data
    # if (sum(transprob[1,]) == 0) {
    #   transprob[1,2] <- nNonZeros / length(trainSeq)
    #   transprob[1,1] <- 1 - transprob[1,2]
    # } else if (i == 4055) { # Check type 6 (specific items), if yes then do the same as for type 2
    #   transprob[1,2] <- nNonZeros / length(trainSeq)
    #   transprob[1,1] <- 1 - transprob[1,2]
    # }
    
    # Repeat steps 2-5 for 1000 times
    if (length(transprob) != 1) {
      mc <- new("markovchain", transitionMatrix = matrix((transprob), nrow = 2, dimnames = list(statesNames,statesNames)))
      for (y in 1:b) {
        # STEP 2: Use the Markov model (mc) to generate zero and nonzero values over the specified forecast horizon
        # In our case (one-step-ahead forecast), we just generate one prediction
        mcprediction <- as.numeric(markovchainSequence(n=1, markovchain=mc,t0 = tail(trainSeq, n=1),include.t0 = FALSE))
        
        if (mcprediction != 0) {
          # STEP 3: Replace the non-zero demands with a random demand value, with replacement, from the already known set of non-zero demands
          nonzeroes <- train[train != 0]
          mcprediction <- sample(nonzeroes, size = 1)
          # STEP 4: Jitter the values of the non-zero demands - a value close to the randomly selected demand value is selected
          # Formula: Jittered = 1 + int{X + Z * sqrt(X)}, where X is the randomly selected demand, Z is a standard normal random deviate
          mcprediction <- 1 + round(mcprediction + rnorm(1)*sqrt(mcprediction))
        }
        
        # STEP 5: Sum the predicted values over the forecast horizon, which is meaningless as we only have one forecast...
        LTD[y,k] <- mcprediction
      }
      predictions[k, i] <- mean(LTD[,k])
    } else { # When there is no zero-demand
      statesNames = c(1)
      mc <- new("markovchain", transitionMatrix=matrix((transprob), byrow=TRUE, nrow=1, dimnames=list(statesNames,statesNames)))
      for (y in 1:b) {
        mcprediction <- as.numeric(markovchainSequence(n=1, markovchain=mc,t0 = tail(trainSeq, n=1),include.t0 = FALSE))
        if (mcprediction != 0) {
          nonzeroes <- train[train != 0]
          mcprediction <- sample(nonzeroes, size = 1)
          mcprediction <- 1 + round(mcprediction + rnorm(1)*sqrt(mcprediction))
        }
        
        LTD[y,k] <- mcprediction
      }
      predictions[k, i] <- mean(LTD[,k])
    }
  }
  
  # Inventory
  for (j in 1:nTFR) {
    r <- targetFillRates[j]
    inventory <- data.frame(TFR = rep(r, length = h), Sigma = NA, Mu = NA,
                            R = NA, Q = NA, ILpreD = NA, D = NA, S = NA, ILpostD = NA, Cost = NA)
    inventory$D <- rep(test_data[,i])
    for (k in 1:h) {
      inventory$Sigma[k] <- sd(data[1:(t+k-1),i])
    }
    inventory$Mu <- predictions[,i]
    
    for (k in 1:h) {
      if (inventory$Mu[k] < 0) {
        inventory$R[k] <- 0
        next
      }
      probLTD <- as.data.frame(table(LTD[,k]))
      colnames(probLTD) <- c("Val", "Freq")
      probLTD$Prob <- probLTD$Freq / 1000
      probLTD$Val <- as.numeric(as.vector(probLTD$Val))
      
      R <- 0
      lossTarget <- (1 - r) * inventory$Mu[k]
      tempQ <- inventory$Mu[k]
      shortagePerCycle <- ESC_WILL(probLTD, inventory$Mu[k], R) #- ESC_NORM(inventory$Mu[k], inventory$Sigma[k], R + tempQ)
      while (shortagePerCycle > lossTarget) {
        R <- R + 1
        shortagePerCycle <- ESC_WILL(probLTD, inventory$Mu[k], R) #- ESC_NORM(inventory$Mu[k], inventory$Sigma[k], R + tempQ)
        if (shortagePerCycle == -1) {
          break
        }
      }
      inventory$R[k] <- R
    }
    
    inventory$Q[1] <- 0                              # Initial order quantity
    inventory$ILpreD[1] <- inventory$R[1]
    inventory$S[1] <- min(inventory$ILpreD[1], inventory$D[1])
    inventory$ILpostD[1] <- inventory$ILpreD[1] - inventory$D[1]
    
    for (k in 2:h) {
      inventory$Q[k] <- max(0, inventory$R[k] - inventory$ILpostD[k-1])
      inventory$ILpreD[k] <- max(inventory$ILpostD[k-1], inventory$R[k])
      inventory$S[k] <- min(inventory$ILpreD[k], inventory$D[k])
      inventory$ILpostD[k] <- inventory$ILpreD[k] - inventory$D[k]
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
  
  listLTD[[i]] <- LTD
}

pathPredictions <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\Non_Params_methods\\Willemain_", dataName, ".Rda", sep = "")
save(predictions, file = pathPredictions)

pathLTD <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\Non_Params_methods\\listLTD_", dataName, ".Rda", sep = "")
save(listLTD, file = pathLTD)


totalAFR <- rowSums(TotalSupply) / sum(TotalDemand)
avgAFR <- rowMeans(FillRates)
holdingCost <- rowSums(HoldingCosts)
ServiceLevel <- data.frame(avgAFR, totalAFR, holdingCost, targetFillRates, Method)
colnames(ServiceLevel) <- c("AchievedFillRates_Avg", "AchievedFillRates_Total", "HoldingCosts", "TargetFillRates", "Method")
file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\Non_Params_methods_Service_level_", Method, "_", dataName, ".Rda", sep = "")
save(ServiceLevel, FillRates, TotalDemand, TotalSupply, file = file_name)

#### SIM1 ####
Willemain_SIM1 <- predictions
Willemain_FAM_SIM1 <- FAM(trainSIM1, testSIM1, Willemain_SIM1)
Willeamin_Service_level_SIM1 <- ServiceLevel
Willemain_IPM_SIM1 <- TotalSupply

#### SIM2 ####
Willemain_SIM2 <- predictions
Willemain_FAM_SIM2 <- FAM(trainSIM2, testSIM2, Willemain_SIM2)
Willemain_Service_level_SIM2 <- ServiceLevel
Willemain_IPM_SIM2 <- TotalSupply

#### SIM3 ####
Willemain_SIM3 <- predictions
Willemain_FAM_SIM3 <- FAM(trainSIM3, testSIM3, Willemain_SIM3)
Willemain_Service_level_SIM3 <- ServiceLevel
Willemain_IPM_SIM3 <- TotalSupply

#### SIM4 ####
Willemain_SIM4 <- predictions
Willemain_FAM_SIM4 <- FAM(trainSIM4, testSIM4, Willemain_SIM4)
Willemain_Service_level_SIM4 <- ServiceLevel
Willemain_IPM_SIM4 <- TotalSupply

#### MAN ####
Willemain_MAN <- predictions
Willemain_FAM_MAN <- FAM(trainMAN, testMAN, Willemain_MAN)
Willemain_Service_level_MAN <- ServiceLevel
Willemain_IPM_MAN <- TotalSupply

#### BRAF ####
Willemain_BRAF <- predictions
Willemain_FAM_BRAF <- FAM(trainBRAF, testBRAF, Willemain_BRAF)
Willemain_Service_level_BRAF <- ServiceLevel
Willemain_IPM_BRAF <- TotalSupply

#### AUTO ####
Willemain_AUTO <- predictions
Willemain_FAM_AUTO <- FAM(trainAUTO, testAUTO, Willemain_AUTO)
Willemain_Service_level_AUTO <- ServiceLevel
Willemain_IPM_AUTO <- TotalSupply

#### OIL ####
Willemain_OIL <- predictions
Willemain_FAM_OIL <- FAM(trainOIL, testOIL, Willemain_OIL)
Willemain_Service_level_OIL <- ServiceLevel
Willemain_IPM_OIL <- TotalSupply


#####################Quantile regression#######################################
library(quantreg)

QuantileRegression_Perf <- function(train_data, test_data, data_name, quantiles = c(0.5, 0.85, 0.9, 0.95, 0.99)) {
  h <- nrow(test_data)                               # Number of periods to be predicted (for each item)
  numItems <- ncol(train_data)                       # Number of items
  predictions <- array(dim = c(h, numItems, length(quantiles)))   # Array storing the predictions
  
  # For each item
  for (j in 1:numItems) {
    x <- rbind(train_data[j], test_data[j])          # All demands (train + test) of item j
    t <- nrow(train_data)                            # Number of periods used for predicting (updated every iteration below)
    
    for (i in 1:h) {
      # Fit quantile regression model for each quantile
      for (q in 1:length(quantiles)) {
        quantile_val <- quantiles[q]
        qr_fit <- rq(x[1:t,] ~ 1, tau = quantile_val, data = data.frame(y = x[1:t,]))
        predictions[i, j, q] <- predict(qr_fit, newdata = data.frame(y = NA))
      }
      t = t + 1
    }
    
    print(paste("Done item", j))
  }
  
  # Store predictions as data frame
  predictions <- as.data.frame(predictions)
  
  # Saving
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\QuantileRegression_", data_name, ".Rda", sep = "")
  save(predictions, file = file_name)
  return(predictions)
}

#### SIM1 ####
QuantileSIM1 <- QuantileRegression_Perf(trainSIM1, testSIM1, "SIM1")
Quantile_FAM_SIM1 <- FAM(trainSIM1, testSIM1, QuantileSIM1)
Quantile_IPM_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1, QuantileSIM1, pricesSIM1, "SIM1", "QR")
Quantile_FAM_SIM1

#### SIM2 ####
QuantileSIM2 <- QuantileRegression_Perf(trainSIM2, testSIM2, "SIM2")
Quantile_FAM_SIM2 <- FAM(trainSIM2, testSIM2, QuantileSIM2)
Quantile_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2, QuantileSIM2, pricesSIM2, "SIM2", "QR")
Quantile_FAM_SIM2

#### SIM3 ####
QuantileSIM3 <- QuantileRegression_Perf(trainSIM3, testSIM3, "SIM3")
Quantile_FAM_SIM3 <- FAM(trainSIM3, testSIM3, QuantileSIM3)
Quantile_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3, QuantileSIM3, pricesSIM3, "SIM3", "QR")
Quantile_FAM_SIM3

#### SIM4 ####
QuantileSIM4 <- QuantileRegression_Perf(trainSIM4, testSIM4, "SIM4")
Quantile_FAM_SIM4 <- FAM(trainSIM4, testSIM4, QuantileSIM4)
Quantile_IPM_SIM4 <- IPM_GAMMA(SIM4, trainSIM4, testSIM4, QuantileSIM4, pricesSIM4, "SIM4", "QR")
Quantile_FAM_SIM4

#### MAN ####
QuantileMAN <- QuantileRegression_Perf(trainMAN, testMAN, "MAN")
Quantile_FAM_MAN <- FAM(trainMAN, testMAN, QuantileMAN)
Quantile_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN, QuantileMAN, pricesMAN, "MAN", "QR")
Quantile_FAM_MAN

#### BRAF ####
QuantileBRAF <- QuantileRegression_Perf(trainBRAF, testBRAF, "BRAF")
Quantile_FAM_BRAF <- FAM(trainBRAF, testBRAF, QuantileBRAF)
Quantile_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF, QuantileBRAF, pricesBRAF, "BRAF", "QR")
Quantile_FAM_BRAF

#### AUTO ####
QuantileAUTO <- QuantileRegression_Perf(trainAUTO, testAUTO, "AUTO")
Quantile_FAM_AUTO <- FAM(trainAUTO, testAUTO, QuantileAUTO)
Quantile_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO, QuantileAUTO, pricesAUTO, "AUTO", "QR")
Quantile_FAM_AUTO

#### OIL ####
QuantileOIL <- QuantileRegression_Perf(trainOIL, testOIL, "OIL")
Quantile_FAM_OIL <- FAM(trainOIL, testOIL, QuantileOIL)
Quantile_IPM_OIL <- IPM_GAMMA(OIL, trainOIL, testOIL, QuantileOIL, pricesOIL, "OIL", "QR")
Quantile_FAM_OIL








