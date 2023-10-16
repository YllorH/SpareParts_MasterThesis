targetFillRates <- c(75:99)
targetFillRates <- targetFillRates / 100
nTFR <- length(targetFillRates)

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
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\IPM", Method, "_", data_name, ".Rda", sep = "")
  save(ServiceLevel, FillRates, TotalDemand, TotalSupply, file = file_name)
  
  return(ServiceLevel)
}

### Remove scientific notations in plots ###
# Function to set numbers with marks and without scientific notation
marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)

#####TRADE OFF CURVES ---- SIM1 ######
servicelevelsSIM1 <- rbind(Croston_IPM_SIM1, SBA_IPM_SIM1, DLP_IPM_SIM1, Willemain_Service_level_SIM1, QR_IPM_SIM1, MLP_IPM_SIM1, LSTM_IPM_SIM1, LightGBM_IPM_SIM1, RF_IPM_SIM1)
# Saving the service level data.
save(servicelevelsSIM1,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_SIM1.Rda")
kable(servicelevelsSIM1, "latex")

ggplot(servicelevelsSIM1, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_x_continuous(limits = c(0.76, 0.99), labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingSIM1.png")

ggplot(servicelevelsSIM1, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_x_continuous(limits = c(0.75, 0.99), labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingSIM1.png")

#####TRADE OFF CURVES ---- SIM2 ######
servicelevelsSIM2 <- rbind(Croston_IPM_SIM2, SBA_IPM_SIM2, DLP_IPM_SIM2, Willemain_Service_level_SIM2, QR_IPM_SIM2, MLP_IPM_SIM2, LSTM_IPM_SIM2, LightGBM_IPM_SIM2, RF_IPM_SIM2)
# Saving the service level data.
save(servicelevelsSIM2,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_SIM2.Rda")
kable(servicelevelsSIM2, "latex")

ggplot(servicelevelsSIM2, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  #scale_x_continuous(limits = c(0.76, 0.98), labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingSIM2.png")

ggplot(servicelevelsSIM2, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  #scale_x_continuous(limits = c(0.70, 0.98), labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingSIM2.png")

#####TRADE OFF CURVES ---- SIM3 ######
servicelevelsSIM3 <- rbind(Croston_IPM_SIM3, SBA_IPM_SIM3, DLP_IPM_SIM3, Willemain_Service_level_SIM3, Quantile_IPM_SIM3, MLP_IPM_SIM3, LSTM_IPM_SIM3, LightGBM_IPM_SIM3, RF_IPM_SIM3)
# Saving the service level data.
save(servicelevelsSIM3,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_SIM3.Rda")
kable(servicelevelsSIM3, "latex")

ggplot(servicelevelsSIM3, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingSIM3.png")

ggplot(servicelevelsSIM3, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingSIM3.png")

#####TRADE OFF CURVES ---- SIM4 ######
servicelevelsSIM4 <- rbind(Croston_IPM_SIM4, SBA_IPM_SIM4, DLP_IPM_SIM4, Willemain_Service_level_SIM4, Quantile_IPM_SIM4, MLP_IPM_SIM4, LSTM_IPM_SIM4, LightGBM_IPM_SIM4, RF_IPM_SIM4)
# Saving the service level data.
save(servicelevelsSIM4,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_SIM4.Rda")
kable(servicelevelsSIM4, "latex")

ggplot(servicelevelsSIM4, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingSIM4.png")

ggplot(servicelevelsSIM4, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5)  +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingSIM4.png")


#####TRADE OFF CURVES ---- MAN with QR######
servicelevelsMAN <- rbind(Croston_IPM_MAN, SBA_IPM_MAN, DLP_IPM_MAN, Willemain_Service_level_MAN, Quantile_IPM_MAN, MLP_IPM_MAN, LSTM_IPM_MAN, LightGBM_IPM_MAN, RF_IPM_MAN)
# Saving the service level data.
save(servicelevelsMAN,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_MAN.Rda")
kable(servicelevelsMAN, "latex")

ggplot(servicelevelsMAN, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingMAN_wQR.png")

ggplot(servicelevelsMAN, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingMAN_wQR.png")

#####TRADE OFF CURVES ---- MAN without QR ######
servicelevelsMAN <- rbind(Croston_IPM_MAN, SBA_IPM_MAN, DLP_IPM_MAN, Willemain_Service_level_MAN, MLP_IPM_MAN, LSTM_IPM_MAN, LightGBM_IPM_MAN, RF_IPM_MAN)
# Saving the service level data.
save(servicelevelsMAN,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_MAN.Rda")
kable(servicelevelsMAN, "latex")

ggplot(servicelevelsMAN, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingMAN.png")

ggplot(servicelevelsMAN, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingMAN.png")


#####TRADE OFF CURVES ---- BRAF with QR######
servicelevelsBRAF <- rbind(Croston_IPM_BRAF, SBA_IPM_BRAF, DLP_IPM_BRAF, Willemain_Service_level_BRAF, Quantile_IPM_BRAF, MLP_IPM_BRAF, LSTM_IPM_BRAF, LightGBM_IPM_BRAF, RF_IPM_BRAF)
# Saving the service level data.
save(servicelevelsBRAF,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_BRAF.Rda")
kable(servicelevelsBRAF, "latex")

ggplot(servicelevelsBRAF, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingBRAF_wQR.png")

ggplot(servicelevelsBRAF, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingBRAF_wQR.png")

#####TRADE OFF CURVES ---- BRAF without QR ######
servicelevelsBRAF <- rbind(Croston_IPM_BRAF, SBA_IPM_BRAF, DLP_IPM_BRAF, Willemain_Service_level_BRAF, MLP_IPM_BRAF, LSTM_IPM_BRAF, LightGBM_IPM_BRAF, RF_IPM_BRAF)
# Saving the service level data.
save(servicelevelsBRAF,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_BRAF.Rda")
kable(servicelevelsBRAF, "latex")

ggplot(servicelevelsBRAF, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingBRAF.png")

ggplot(servicelevelsBRAF, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingBRAF.png")

#####TRADE OFF CURVES ---- AUTO ######
servicelevelsAUTO <- rbind(Croston_IPM_AUTO, SBA_IPM_AUTO, DLP_IPM_AUTO, Willemain_Service_level_AUTO, Quantile_IPM_AUTO, MLP_IPM_AUTO, LSTM_IPM_AUTO, LightGBM_IPM_AUTO, RF_IPM_AUTO)
# Saving the service level data.
save(servicelevelsAUTO,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_AUTO.Rda")
kable(servicelevelsAUTO, "latex")

ggplot(servicelevelsAUTO, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingAUTO.png")

ggplot(servicelevelsAUTO, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingAUTO.png")

#####TRADE OFF CURVES ---- OIL with QR######
servicelevelsOIL <- rbind(Croston_IPM_OIL, SBA_IPM_OIL, DLP_IPM_OIL, Willemain_Service_level_OIL, Quantile_IPM_OIL, MLP_IPM_OIL, LSTM_IPM_OIL, LightGBM_IPM_OIL, RF_IPM_OIL)
# Saving the service level data.
save(servicelevelsOIL,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_OIL.Rda")
kable(servicelevelsOIL, "latex")

ggplot(servicelevelsOIL, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingOIL_wQR.png")

ggplot(servicelevelsOIL, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingOIL_wQR.png")

#####TRADE OFF CURVES ---- OIL without QR######
servicelevelsOIL <- rbind(Croston_IPM_OIL, SBA_IPM_OIL, DLP_IPM_OIL, Willemain_Service_level_OIL, MLP_IPM_OIL, LSTM_IPM_OIL, LightGBM_IPM_OIL, RF_IPM_OIL)
# Saving the service level data.
save(servicelevelsOIL,file="C:\\Users\\yllor\\OneDrive\\Bureau\\IPM_OIL.Rda")
kable(servicelevelsOIL, "latex")

ggplot(servicelevelsOIL, aes(x = AchievedFillRates_Avg, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Average achieved fill rate") + ylab("Inventory holding costs")
ggsave("Avg_achievedvsholdingOIL.png")

ggplot(servicelevelsOIL, aes(x = AchievedFillRates_Total, y = HoldingCosts, group = Method)) +
  geom_line(aes(color = Method), size = 0.5) +
  geom_point(size = 0.5) +
  scale_y_continuous(labels = marks_no_sci)+
  xlab("Total achieved fill rate") + ylab("Inventory holding costs")
ggsave("Tot_achievedvsholdingOIL.png")
