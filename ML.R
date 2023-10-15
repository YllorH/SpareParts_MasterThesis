###############################################ML METHODS###################################
################### MLP #######################
normalize <- function(x) {return ((x - min(x)) / (max(x) - min(x)))} #Supporting functions for MLP_Perf

MLP_Perf <- function(data, data_name) {
  # 1, Preparing training data (same)
  # 1.1, Normalize data 
  mldata <- data
  for (i in 1:ncol(data)){
    mldata[,i]  <- normalize(data[,i]) 
  }
  
  # 1.2, Splitting the data again into test and train data.
  
  sample = round(nrow(data)*.70) 
  train_data <- mldata[1:(sample), ]
  test_data <- mldata[-(1:(sample)), ]
  
  # mldata variable, length 6 with 5 input and 1 output.
  n <- nrow(train_data)
  mldata <- matrix(ncol = 6, nrow = ncol(train_data) * (n - 5))
  
  # For loop that creates a 6 column data set with the first 5 as inputs to the neural net and the 6th column as the target variable.
  for (i in 1:ncol(train_data)){
    for (t in 1:6) {
      startRow <- ((i - 1) * (n - 5) + 1)
      endRow <- i * (n - 5)
      mldata[startRow:endRow, t] <- t(train_data[t:(t + n - 6),i])
    }
  }
  
  mldata <- as.data.frame(mldata)
  
  
  # 2, Train the model
  model <- mlp(mldata[,1:5], mldata[,6], size=6, learnFuncParams=c(0.01, 0.025, 0.1, 0.25, 0.3), maxit= 200, early_stopping = TRUE ,verbose = TRUE)
  
  print("Done training model")
  
  
  # 3, Make predictions (same)
  
  # Creating the prediction matrix.
  h <- nrow(test_data)                               # Number of periods to be predicted (for each item)
  numItems <- ncol(train_data)                       # Number of items
  predictions <- matrix(nrow = h, ncol = numItems)   # Matrix storing the predictions
  
  # Predicting the forecast horizon with as input the last five values of a series and saving it as the predictions.
  for (j in 1:numItems) {
    input <- tail(train_data[,j], n=5)
    for (i in 1:nrow(test_data)) {
      predictions[i,j] <- predict(model, t(input))
      input[1:4] <- input[2:5]
      input[5] <- test_data[i,j]
    }
  }
  
  # Denormalizing the data for evaluation of forecasts.
  for (i in 1:ncol(predictions)) {
    predictions[,i] <- predictions[,i] * (max(data[,i]) - min(data[,i])) + min(data[,i])
  }
  
  file_name <- paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\ML_methods\\MLP_", data_name, ".Rda", sep = "")
  save(predictions, file = file_name)
  return(predictions)
}

###############Forecasting Accuracy Measures (FAM)########################
##GMAE function
GMAE_fun <- function(predicted, actual) {
  gmae <- prod(abs(predicted - actual))^(1 / length(predicted))
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


set.seed(8733)
##### SIM1 #####
MLP_SIM1 <- MLP_Perf(SIM1, "SIM1")
MLP_FAM_SIM1 <- FAM(trainSIM1, testSIM1, MLP_SIM1)
MLP_IPM_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1, MLP_SIM1, pricesSIM1, "SIM1", "MLP" )
MLP_FAM_SIM1

#### SIM2 ####
MLP_SIM2 <- MLP_Perf(SIM2, "SIM2")
MLP_FAM_SIM2 <- FAM(trainSIM2, testSIM2, MLP_SIM2)
MLP_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2, MLP_SIM2, pricesSIM2, "SIM2", "MLP" )
MLP_FAM_SIM2

#### SIM3 ####
MLP_SIM3 <- MLP_Perf(SIM3, "SIM3")
MLP_FAM_SIM3 <- FAM(trainSIM3, testSIM3, MLP_SIM3)
MLP_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3, MLP_SIM3, pricesSIM3, "SIM3", "MLP" )
MLP_FAM_SIM3

#### SIM4 ####
MLP_SIM4 <- MLP_Perf(SIM4, "SIM4")
MLP_FAM_SIM4 <- FAM(trainSIM4, testSIM4, MLP_SIM4)
MLP_IPM_SIM4 <- IPM_GAMMA(SIM4, trainSIM4, testSIM4, MLP_SIM4, pricesSIM4, "SIM4", "MLP" )
MLP_FAM_SIM4
#### MAN ####
MLP_MAN <- MLP_Perf(MAN, "MAN")
MLP_FAM_MAN <- FAM(trainMAN, testMAN, MLP_MAN)
MLP_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN, MLP_MAN, pricesMAN, "MAN", "MLP" )
MLP_FAM_MAN

#### BRAF ####
MLP_BRAF <- MLP_Perf(BRAF, "BRAF")
MLP_FAM_BRAF <- FAM(trainBRAF, testBRAF, MLP_BRAF)
MLP_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF, MLP_BRAF, pricesBRAF, "BRAF", "MLP" )
MLP_FAM_BRAF

#### AUTO ####
MLP_AUTO <- MLP_Perf(AUTO, "AUTO")
MLP_FAM_AUTO <- FAM(trainAUTO, testAUTO, MLP_AUTO)
MLP_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO, MLP_AUTO, pricesAUTO, "AUTO", "MLP" )
MLP_FAM_AUTO

#### OIL ####
MLP_OIL <- MLP_Perf(OIL, "OIL")
MLP_FAM_OIL <- FAM(trainOIL, testOIL, MLP_OIL)
MLP_IPM_OIL <- IPM_GAMMA(OIL, trainOIL, testOIL, MLP_OIL, pricesOIL, "OIL", "MLP")
MLP_FAM_OIL
MLP_IPM_OIL


######################LightGBM######################################
normalize <- function(x) {return ((x - min(x)) / (max(x) - min(x)))}

LightGBM_Perf <- function(data, data_name) {
  # 1, Preparing training data (same)
  # 1.1, Normalize data
  data <- SIM1
  mldata <- data
  for (i in 1:ncol(data)){
    mldata[,i]  <- normalize(data[,i])
  }
  
  # 1.2, Splitting the data again into test and train data.
  
  sample = round(nrow(data)*.70) 
  train_data <- mldata[1:(sample), ]
  test_data <- mldata[-(1:(sample)), ]
  
  # mldata variable, length 6 with 5 input and 1 output.
  n <- nrow(train_data)
  mldata <- matrix(ncol = 6, nrow = ncol(train_data) * (n - 5))
  
  # For loop that creates a 6 column data set with the first 5 as inputs to the neural net and the 6th column as the target variable.
  for (i in 1:ncol(train_data)){
    for (t in 1:6) {
      startRow <- ((i - 1) * (n - 5) + 1)
      endRow <- i * (n - 5)
      mldata[startRow:endRow, t] <- t(train_data[t:(t + n - 6),i])
    }
  }
  
  mldata <- as.data.frame(mldata)
  
  #2, Train the model
  
  #Setting the hyper-parameters for the LightGBM model
  
  #Grid search 
  num_leaves = seq(10,150,10)
  max_depth = round(log(num_leaves)/log(2),0)
  num_iterations = seq(500,15000,500)
  early_stopping_rounds = round(num_iterations*.1,0)
  learning_rate = seq(0.01,0.5,0.05)
  
  p <- expand.grid(objective = "regression", 
                   boosting = "gbdt",
                   force_row_wise = TRUE,
                   learning_rate = learning_rate,
                   max_depth = max_depth,
                   num_leaves = num_leaves,
                   num_iterations = num_iterations,
                   early_stopping_rounds = early_stopping_rounds,
                   bagging_freq =1,
                   lambda_l2=0.1)
  
  p <- unique(p) #the unique function creates multiple grids, such that all the given values for the hyperparams are used
  
  #Training the LightGBM algorithm.
  

  cv_results <- lgb.cv(
    params = list(
      objective = "regression",
      metric = "rmse",
      max_depth = p$max_depth[l],
      num_leaves = p$num_leaves[l],
      num_iterations = p$num_iterations[l],
      early_stopping_rounds=p$early_stopping_rounds[l],
      learning_rate=p$learning_rate[l],
      seed=8733), 
      data=as.matrix(mldata[,(1:5)]),
      label=mldata$V6)
  
  #get best iteration
  best_iter <- cv_results$best_iter
  
  # create lgb.Dataset object from training data
  dtrain <- lgb.Dataset(as.matrix(mldata[,(1:5)]), label = mldata$V6)
  
  # train model using lgb.Dataset object
  light_gbm_tuned <- lgb.train(
    params = list(
      objective = "regression",
      metric = "rmse",
      max_depth = p$max_depth[l],
      num_leaves = p$num_leaves[l],
      num_iterations = best_iter,
      learning_rate = p$learning_rate[l]),
    data = dtrain)
  
  print("Done training model")
  
  #3, Make predictions (same)
  
  #Creating the prediction matrix.
  h<-nrow(test_data)                               #Number of periods to be predicted (for each item)
  numItems<-ncol(train_data)                       #Number of items
  predictions<-matrix(nrow=h,ncol=numItems)   #Matrix storing the predictions
  
  #Predicting the forecast horizon with as input the last five values of a series and saving it as the predictions.
  for(j in seq_len(numItems)) {
    input<-tail(train_data[,j],n=5)
    for(i in seq_len(nrow(test_data))) {
      predictions[i,j]<-predict(light_gbm_tuned,t(input))
      input[1:4]<-input[2:5]
      input[5]<-test_data[i,j]
    }
  }
  
  #Denormalizing the data for evaluation of forecasts.
  for(i in seq_len(ncol(predictions))) {
    predictions[,i]<-predictions[,i]*(max(data[,i])-min(data[,i]))+min(data[,i])
  }
  
  file_name<-paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\ML_methods\\LightGBM_",data_name,".Rda",sep="")
  save(predictions,file=file_name)
  
  return(predictions)
}


#### SIM1 ####
LightGBM_SIM1 <- LightGBM_Perf(SIM1, "SIM1")
LightGBM_FAM_SIM1 <- FAM(trainSIM1, testSIM1, LightGBM_SIM1)
LigthGBM_IPM_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1, LightGBM_SIM1, pricesSIM1, "SIM1", "LightGBM")
LightGBM_FAM_SIM1

#### SIM2 ####
LightGBM_SIM2 <- LightGBM_Perf(SIM2, "SIM2")
LightGBM_FAM_SIM2 <- FAM(trainSIM2, testSIM2, LightGBM_SIM2)
LigthGBM_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2, LightGBM_SIM2, pricesSIM2, "SIM2", "LightGBM")
LightGBM_FAM_SIM2

#### SIM3 ####
LightGBM_SIM3 <- LightGBM_Perf(SIM3, "SIM3")
LightGBM_FAM_SIM3 <- FAM(trainSIM3, testSIM3, LightGBM_SIM3)
LightGBM_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3, LightGBM_SIM3, pricesSIM3, "SIM3", "LightGBM")
LightGBM_FAM_SIM3

#### SIM4 ####
LightGBM_SIM4 <- LightGBM_Perf(SIM4, "SIM4")
LightGBM_FAM_SIM4 <- FAM(trainSIM4, testSIM4, LightGBM_SIM4)
LightGBM_IPM_SIM4 <- IPM_GAMMA(SIM4, trainSIM4, testSIM4, LightGBM_SIM4, pricesSIM4, "SIM4", "LightGBM")
LightGBM_FAM_SIM4

#### MAN ####
LightGBM_MAN <- LightGBM_Perf(MAN, "MAN")
LightGBM_FAM_MAN <- FAM(trainMAN, testMAN, LightGBM_MAN)
LightGBM_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN, LightGBM_MAN, pricesMAN, "MAN", "LightGBM")
LightGBM_FAM_MAN

#### BRAF ####
LightGBM_BRAF <- LightGBM_Perf(BRAF, "BRAF")
LightGBM_FAM_BRAF <- FAM(trainBRAF, testBRAF, LightGBM_BRAF)
LightGBM_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF, LightGBM_BRAF, pricesBRAF, "BRAF", "LightGBM")
LightGBM_FAM_BRAF

#### AUTO ####
LightGBM_AUTO <- LightGBM_Perf(AUTO, "AUTO")
LightGBM_FAM_AUTO <- FAM(trainAUTO, testAUTO, LightGBM_AUTO)
LightGBM_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO, LightGBM_AUTO, pricesAUTO, "AUTO", "LightGBM")
LightGBM_FAM_AUTO

#### OIL ####
LightGBM_OIL <- LightGBM_Perf(OIL, "OIL")
LightGBM_FAM_OIL <- FAM(trainOIL, testOIL, LightGBM_OIL)
LightGBM_IPM_OIL <- IPM_GAMMA(OIL, trainOIL, testOIL, LightGBM_OIL, pricesOIL, "OIL", "LightGBM")
LightGBM_FAM_OIL


###########RANDOM FOREST#####################################
normalize <- function(x) {return ((x - min(x)) / (max(x) - min(x)))}
RandomForest_Perf <- function(data, data_name) {
  set.seed(8733)
  # 1, Preparing training data (same)
  # 1.1, Normalize data
  mldata <- data
  for (i in 1:ncol(data)){
    mldata[,i]  <- normalize(data[,i])
  }  
  # 1.2, Splitting the data again into test and train data.  
  sample = round(nrow(data)*.70) 
  train_data <- mldata[1:(sample), ]
  test_data <- mldata[-(1:(sample)), ]  
  # mldata variable, length 6 with 5 input and 1 output.
  n <- nrow(train_data)
  mldata <- matrix(ncol = 6, nrow = ncol(train_data) * (n - 5))  
  # For loop that creates a 6 column data set with the first 5 as inputs to the neural net and the 6th column as the target variable.
  for (i in 1:ncol(train_data)){
    for (t in 1:6) {
      startRow <- ((i - 1) * (n - 5) + 1)
      endRow <- i * (n - 5)
      mldata[startRow:endRow, t] <- t(train_data[t:(t + n - 6),i])
    }
  }  
  mldata <- as.data.frame(mldata)  
  
  #2, Train the model  

  #Find the best hyper parameters
  rf_model <- random_forest_parameters(base = mldata,
                                       target = "V6",
                                       model_type = "regression",
                                       ntree = c(100),
                                       mtry = c(1,2),
                                       nodesize = 5)
  #Train the best model
  best.ntree <- as.numeric(rf_model$best_ntree)
  best.mtry <- as.numeric(rf_model$best_mtry)
  best.nodesize <- as.numeric(rf_model$best_nodesize)
  
  final_rf_model <- randomForest(mldata$V6 ~ ., data = mldata,
                                 ntree = best.ntree,
                                 mtry = best.mtry,
                                 nodesize = best.nodesize)
  
  #3, Make predictions
  
  #Creating the prediction matrix.
  h<-nrow(test_data)                               #Number of periods to be predicted (for each item)
  numItems<-ncol(train_data)                       #Number of items
  predictions<-matrix(nrow=h,ncol=numItems)   #Matrix storing the predictions
  
  #Predicting the forecast horizon with as input the last five values of a series and saving it as the predictions.
  for(j in seq_len(numItems)) {
    input<-tail(train_data[,j],n=5)
    for(i in seq_len(nrow(test_data))) {
      predictions[i,j]<-predict(final_rf_model,t(input))
      input[1:4]<-input[2:5]
      input[5]<-test_data[i,j]
    }
  }
  
  #Denormalizing the data for evaluation of forecasts.
  for(i in seq_len(ncol(predictions))) {
    predictions[,i]<-predictions[,i]*(max(data[,i])-min(data[,i]))+min(data[,i])
  }
  
  file_name<-paste("C:\\Users\\yllor\\OneDrive\\Bureau\\Thesis\\SpareParts\\RF3.R",data_name,".Rda",sep="")
  save(predictions,file=file_name)
  
  return(predictions)
}

#### SIM1 ####
RF_SIM1 <- RandomForest_Perf(SIM1, "SIM1")
RF_FAM_SIM1 <- FAM(trainSIM1, testSIM1, RF_SIM1)
RF_IPM_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1, RF_SIM1, pricesSIM1, "SIM1", "RF")
RF_FAM_SIM1

#### SIM2 ####
RF_SIM2 <- RandomForest_Perf(SIM2, "SIM2")
RF_FAM_SIM2 <- FAM(trainSIM2, testSIM2, RF_SIM2)
RF_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2, RF_SIM2, pricesSIM2, "SIM2", "RF")
RF_FAM_SIM2


#### SIM3 ####
RF_SIM3 <- RandomForest_Perf(SIM3, "SIM3")
RF_FAM_SIM3 <- FAM(trainSIM3, testSIM3, RF_SIM3)
RF_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3, RF_SIM3, pricesSIM3, "SIM3", "RF")
RF_FAM_SIM3

#### SIM4 ####
RF_SIM4 <- RandomForest_Perf(SIM4, "SIM4")
RF_FAM_SIM4 <- FAM(trainSIM4, testSIM4, RF_SIM4)
RF_IPM_SIM4 <- IPM_GAMMA(SIM4, trainSIM4, testSIM4, RF_SIM4, pricesSIM4, "SIM4", "RF")
RF_FAM_SIM4

#### MAN ####
RF_MAN <- RandomForest_Perf(MAN, "MAN")
RF_FAM_MAN <- FAM(trainMAN, testMAN, RF_MAN)
RF_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN, RF_MAN, pricesMAN, "MAN", "RF")
RF_FAM_MAN

#### BRAF ####
RF_BRAF <- RandomForest_Perf(BRAF, "BRAF")
RF_FAM_BRAF <- FAM(trainBRAF, testBRAF, RF_BRAF)
RF_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF, RF_BRAF, pricesBRAF, "BRAF", "RF")
RF_FAM_BRAF

#### AUTO ####
RF_AUTO <- RandomForest_Perf(AUTO, "AUTO")
RF_FAM_AUTO <- FAM(trainAUTO, testAUTO, RF_AUTO)
RF_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO, RF_AUTO, pricesAUTO, "AUTO", "RF")
RF_FAM_AUTO

#### OIL ####
RF_OIL <- RandomForest_Perf(OIL, "OIL")
RF_FAM_OIL <- FAM(trainOIL, testOIL, RF_OIL)
RF_IPM_OIL <- IPM_GAMMA(OIL, trainOIL, testOIL, RF_OIL, pricesOIL, "OIL", "RF")
RF_FAM_OIL


#####################LSTM####################################
set.seed(8733)
# Define normalization function
normalize <- function(x) {return ((x - min(x)) / (max(x) - min(x)))}

# Define LSTM_Perf function
LSTM_Perf <- function(data, data_name) {
  # Normalize the data
  mldata <- data
  for (i in 1:ncol(data)){
    mldata[,i]  <- normalize(data[,i])
  }
  
  # 1.2 Splitting the data again into test and train data.
  sample = round(nrow(data)*.70)
  train_data <- mldata[1:(sample), ]
  test_data <- mldata[-(1:(sample)), ]
  
  # Data preparation
  #lag  set to 2 for AUTO
  lag <- 5 # Nr of previous time steps used as input variables for each sequence to predict next time period
  delay <- 1 # How far in the future the model should predict (prediction horizon set to 1, predict the next value (time step) into the future)
  n <- 1  # Next steps: how many future time steps to forecast
  
  # Get the number of samples (N) and the dimension of data (d)
  N <- nrow(train_data)
  d <- ncol(train_data)
  
  # Calculate k (how many complete samples (input - output pairs) can be created from the data set)
  k <- N - (lag + delay + n)
  
  ## Input and output slices: they represent the sequences of data that will be used as input into the LSTM and as output
  # Creating matrix for the input slices
  input_slices <- matrix(nrow = k, ncol = lag)
  for (i in 1:k) {
    input_slices[i, ] <- i:(i + lag - 1)
  }
  
  # Creating matrix for the output slices
  output_slices <- matrix(nrow = k, ncol = n)
  for (i in 1:k) {
    output_slices[i, ] <- (i + lag + delay):(i + lag + delay + n - 1)
  }
  
  
  # Creating the X and Y variable to use in LSTM
  X <- lapply(1:ncol(data), function(i) {
    sapply(1:k, function(j) {
      data[(j):(j+lag-1), i]
    })
  })
  
  Y <- lapply(1:ncol(data), function(i) {
    sapply(1:k, function(j) {
      data[(j+lag+delay):(j+lag+delay+n-1), i]
    })
  })
  
  # Convert X and Y to arrays
  X <- array(unlist(X), dim = c(k, lag, d))
  Y <- array(unlist(Y), dim = c(k, n, d))
  
  # Same pre-processing for test data
  
  N_test <- nrow(test_data)
  k_test <- N_test - (lag + delay + n)
  
  input_slices_test <- matrix(nrow = k_test, ncol = lag)
  for (i in 1:k_test) {
    input_slices_test[i, ] <- i:(i + lag - 1)
  }
  
  output_slices_test <- matrix(nrow = k_test, ncol = n)
  for (i in 1:k_test) {
    output_slices_test[i, ] <- (i + lag + delay):(i + lag + delay + n - 1)
  }
  
  X_test <- lapply(1:ncol(test_data), function(i) {
    sapply(1:k_test, function(j) {
      test_data[(j):(j+lag-1), i]
    })
  })
  
  Y_test <- lapply(1:ncol(test_data), function(i) {
    sapply(1:k_test, function(j) {
      test_data[(j+lag+delay):(j+lag+delay+n-1), i]
    })
  })
  
  X_test <- array(unlist(X_test), dim = c(k_test, lag, d))
  Y_test <- array(unlist(Y_test), dim = c(k_test, n, d))
  
  # Define model (2 LSTM layers with 50 units each, dropout layer after each
  #LSTM layer to prevent overfitting, a dense layer with units (=cells) equal to the number of variables in the data)
  model <- keras_model_sequential() %>%
    layer_lstm(units = 50, return_sequences = TRUE, input_shape = c(lag, d)) %>%
    layer_dropout(rate = 0.2) %>%
    layer_lstm(units = 50, return_sequences = FALSE) %>%
    layer_dropout(rate = 0.2) %>%
    layer_dense(units = d)
  
  # Compile model
  model %>% compile(
    optimizer = optimizer_adam(lr = 0.001),
    loss = 'mean_squared_error'
  )
  
  # Train model
  history <- model %>% fit(
    X, Y,
    epochs = 10,
    batch_size = 32,
    validation_split = 0.2
  )
  
  # Evaluate model on test data
  score <- model %>% evaluate(X_test, Y_test, batch_size = 32)
  
  # Make predictions on test data
  predictions <- model %>% predict(X_test)
  
  # Denormalize predictions
  for (i in seq_len(ncol(predictions))) {
    predictions[,i] <- predictions[,i] * (max(data[,i]) - min(data[,i])) + min(data[,i])
  }
  
  return(predictions)
}

#### SIM1 ####
LSTM_SIM1 <- LSTM_Perf(SIM1, "SIM1")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testSIM1_adjusted <- testSIM1[1:11, ]
train
LSTM_FAM_SIM1 <- FAM(trainSIM1, testSIM1_adjusted, LSTM_SIM1)
LSTM_IPM_SIM1 <- IPM_GAMMA(SIM1, trainSIM1, testSIM1_adjusted, LSTM_SIM1, pricesSIM1, "SIM1", "LSTM")
LSTM_FAM_SIM1

#### SIM2 ####
LSTM_SIM2 <- LSTM_Perf(SIM2, "SIM2")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testSIM2_adjusted <- testSIM2[1:11, ]
LSTM_FAM_SIM2 <- FAM(trainSIM2, testSIM2_adjusted, LSTM_SIM2)
LSTM_IPM_SIM2 <- IPM_GAMMA(SIM2, trainSIM2, testSIM2_adjusted, LSTM_SIM2, pricesSIM2, "SIM2", "LSTM")
LSTM_FAM_SIM2

#### SIM3 ####
LSTM_SIM3 <- LSTM_Perf(SIM3, "SIM3")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testSIM3_adjusted <- testSIM3[1:11, ]
LSTM_FAM_SIM3 <- FAM(trainSIM3, testSIM3_adjusted, LSTM_SIM3)
LSTM_IPM_SIM3 <- IPM_GAMMA(SIM3, trainSIM3, testSIM3_adjusted, LSTM_SIM3, pricesSIM3, "SIM3", "LSTM")
LSTM_FAM_SIM3

#### SIM4 ####
LSTM_SIM4 <- LSTM_Perf(SIM4, "SIM4")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testSIM4_adjusted <- testSIM4[1:11, ]
LSTM_FAM_SIM4 <- FAM(trainSIM4, testSIM4_adjusted, LSTM_SIM4)
LSTM_IPM_SIM4 <- IPM_GAMMA(SIM4, trainSIM4, testSIM4_adjusted, LSTM_SIM4, pricesSIM4, "SIM4", "LSTM")
LSTM_FAM_SIM4

#### MAN ####
LSTM_MAN <- LSTM_Perf(MAN, "MAN")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testMAN_adjusted <- testMAN[1:38, ]
LSTM_FAM_MAN <- FAM(trainMAN, testMAN_adjusted, LSTM_MAN)
LSTM_IPM_MAN <- IPM_GAMMA(MAN, trainMAN, testMAN_adjusted, LSTM_MAN, pricesMAN, "MAN", "LSTM")
LSTM_FAM_MAN

#### BRAF ####
LSTM_BRAF <- LSTM_Perf(BRAF, "BRAF")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testBRAF_adjusted <- testBRAF[1:18, ]
LSTM_FAM_BRAF <- FAM(trainBRAF, testBRAF_adjusted, LSTM_BRAF)
LSTM_IPM_BRAF <- IPM_GAMMA(BRAF, trainBRAF, testBRAF_adjusted, LSTM_BRAF, pricesBRAF, "BRAF", "LSTM")
LSTM_FAM_BRAF

#### AUTO ####
LSTM_AUTO <- LSTM_Perf(AUTO, "AUTO")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testAUTO_adjusted <- testAUTO[1:2, ]
LSTM_FAM_AUTO <- FAM(trainAUTO, testAUTO_adjusted, LSTM_AUTO)
LSTM_IPM_AUTO <- IPM_GAMMA(AUTO, trainAUTO, testAUTO_adjusted, LSTM_AUTO, pricesAUTO, "AUTO", "LSTM")
LSTM_FAM_AUTO

#### OIL ####
LSTM_OIL <- LSTM_Perf(OIL, "OIL")
# Adjusting the length of the test data to match the length of predictions created by LSTM_Function
testOIL_adjusted <- testOIL[1:10, ]
LSTM_FAM_OIL <- FAM(trainOIL, testOIL_adjusted, LSTM_OIL)
LSTM_IPM_OIL <- IPM_GAMMA(OIL, trainOIL, testOIL_adjusted, LSTM_OIL, pricesOIL, "OIL", "LSTM")
LSTM_FAM_OIL
