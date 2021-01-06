#Load required library
library(qdapRegex)
library(varhandle)
library(tmle)
library(SuperLearner)
library(corrplot)
R <- 1000
set.seed(2310)
#Empty vectors 
Gcomputation_Case1  <- rep(NA,R)
IPTW_Case1  <- rep(NA,R)
AIPW_Case1 <- rep(NA,R)
TMLE_Case1  <- rep(NA,R)
Gcomputation_Case2     <- rep(NA,R)
IPTW_Case2 <- rep(NA,R)
AIPW_Case2 <- rep(NA,R)
TMLE_Case2  <- rep(NA,R)
Gcomputation_Case3  <- rep(NA,R)
IPTW_Case3  <- rep(NA,R)
AIPW_Case3 <- rep(NA,R)
TMLE_Case3  <- rep(NA,R)
Gcomputation_Case4  <- rep(NA,R)
IPTW_Case4  <- rep(NA,R)
AIPW_Case4 <- rep(NA,R)
TMLE_Case4 <- rep(NA,R)
Gcomputation_Case5  <- rep(NA,R)
IPTW_Case5  <- rep(NA,R)
AIPW_Case5 <- rep(NA,R)
TMLE_Case5 <- rep(NA,R)
TMLE_Case5 <- rep(NA,R)
Gcomputation_Case6  <- rep(NA,R)
IPTW_Case6  <- rep(NA,R)
AIPW_Case6 <- rep(NA,R)
TMLE_Case6 <- rep(NA,R)
TMLE_Case6 <- rep(NA,R)
#Define function for comparision
for(r in 1:R){  
  print(paste("This is simulation run number",r)) 
  ObsData <- generateData(n= 200) 
  ##------------------------------------------------------------------------------------------------------##
  ##--------------------------------------Case1 : Both models misspecified---------------------------------##
  ##------------------------------------------------------------------------------------------------------##
  gm <- glm(Y ~ A + w1 + w2 + w3 + w4, family = binomial, data = ObsData)
  # Prediction for A, A = 1 and, A = 0
  QAW_0 <- predict(gm, type = "response")
  Q1W_0 = predict(gm, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")
  Q0W_0 = predict(gm, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")
  # Plug in G-model
  Gcomputation_Case1[r]<-mean(Q1W_0 - Q0W_0)
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  # Step 2 estimation and prediction of the propensity score (ps)
  psm <- glm(A ~ w1 + w2 + w3 + w4,  family = binomial, data = ObsData)
  P.1W = predict(psm, type = "response")
  #For each individual, assign the predicted probability of being not a smoker to P.0W
  P.0W<- 1- P.1W
  
  #calculate the Clever Covariate for each individual
  H.AW<- as.numeric(ObsData$A==1)/P.1W - as.numeric(ObsData$A==0)/P.0W
  #Clever Covariate for those who are smoker - the interpretation of the inverse
  #of the probability of being a smoker (H.1W= H given Smoker==1 and covariates W)
  H.1W<- 1/P.1W
  #And for those who are not
  H.0W<- -1/P.0W
  #This party is basically calculating the IPTW_ATE for later comparison
  IPTW_Case1[r] <-mean(as.numeric(ObsData$A==1)/P.1W*ObsData$Y) - mean(as.numeric(ObsData$A==0)/P.0W*ObsData$Y)
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  Q1W_AIPW <- mean((as.numeric(ObsData$A==1)) * (ObsData$Y - Q1W_0) / P.1W + Q1W_0)
  Q0W_AIPW <- mean((as.numeric(ObsData$A==0)) * (ObsData$Y - Q0W_0) / P.0W + Q0W_0)
  AIPW_Case1[r]<-Q1W_AIPW-Q0W_AIPW
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################
  logitUpdate<- glm(ObsData$Y ~ -1 +offset(qlogis(QAW_0)) + H.AW, family='binomial')
  summary(logitUpdate)
  epsilon <- logitUpdate$coef
  epsilon
  # Then update the initial expected outcome :
  # logit(E(Y|A,W)_updated) = logit(E(Y|A,W))+ epsilon*H(A|W)
  QAW_0.updated<- plogis(qlogis(QAW_0)+ epsilon*H.AW)
  Q1W_0.updated<- plogis(qlogis(Q1W_0)+ epsilon*H.1W)
  Q0W_0.updated<- plogis(qlogis(Q0W_0)+ epsilon*H.0W)
  TMLE_Case1[r]<- mean(Q1W_0.updated - Q0W_0.updated)
  
  ##------------------------------------------------------------------------------------------------------##
  ##----------------------Case2 : Outcome model correctly specified--------------##
  ##------------------------------------------------------------------------------------------------------##
  gm <- glm(Y ~ A + w1 + w2 + w3+w4+w3*w4, family = binomial, data = ObsData)
  # Prediction for A, A = 1 and, A = 0
  QAW_0 <- predict(gm, type = "response")
  Q1W_0 = predict(gm, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")
  Q0W_0 = predict(gm, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")
  # Plug in G-model
  Gcomputation_Case2[r]<-mean(Q1W_0 - Q0W_0)
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  # Step 2 estimation and prediction of the propensity score (ps)
  psm <- glm(A ~ w1 + w2 + w3,  family = binomial, data = ObsData)
  P.1W = predict(psm, type = "response")
  #For each individual, assign the predicted probability of being not a smoker to P.0W
  P.0W<- 1- P.1W
  #calculate the Clever Covariate for each individual
  H.AW<- as.numeric(ObsData$A==1)/P.1W - as.numeric(ObsData$A==0)/P.0W
  #Clever Covariate for those who are smoker - the interpretation of the inverse
  #of the probability of being a smoker (H.1W= H given Smoker==1 and covariates W)
  H.1W<- 1/P.1W
  #And for those who are not
  H.0W<- -1/P.0W
  #This party is basically calculating the IPTW_ATE for later comparison
  IPTW_Case2[r] <-mean(as.numeric(ObsData$A==1)/P.1W*ObsData$Y) - mean(as.numeric(ObsData$A==0)/P.0W*ObsData$Y)
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  Q1W_AIPW <- mean((as.numeric(ObsData$A==1)) * (ObsData$Y - Q1W_0) / P.1W + Q1W_0)
  Q0W_AIPW <- mean((as.numeric(ObsData$A==0)) * (ObsData$Y - Q0W_0) / P.0W + Q0W_0)
  AIPW_Case2[r]<-Q1W_AIPW-Q0W_AIPW
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  logitUpdate<- glm(ObsData$Y ~ -1 +offset(qlogis(QAW_0)) + H.AW, family='binomial')
  summary(logitUpdate)
  epsilon <- logitUpdate$coef
  epsilon
  # Then update the initial expected outcome :
  # logit(E(Y|A,W)_updated) = logit(E(Y|A,W))+ epsilon*H(A|W)
  QAW_0.updated<- plogis(qlogis(QAW_0)+ epsilon*H.AW)
  Q1W_0.updated<- plogis(qlogis(Q1W_0)+ epsilon*H.1W)
  Q0W_0.updated<- plogis(qlogis(Q0W_0)+ epsilon*H.0W)
  TMLE_Case2[r]<- mean(Q1W_0.updated - Q0W_0.updated)
  
  
  ##----------------------------Case3: Using Super Learner-----------------------##
  ##-----------------------------------Step 1:-----------------------------------##
  ##-------------------------Estimate the initial E0(Y |A, W)--------------------##
  ##-----------------------------------------------------------------------------##
  #Define the predictor for the Superlearner (Exposure + Baseline Covariates)
  X<- subset(ObsData, select=c(A, w1, w2, w3,w4))
  
  #Specify the Superlearner Library
  set.seed(2310)
  SL.library <- c("SL.glm", "SL.step","SL.glm.interaction")
  #Estimate the conditional expectation of outcome, given exposure(A) and covariates(w)
  EY.AW_SL<- SuperLearner(Y=ObsData$Y, X=X, SL.library=SL.library, family="binomial")
  
  #Obtain the initial estimate of the expected outcome, given exposure(Smoker), and covariates(W)
  EY.AW <- predict(EY.AW_SL, newdata =ObsData)$pred
  
  #Similarly, obtain the expected outcome for those who under exposure
  EY.1W<- predict(EY.AW_SL, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")$pred
  
  #And those who not
  EY.0W<- predict(EY.AW_SL, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")$pred
  
  
  #Calculating the ATE following G-Computation method (simple subtitution)
  #as the difference in potential outcomes Y1, Y0, corresponding to A=1,0 respectively
  Gcomputation_Case3[r]<- mean(EY.1W- EY.0W)
  
  ##########################################################################################################
  ##########################################################################################################
  ########################################################################################################## 
  
  #Using Superlearner to predict the propensity score
  #which is the probability of the exposure - mother being smoker denoted P0(A=| 1 X)
  #given a set of covariates (W- our predictors)
  P.AW_SL<- SuperLearner(Y=ObsData$A, X=subset(X, select=c(w1, w2, w3, w4)),
                         SL.library=SL.library, family="binomial")
  
  #For each individual, assign the predicted probability of being a smoker to P.1W
  P.1W<- P.AW_SL$SL.predict
  
  #For each individual, assign the predicted probability of being not a smoker to P.0W
  P.0W<- 1- P.1W
  #calculate the Clever Covariate for each individual
  H.AW<- as.numeric(ObsData$A==1)/P.1W - as.numeric(ObsData$A==0)/P.0W
  
  #Clever Covariate for those who are smoker - the interpretation of the inverse
  #of the probability of being a smoker (H.1W= H given Smoker==1 and covariates W)
  H.1W<- 1/P.1W
  
  #And for those who are not
  H.0W<- -1/P.0W
  
  #This party is basically calculating the IPTW_ATE for later comparison
  IPTW_Case3[r] <-mean(as.numeric(ObsData$A==1)/P.1W*ObsData$Y) - mean(as.numeric(ObsData$A==0)/P.0W*ObsData$Y)
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  Q1W_AIPW <- mean((as.numeric(ObsData$A==1)) * (ObsData$Y - EY.1W) / P.1W + EY.1W)
  Q0W_AIPW <- mean((as.numeric(ObsData$A==0)) * (ObsData$Y - EY.0W) / P.0W + EY.0W)
  AIPW_Case3[r]<-Q1W_AIPW-Q0W_AIPW
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################
  logitUpdate<- glm(ObsData$Y ~ -1 +offset(qlogis(EY.AW)) + H.AW, family='binomial')
  summary(logitUpdate)
  
  # epsilon denote the resulting maximum likelihood estimate of the coefficient 
  # on the clever covariate H.AW
  epsilon <- logitUpdate$coef
  epsilon
  
  # Then update the initial expected outcome :
  # logit(E(Y|A,W)_updated) = logit(E(Y|A,W))+ epsilon*H(A|W)
  EY.AW.updated<- plogis(qlogis(EY.AW)+ epsilon*H.AW)
  EY.1W.updated<- plogis(qlogis(EY.1W)+ epsilon*H.1W)
  EY.0W.updated<- plogis(qlogis(EY.0W)+ epsilon*H.0W)
  
  # Estimate the statistical parameter by substituting the targeted predictions into
  # G-Computation's formula
  TMLE_Case3[r]<-mean(EY.1W.updated-EY.0W.updated)
  
  
  ##----------------------------Case4: Using Super Learner-----------------------##
  ##-----------------------------------Step 1:-----------------------------------##
  ##-------------------------Estimate the initial E0(Y |A, W)--------------------##
  ##-----------------------------------------------------------------------------##
  #Define the predictor for the Superlearner (Exposure + Baseline Covariates)
  X<- subset(ObsData, select=c(A, w1, w2, w3, w4))
  
  #Specify the Superlearner Library
  set.seed(2310)
  SL.library <- c("SL.gam","SL.randomForest","SL.rpart")
  #Estimate the conditional expectation of outcome, given exposure(A) and covariates(w)
  EY.AW_SL<- SuperLearner(Y=ObsData$Y, X=X, SL.library=SL.library, family="binomial")
  EY.AW_SL
  
  #Obtain the initial estimate of the expected outcome, given exposure(Smoker), and covariates(W)
  EY.AW <- predict(EY.AW_SL, newdata =ObsData)$pred
  EY.AW
  
  #Similarly, obtain the expected outcome for those who under exposure
  EY.1W<- predict(EY.AW_SL, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")$pred
  
  #And those who not
  EY.0W<- predict(EY.AW_SL, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")$pred
  
  
  #Calculating the ATE following G-Computation method (simple subtitution)
  #as the difference in potential outcomes Y1, Y0, corresponding to A=1,0 respectively
  Gcomputation_Case4[r]<- mean(EY.1W- EY.0W)
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################
  
  #Using Superlearner to predict the propensity score
  #which is the probability of the exposure - mother being smoker denoted P0(A=| 1 X)
  #given a set of covariates (W- our predictors)
  set.seed(2310)
  P.AW_SL<- SuperLearner(Y=ObsData$A, X=subset(X, select=c(w1, w2, w3, w4)),
                         SL.library=SL.library, family="binomial")
  
  #For each individual, assign the predicted probability of being a smoker to P.1W
  P.1W<- P.AW_SL$SL.predict
  
  #For each individual, assign the predicted probability of being not a smoker to P.0W
  P.0W<- 1- P.1W
  #calculate the Clever Covariate for each individual
  H.AW<- as.numeric(ObsData$A==1)/P.1W - as.numeric(ObsData$A==0)/P.0W
  
  #Clever Covariate for those who are smoker - the interpretation of the inverse
  #of the probability of being a smoker (H.1W= H given Smoker==1 and covariates W)
  H.1W<- 1/P.1W
  
  #And for those who are not
  H.0W<- -1/P.0W
  
  #This party is basically calculating the IPTW_ATE for later comparison
  IPTW_Case4[r] <-mean(as.numeric(ObsData$A==1)/P.1W*ObsData$Y) - mean(as.numeric(ObsData$A==0)/P.0W*ObsData$Y)
  
  ##########################################################################################################
  ##########################################################################################################
  ########################################################################################################## 
  
  Q1W_AIPW <- mean((as.numeric(ObsData$A==1)) * (ObsData$Y - EY.1W) / P.1W + EY.1W)
  Q0W_AIPW <- mean((as.numeric(ObsData$A==0)) * (ObsData$Y - EY.0W) / P.0W + EY.0W)
  AIPW_Case4[r]<-Q1W_AIPW-Q0W_AIPW
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################
  
  logitUpdate<- glm(ObsData$Y ~ -1 +offset(qlogis(EY.AW)) + H.AW, family='binomial')
  summary(logitUpdate)
  
  # epsilon denote the resulting maximum likelihood estimate of the coefficient 
  # on the clever covariate H.AW
  epsilon <- logitUpdate$coef
  epsilon
  
  # Then update the initial expected outcome :
  # logit(E(Y|A,W)_updated) = logit(E(Y|A,W))+ epsilon*H(A|W)
  EY.AW.updated<- plogis(qlogis(EY.AW)+ epsilon*H.AW)
  EY.1W.updated<- plogis(qlogis(EY.1W)+ epsilon*H.1W)
  EY.0W.updated<- plogis(qlogis(EY.0W)+ epsilon*H.0W)
  
  # Estimate the statistical parameter by substituting the targeted predictions into
  # G-Computation's formula
  
  TMLE_Case4[r]<-mean(EY.1W.updated-EY.0W.updated)
  
  
  
  ##------------------------------------------------------------------------------------------------------##
  ##--------------------------------Case5 : Treatment model correctly specified---------------------------##
  ##------------------------------------------------------------------------------------------------------##
  gm <- glm(Y ~ A + w1 + w2 + w3, family = binomial, data = ObsData)
  # Prediction for A, A = 1 and, A = 0
  QAW_0 <- predict(gm, type = "response")
  Q1W_0 = predict(gm, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")
  Q0W_0 = predict(gm, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")
  # Plug in G-model
  Gcomputation_Case5[r]<-mean(Q1W_0 - Q0W_0)
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  # Step 2 estimation and prediction of the propensity score (ps)
  psm <- glm(A ~ w1 + w2 + w3 +w4+w3*w4,  family = binomial, data = ObsData)
  P.1W = predict(psm, type = "response")
  #For each individual, assign the predicted probability of being not a smoker to P.0W
  P.0W<- 1- P.1W
  #calculate the Clever Covariate for each individual
  H.AW<- as.numeric(ObsData$A==1)/P.1W - as.numeric(ObsData$A==0)/P.0W
  #Clever Covariate for those who are smoker - the interpretation of the inverse
  #of the probability of being a smoker (H.1W= H given Smoker==1 and covariates W)
  H.1W<- 1/P.1W
  #And for those who are not
  H.0W<- -1/P.0W
  #This party is basically calculating the IPTW_ATE for later comparison
  IPTW_Case5[r] <-mean(as.numeric(ObsData$A==1)/P.1W*ObsData$Y) - mean(as.numeric(ObsData$A==0)/P.0W*ObsData$Y)
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  Q1W_AIPW <- mean((as.numeric(ObsData$A==1)) * (ObsData$Y - Q1W_0) / P.1W + Q1W_0)
  Q0W_AIPW <- mean((as.numeric(ObsData$A==0)) * (ObsData$Y - Q0W_0) / P.0W + Q0W_0)
  AIPW_Case5[r]<-Q1W_AIPW-Q0W_AIPW
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################  
  logitUpdate<- glm(ObsData$Y ~ -1 +offset(qlogis(QAW_0)) + H.AW, family='binomial')
  summary(logitUpdate)
  epsilon <- logitUpdate$coef
  epsilon
  # Then update the initial expected outcome :
  # logit(E(Y|A,W)_updated) = logit(E(Y|A,W))+ epsilon*H(A|W)
  QAW_0.updated<- plogis(qlogis(QAW_0)+ epsilon*H.AW)
  Q1W_0.updated<- plogis(qlogis(Q1W_0)+ epsilon*H.1W)
  Q0W_0.updated<- plogis(qlogis(Q0W_0)+ epsilon*H.0W)
  TMLE_Case5[r]<- mean(Q1W_0.updated - Q0W_0.updated)
  
  
  ##----------------------------Case4: Using Super Learner-----------------------##
  ##-----------------------------------Step 1:-----------------------------------##
  ##-------------------------Estimate the initial E0(Y |A, W)--------------------##
  ##-----------------------------------------------------------------------------##
  #Define the predictor for the Superlearner (Exposure + Baseline Covariates)
  X<- subset(ObsData, select=c(A, w1, w2, w3, w4))
  
  #Specify the Superlearner Library
  set.seed(2310)
  SL.library <- c("SL.glm", "SL.step","SL.glm.interaction","SL.gam","SL.randomForest","SL.rpart")
  #Estimate the conditional expectation of outcome, given exposure(A) and covariates(w)
  EY.AW_SL<- SuperLearner(Y=ObsData$Y, X=X, SL.library=SL.library, family="binomial")
  EY.AW_SL
  
  #Obtain the initial estimate of the expected outcome, given exposure(Smoker), and covariates(W)
  EY.AW <- predict(EY.AW_SL, newdata =ObsData)$pred
  EY.AW
  
  #Similarly, obtain the expected outcome for those who under exposure
  EY.1W<- predict(EY.AW_SL, newdata=data.frame(A = 1, ObsData [,c("w1","w2","w3","w4")]), type = "response")$pred
  
  #And those who not
  EY.0W<- predict(EY.AW_SL, newdata=data.frame(A = 0, ObsData [,c("w1","w2","w3","w4")]), type = "response")$pred
  
  
  #Calculating the ATE following G-Computation method (simple subtitution)
  #as the difference in potential outcomes Y1, Y0, corresponding to A=1,0 respectively
  Gcomputation_Case6[r]<- mean(EY.1W- EY.0W)
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################
  
  #Using Superlearner to predict the propensity score
  #which is the probability of the exposure - mother being smoker denoted P0(A=| 1 X)
  #given a set of covariates (W- our predictors)
  set.seed(2310)
  P.AW_SL<- SuperLearner(Y=ObsData$A, X=subset(X, select=c(w1, w2, w3, w4)),
                         SL.library=SL.library, family="binomial")
  
  #For each individual, assign the predicted probability of being a smoker to P.1W
  P.1W<- P.AW_SL$SL.predict
  
  #For each individual, assign the predicted probability of being not a smoker to P.0W
  P.0W<- 1- P.1W
  #calculate the Clever Covariate for each individual
  H.AW<- as.numeric(ObsData$A==1)/P.1W - as.numeric(ObsData$A==0)/P.0W
  
  #Clever Covariate for those who are smoker - the interpretation of the inverse
  #of the probability of being a smoker (H.1W= H given Smoker==1 and covariates W)
  H.1W<- 1/P.1W
  
  #And for those who are not
  H.0W<- -1/P.0W
  
  #This party is basically calculating the IPTW_ATE for later comparison
  IPTW_Case6[r] <-mean(as.numeric(ObsData$A==1)/P.1W*ObsData$Y) - mean(as.numeric(ObsData$A==0)/P.0W*ObsData$Y)
  
  ##########################################################################################################
  ##########################################################################################################
  ########################################################################################################## 
  
  Q1W_AIPW <- mean((as.numeric(ObsData$A==1)) * (ObsData$Y - EY.1W) / P.1W + EY.1W)
  Q0W_AIPW <- mean((as.numeric(ObsData$A==0)) * (ObsData$Y - EY.0W) / P.0W + EY.0W)
  AIPW_Case6[r]<-Q1W_AIPW-Q0W_AIPW
  
  ##########################################################################################################
  ##########################################################################################################
  ##########################################################################################################
  
  logitUpdate<- glm(ObsData$Y ~ -1 +offset(qlogis(EY.AW)) + H.AW, family='binomial')
  summary(logitUpdate)
  
  # epsilon denote the resulting maximum likelihood estimate of the coefficient 
  # on the clever covariate H.AW
  epsilon <- logitUpdate$coef
  epsilon
  
  # Then update the initial expected outcome :
  # logit(E(Y|A,W)_updated) = logit(E(Y|A,W))+ epsilon*H(A|W)
  EY.AW.updated<- plogis(qlogis(EY.AW)+ epsilon*H.AW)
  EY.1W.updated<- plogis(qlogis(EY.1W)+ epsilon*H.1W)
  EY.0W.updated<- plogis(qlogis(EY.0W)+ epsilon*H.0W)
  
  # Estimate the statistical parameter by substituting the targeted predictions into
  # G-Computation's formula
  
  TMLE_Case6[r]<-mean(EY.1W.updated-EY.0W.updated)
  
  
}

######################################################################################
################################### Inperpret ########################################
######################################################################################
### Calculate mean estimates
result_Gcomputation_Case1 <- mean(Gcomputation_Case1)
result_Gcomputation_Case2 <- mean(Gcomputation_Case2)
result_Gcomputation_Case3 <- mean(Gcomputation_Case3)
result_Gcomputation_Case4 <- mean(Gcomputation_Case4)
result_IPTW_Case1 <- mean(IPTW_Case1)
result_IPTW_Case2 <- mean(IPTW_Case2)
result_IPTW_Case3 <- mean(IPTW_Case3)
result_IPTW_Case4 <- mean(IPTW_Case4)
result_AIPW_Case1 <- mean(AIPW_Case1)
result_AIPW_Case2 <- mean(AIPW_Case2)
result_AIPW_Case3 <- mean(AIPW_Case3)
result_AIPW_Case4 <- mean(AIPW_Case4)
result_TMLE_Case1 <- mean(TMLE_Case1)
result_TMLE_Case2 <- mean(TMLE_Case2)
result_TMLE_Case3 <- mean(TMLE_Case3)
result_TMLE_Case4 <- mean(TMLE_Case4)
result_TMLE_Case5<-mean(TMLE_Case5)
results <- rbind(result_Gcomputation_Case1, result_Gcomputation_Case2, result_Gcomputation_Case3, result_Gcomputation_Case4,
                 result_IPTW_Case1, result_IPTW_Case2, result_IPTW_Case3, result_IPTW_Case4, 
                 result_AIPW_Case1, result_AIPW_Case2, result_AIPW_Case3, result_AIPW_Case4,
                 result_TMLE_Case1, result_TMLE_Case2, result_TMLE_Case3, result_TMLE_Case4,result_TMLE_Case5)
evaluate.Performance<-function(estimation, true_ate){
  ave<-mean(estimation)
  biaspercent<-abs(mean(estimation-true_ate)*100)
  bias <-mean(estimation-true_ate)
  var<-var(estimation)
  mse<-mean((estimation-true_ate)^2)
  z <- bias/sqrt(var)
  cover<-(pnorm(1.96-z)-pnorm(-1.96-z))*100
  data.frame(ave,biaspercent,bias,mse,cover)
}
df4<-rbind(evaluate.Performance(Gcomputation_Case1, True_ATE),
           evaluate.Performance(IPTW_Case1, True_ATE),
           evaluate.Performance(AIPW_Case1, True_ATE),
           evaluate.Performance(TMLE_Case1, True_ATE),
           evaluate.Performance(Gcomputation_Case2, True_ATE),
           evaluate.Performance(IPTW_Case2, True_ATE),
           evaluate.Performance(AIPW_Case2, True_ATE),
           evaluate.Performance(TMLE_Case2, True_ATE),
           evaluate.Performance(Gcomputation_Case5, True_ATE),
           evaluate.Performance(IPTW_Case5, True_ATE),
           evaluate.Performance(AIPW_Case5, True_ATE),
           evaluate.Performance(TMLE_Case5, True_ATE),
           evaluate.Performance(Gcomputation_Case3, True_ATE),
           evaluate.Performance(IPTW_Case3, True_ATE),
           evaluate.Performance(AIPW_Case3, True_ATE),
           evaluate.Performance(TMLE_Case3, True_ATE),
           evaluate.Performance(Gcomputation_Case4, True_ATE),
           evaluate.Performance(IPTW_Case4, True_ATE),
           evaluate.Performance(AIPW_Case4, True_ATE),
           evaluate.Performance(TMLE_Case4, True_ATE),
           evaluate.Performance(Gcomputation_Case6, True_ATE),
           evaluate.Performance(IPTW_Case6, True_ATE),
           evaluate.Performance(AIPW_Case6, True_ATE),
           evaluate.Performance(TMLE_Case6, True_ATE))
quantile(Gcomputation_Case1, c(.025, .975))
quantile(IPTW_Case1, c(.025, .975))
quantile(AIPW_Case1, c(.025, .975))
quantile(TMLE_Case1, c(.025, .975))
quantile(Gcomputation_Case2, c(.025, .975))
quantile(IPTW_Case2, c(.025, .975))
quantile(AIPW_Case2, c(.025, .975))
quantile(TMLE_Case2, c(.025, .975))
quantile(Gcomputation_Case5, c(.025, .975))
quantile(IPTW_Case5, c(.025, .975))
quantile(AIPW_Case5, c(.025, .975))
quantile(TMLE_Case5, c(.025, .975))
quantile(Gcomputation_Case3, c(.025, .975))
quantile(IPTW_Case3, c(.025, .975))
quantile(AIPW_Case3, c(.025, .975))
quantile(TMLE_Case3, c(.025, .975))
quantile(Gcomputation_Case4, c(.025, .975))
quantile(IPTW_Case4, c(.025, .975))
quantile(AIPW_Case4, c(.025, .975))
quantile(TMLE_Case4, c(.025, .975))


