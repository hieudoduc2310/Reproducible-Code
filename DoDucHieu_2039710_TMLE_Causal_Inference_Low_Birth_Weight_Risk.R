### Title:    TMLE_ATE Thesis project
### Author:   DO Duc Hieu - 2039710
### Created:  2020-09-03
### Modified: 2020-09-25

#Load required library
library(qdapRegex)
library(varhandle)
library(TMLE_ATE)
library(SuperLearner)
set.seed(2310)


#Load the data
data<- read.csv('Birthweigth.csv')
str(data)
head(data)
#Function to extract value from factor and return as numeric (BWT, Age, LWT)
extract_text <- function(text) {
  temp <- qdapRegex::ex_between(text, "[", ",")[[1]]
  return(as.numeric(temp))
}
data[,"BWT"]<-sapply(data[,"BWT"], extract_text)
data[,"Age"]<-sapply(data[,"Age"], extract_text)
data[,"LWT"]<-sapply(data[,"LWT"], extract_text)
#Function to convert factor "False", "True" to 0,1
convert_logical <- function(x) {
  temp <- unfactor(x)
  temp[temp=="True"]= 1
  temp[temp=="False"]= 0
  x<- as.integer(unlist(temp))
  return(x)
}
data$Low <- convert_logical(data$Low)
data$Smoker <- convert_logical(data$Smoker)
data$Hypertension <- convert_logical(data$Hypertension)
data$UI <- convert_logical(data$UI)
#recheck the variable type
str(data)
#Normalize the data
data$BWT<- (data$BWT-min(data$BWT))/(max(data$BWT)-min(data$BWT))
#recheck the data
#str(data)
head(data)

##Checking for Multicolinearity in the data:
m_col <- cor(data)
corrplot(m_col)
data$MotherID<-NULL
data$BWT<-NULL
#Clone the dataset to the future extended dataset
ex_data<-data
##-----------------------------------------------------------------------------##
##------------------------------------EDA--------------------------------------##
##-----------------------------------------------------------------------------##
tabletemp = table(data.frame(data$Smoker,data$Low))
tabletemp
100*tabletemp/rowSums(tabletemp)
chisq.test(tabletemp, correct=FALSE)
#The test statistic is 4.924, and the p-value is .026. Since .026 < .05 we reject H0.
#The exploratory analysis revealed that low birth weight babies are more prevalent among mothers who smoke during pregnancy than among mothers who don't (40.54 vs. 25.22%).
#The formal statistical test produced a p-value of .026 indicating that the data provide enough evidence to conclude that low birth weight is significantly related to smoking.
#The risks associated with smoking, specifically smoking during pregnancy, should be communicated to pregnant women in as many ways as possible (including pamphlets, Internet, prenatal caregivers, etc.).


##-----------------------------------------------------------------------------##
##------------------------------Implement the algorithm------------------------##
##-----------------------------------------------------------------------------##

##-----------------------------------Step 1:-----------------------------------##
##-------------------------Estimate the initial E0(Y |A, W)--------------------##
##-----------------------------------------------------------------------------##
##With Y is the outcome, A is exposure and W is set of covariates
#Prepare the data for the algorithm
#Define the predictor for the Superlearner (Exposure + Baseline Covariates)
X<- subset(data, select=c(Smoker, LWT, Race, Age, PTL, Hypertension, UI, FTV))

#Cloning X1,X0 from the total population
X1<-X0<-X

#DF X1 set individuals who under the exposure (those who smoke)
X1$Smoker<-1

#DF X0 set individuals who not under the exposure (those who not smoke)
X0$Smoker<-0

#Specify the Superlearner Library
SL.library <- c("SL.glm")
SL.library <- c("SL.glm","SL.step","SL.gam","SL.glm.interaction","SL.ranger")
SL.library <- c("SL.glm")


#Estimate the conditional expectation of outcome, given exposure(A) and covariates(w)
EY.AW_SL<- SuperLearner(Y=data$Low, X=X, SL.library=SL.library, family="binomial")
EY.AW_SL

#Obtain the initial estimate of the expected outcome, given exposure(Smoker), and covariates(W)
EY.AW <- predict(EY.AW_SL, newdata=data)$pred
EY.AW
ex_data$EY.AW <- EY.AW

#Similarly, obtain the expected outcome for those who under exposure
EY.1W<- predict(EY.AW_SL, newdata=X1)$pred
ex_data$EY.1W <- EY.1W

#And those who not
EY.0W<- predict(EY.AW_SL, newdata=X0)$pred
ex_data$EY.0W <- EY.0W

#Calculating the ATE following G-Computation method (simple subtitution)
#as the difference in potential outcomes Y1, Y0, corresponding to A=1,0 respectively
Gcomputation_ATE<- mean(EY.1W- EY.0W)
Gcomputation_ATE


# *NOTE: In this step we estimate the inial conditional expectation of outcome under
# exposure and covariate. The different between TMLE_ATE and G-computation method (or Simple estimator
# subtition) is that: G-com will imediately calculate the ATE base on these potential
# outcomes (EY.1W and EY.0W)

##----------------------------------Step 2:---------------------------------##
##---------------Estimate the  P0(A = 1|W) (or Propensity score)------------##
##--------------------------------------------------------------------------##
1/0.-5

#Using Superlearner to predict the propensity score
#which is the probability of the exposure - mother being smoker denoted P0(A=| 1 X)
#given a set of covariates (W- our predictors)
P.AW_SL<- SuperLearner(Y=data$Smoker, X=subset(X, select=c(Age, LWT, Race, PTL, Hypertension, UI, FTV)),
                       SL.library=SL.library, family="binomial")

#For each individual, assign the predicted probability of being a smoker to P.1W
P.1W<- P.AW_SL$SL.predict
ex_data$P.1W<-P.1W
summary(P.1W)
#For each individual, assign the predicted probability of being not a smoker to P.0W
P.0W<- 1- P.1W
ex_data$P.0W<-P.0W
summary(P.1W)
##-----------------------------------Step 4.1:------------------------------##
##-------------------Update initial estimate of E0(Y |A, W)-----------------##
##--------------------------------------------------------------------------##


#calculate the Clever Covariate for each individual
H.AW<- as.numeric(data$Smoker==1)/P.1W - as.numeric(data$Smoker==0)/P.0W
ex_data$H.AW<-H.AW

#Clever Covariate for those who are smoker - the interpretation of the inverse
#of the probability of being a smoker (H.1W= H given Smoker==1 and covariates W)
H.1W<- 1/P.1W
ex_data$H.1W<-H.1W

#And for those who are not
H.0W<- -1/P.0W
ex_data$H.0W<-H.0W

#This party is basically calculating the IPTW_ATE for later comparison
IPTW_ATE <-mean(as.numeric(data$Smoker==1)/P.1W*data$Low) - mean(as.numeric(data$Smoker==0)/P.0W*data$Low)
IPTW_ATE

##----------------------------------Step 4.2:------------------------------##
##----------------------------Target the estimator-------------------------##
##-------------------------------------------------------------------------##

# Note*: we generate updated (???targeted???) estimates of the
# set of potential outcomes, incorporating information from
# the propensity score (denominator part of the clever covariate) in order 
# to reduce the bias.
# Conceptually, this step regress residuals from initial outcome E(Y|A,W) onto 
# a covariates that is the function of the propensity score H(A,W):
# Y= E(Y|A,W)+ epsilon*H(A,W)
# Run a logistic regression of the outcome Y on the clever covariate H(A, W) with the logit (offset) of the
# initial estimator. Ignore the intercept with -1. Our aim is to estimate delta -
# the resulting coefficient infront of the clever covariate by fitting the 
# following logistic regression model.

logitUpdate<- glm(data$Low ~ -1 +offset(qlogis(EY.AW)) + H.AW, family='binomial')
summary(logitUpdate)
logitUpdate<- glm(data$Low ~ -1 +offset(qlogis(EY.AW)) + H.1W, family='binomial')
logitUpdate$coef
logitUpdate<- glm(data$Low ~ -1 +offset(qlogis(EY.AW)) + H.0W, family='binomial')
logitUpdate$coef
# delta denote the resulting maximum likelihood estimate of the coefficient 
# on the clever covariate H.AW
delta <- logitUpdate$coef
delta
ex_data$epsilon<-delta
# Then update the initial expected outcome :
# logit(E(Y|A,W)_updated) = logit(E(Y|A,W))+ epsilon*H(A|W)
EY.AW.updated<- plogis(qlogis(EY.AW)+ delta*H.AW)
ex_data$EY.AW.utd<- EY.AW.updated
EY.1W.updated<- plogis(qlogis(EY.1W)+ delta*H.1W)
ex_data$EY.1W.updated<-EY.1W.updated
EY.0W.updated<- plogis(qlogis(EY.0W)+ delta*H.0W)
ex_data$EY.0W.updated<-EY.0W.updated

# Estimate the statistical parameter by substituting the targeted predictions into
# G-Computation's formula
TMLE_ATE<- mean(EY.1W.updated - EY.0W.updated)

#Simple Comparision
c(Gcomputation_ATE, IPTW_ATE, TMLE_ATE)

#Under the causal
# under our assumptions, the risk of having a low-weigth baby is essentially 11% higher if the 
# mother was a smoker when having the baby than who is not.
# Checking the extended data_set
colnames(ex_data)

##----------------------------------Step 5:------------------------------##
##------------------Inference on the simulation's result-----------------##
##-----------------------------------------------------------------------##

#calculating the influence curve (book p.96):
#Basically, the fomular of influence curve can be represent as:
#IC = The clever covariate * the residual between observed outcome and targeted prediction
#plus the diference between targeted prediction given A=1 and A=0
#minus the target parameter

#Put what we have estimated into the formula
IC <- H.AW*(data$Low - EY.AW.updated) + EY.1W.updated - EY.0W.updated - TMLE_ATE

#Now, we can calculate the variance of TMLE with the sample variance of the estimated
#influence curve, scaled by the sample size

#estimate sigma^2 with the variance of the IC divided by n
varHat.IC<-var(IC)/nrow(data)

varHat.IC
#standard error estimate
se <- sqrt(varHat.IC)
se
#obtain 95% two-sided confidence intervals:
c(TMLE_ATE+qnorm(0.05/2, lower.tail=T)*se,
  TMLE_ATE+qnorm(0.05/2, lower.tail=F)*se)

#calculate the pvalue
2* pnorm( abs(TMLE_ATE/se), lower.tail=F )

