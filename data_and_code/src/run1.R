
### run1.R
### 
### Meridith Bartley, Ephraim Hanks, Erin Schliep, Tyler Wagner, Nathan Wikle
###
### This script performs the analysis for the first holdout set, (random25)

# source in functions
source("./src/functions.R")

# load data
lake.data.180627 <- readRDS("./data/lake_data_180627.rds")

################################################################################
###### 1. Create appropriate data subsets, based on holdouts
################################################################################

# store subsetted data in 
sub.data <- subsetData(lake.data.180627)[[1]] 

################################################################################
###### 2. Create initial input structures for MCMC function
################################################################################

## list to store important constants (like X, n, K, etc.)
data <- list()
data$Y <- as.matrix(sub.data$Ytrain)
#number of lakes
data$n <- nrow(data$Y) 
#number of variables
data$K <- ncol(data$Y)  
# covariate matrix (remove lat and lon)
data$X <- as.matrix(cbind(rep(1,data$n), sub.data$Xtrain[,-c(1,2)])) 
# number of covariates (including intercept)
data$p <- ncol(data$X) 
# compute which nutrients are observed for each lake
data$obs.mat <- matrix(0, nrow=data$n, ncol=data$K) 
for(k in 1:data$K){
  data$obs.mat[which(is.na(data$Y[,k]) != T), k] <- 1
}

# put the lat and lons here:
data$LLxtest <- sub.data$Xtest[,1:2] 
data$LLxtrain <- sub.data$Xtrain[,1:2]

# create X and Y for test sets
data$Ytest <- as.matrix(sub.data$Ytest)
data$Xtest <- as.matrix(cbind(rep(1,nrow(sub.data$Xtest)), 
                              sub.data$Xtest[,-c(1,2)]))

## save hyperparameters for beta and Sigma priors to a list
priors <- list()
priors$Sigma.beta <- diag(1000, data$K * data$p)
priors$Tau <- diag(1, data$K)
priors$nu <- data$K + 1

# initialize beta values to linear regression parameter estimates
Beta <- matrix(0,nrow=data$K,ncol=data$p)
for(k in 1:data$K){
  Beta[k,] <- as.vector(lm(data$Y[,k] ~ -1 + data$X)$coef)
}

## save initial values for beta, Sigma, and missing Zs
pars <- list()
pars$beta <- as.vector(t(Beta))
pars$Sigma <- diag(.5, data$K)
pars$Z <- matrix(0,ncol=data$K, nrow=data$n)
pars$Z <- data$Y
pars$Z <- Z.update(data, pars, priors)$Z

################################################################################
###### 3. Run MCMC on appoproriate subset
################################################################################

# run mcmc for 2000 iterations
run1 <- drive.lakes(data, pars, priors, 2000, 1000) 
saveRDS(run1, file = "./data/run1_180627.rds")				

################################################################################
###### 4. Generate samples from the posterior predictive distribution
################################################################################

# predict all responses (including secchi)
pred1 <- predict.all(run1, 200, 2000, 1800)	
saveRDS(pred1, file = "./data/pred1_180627.rds")			

# predict nutrient responses, conditional on secchi
pred1.cond <- predict.cond(run1, 200, 2000, 1800)
saveRDS(pred1.cond, file = "./data/pred1_cond_180627.rds")











