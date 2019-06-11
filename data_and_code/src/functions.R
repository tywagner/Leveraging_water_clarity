
### functions.R
### 
### Meridith Bartley, Ephraim Hanks, Erin Schliep, Tyler Wagner, Nathan Wikle
###
### Functions used in the data analysis from the manuscript,
###   "Increasing accuracy of lake nutrient predictions in thousands of 
###     lakes by leveraging water clarity data," 2019.

### Load necessary functions:
library(mvtnorm)
library(MCMCpack)

################################################################################
###### 1. Data Cleanup
################################################################################

subsetData <- function(data.list){
  # Subset data into lists of training and holdout sets, based on previously
  #   chosen holdout sets.
  # Input:
  #   data.list: a large list of training and holdout sets
  # Output:
  #   A list with X.test, Y.test, X.train, and Y.train structures for 
  # the various holdout sets.

  # grab holdout 
  holdout <- data.list$Holdout
  
  ### random25_holdout (hold out TRUE values)
  random25 <- list()
  random25$Ytest <- data.list$logY[holdout[, 'random25_holdout'], ]
  random25$Ytrain <- data.list$logY[!holdout[, 'random25_holdout'], ]
  random25$Xtest <- data.list$X[holdout[, 'random25_holdout'], ]
  random25$Xtrain <- data.list$X[!holdout[, 'random25_holdout'], ]
  
  ### random75_holdout (hold out TRUE values)
  random75 <- list()
  random75$Ytest <- data.list$logY[holdout[, 'random75_holdout'], ]
  random75$Ytrain <- data.list$logY[!holdout[, 'random75_holdout'], ]
  random75$Xtest <- data.list$X[holdout[, 'random75_holdout'], ]
  random75$Xtrain <- data.list$X[!holdout[, 'random75_holdout'], ]
  
  ### hu4_ag50_holdout (hold out TRUE values)
  hu4_ag50 <- list()
  hu4_ag50$Ytest <- data.list$logY[holdout[, 'hu4_ag50_holdout'], ]
  hu4_ag50$Ytrain <- data.list$logY[!holdout[, 'hu4_ag50_holdout'], ]
  hu4_ag50$Xtest <- data.list$X[holdout[, 'hu4_ag50_holdout'], ]
  hu4_ag50$Xtrain <- data.list$X[!holdout[, 'hu4_ag50_holdout'], ]
  
  ### hu4_random50_holdout (hold out TRUE values)
  hu4_random50 <- list()
  hu4_random50$Ytest <- data.list$logY[holdout[,'hu4_random50_holdout'], ]
  hu4_random50$Ytrain <- data.list$logY[!holdout[,'hu4_random50_holdout'], ]
  hu4_random50$Xtest <- data.list$X[holdout[,'hu4_random50_holdout'], ]
  hu4_random50$Xtrain <- data.list$X[!holdout[,'hu4_random50_holdout'], ]
  
  ### hu4_strat75_holdout (hold out TRUE values)
  hu4_strat75 <- list()
  hu4_strat75$Ytest <- data.list$logY[holdout[, 'hu4_strat75_holdout'], ]
  hu4_strat75$Ytrain <- data.list$logY[!holdout[, 'hu4_strat75_holdout'], ]
  hu4_strat75$Xtest <- data.list$X[holdout[, 'hu4_strat75_holdout'], ]
  hu4_strat75$Xtrain <- data.list$X[!holdout[, 'hu4_strat75_holdout'], ]
  
  ### cluster_strat75_holdout (hold out TRUE values)
  cluster_strat75 <- list()
  cluster_strat75$Ytest <- data.list$logY[holdout[, 'cluster_strat75_holdout'], ]
  cluster_strat75$Ytrain <- data.list$logY[!holdout[, 'cluster_strat75_holdout'], ]
  cluster_strat75$Xtest <- data.list$X[holdout[, 'cluster_strat75_holdout'], ]
  cluster_strat75$Xtrain <- data.list$X[!holdout[, 'cluster_strat75_holdout'], ]
  
  ### cluster_random50_holdout (hold out TRUE values)
  cluster_random50 <- list()
  cluster_random50$Ytest <- data.list$logY[holdout[, 'cluster_random50_holdout'], ]
  cluster_random50$Ytrain <- data.list$logY[!holdout[, 'cluster_random50_holdout'], ]
  cluster_random50$Xtest <- data.list$X[holdout[, 'cluster_random50_holdout'], ]
  cluster_random50$Xtrain <- data.list$X[!holdout[, 'cluster_random50_holdout'], ]
  
  final.list <- list(random25 = random25, 
                     random75 = random75, 
                     hu4_ag50 = hu4_ag50, 
                     hu4_random50 = hu4_random50, 
                     hu4_strat75 = hu4_strat75,
                     cluster_strat75 = cluster_strat75,
                     cluster_random50 = cluster_random50)
  
  final.list
}

################################################################################
###### 2. MCMC Code
################################################################################

beta.update <- function(data,pars,priors){
  # Fully conditional conjugate updates of beta (mean parameters).
  # Input:
  #   data: list containing data variables (X, K, p)
  #   pars: list with parameter values (beta, Sigma, Z)
  #   priors: list with prior hyperparameters (Sigma.beta, 
  # Output:
  #   Updated beta parameters (within pars list).
  
  A <- rep(0,data$p*data$K)	
  B <- solve(priors$Sigma.beta)
  Sigma.Inv <- solve(pars$Sigma)
  for(i in 1:data$n){
    B <- B + kronecker(diag(1, data$K), data$X[i,]) %*% Sigma.Inv %*% 
      t(kronecker(diag(1, data$K), data$X[i,]))
    A <- A + kronecker(diag(1, data$K), data$X[i,]) %*% Sigma.Inv %*% pars$Z[i,]
  }
  pars$beta <- as.vector(rmvnorm(1, solve(B) %*% A, solve(B)))
  return(pars)	
}

Sigma.update <- function(data,pars,priors){
  # Fully conditional conjugate update of Sigma (KxK covariance matrix).
  # Input:
  #   data: list containing data variables (X, K, p)
  #   pars: list with parameter values (beta, Sigma, Z)
  #   priors: list with prior hyperparameters (Sigma.beta, 
  # Output:
  #   Updated Sigma matrix (within pars list).
  
  TmXB <- t(pars$Z - matrix(kronecker(diag(data$K), data$X) %*% pars$beta, 
                            ncol = data$K))
  pars$Sigma <- riwish(data$n + priors$nu, TmXB %*% t(TmXB) + priors$Tau)	
  return(pars)	
}

Z.update <- function(data,pars,priors){
  # Impute missing Z values using conditional normal updates.
  # Input:
  #   data: list containing data variables (X, K, p)
  #   pars: list with parameter values (beta, Sigma, Z)
  #   priors: list with prior hyperparameters  
  # Output:
  #   Updated missing Z values (within pars list).
  
  C <- pars$Sigma
  
  for(i in which(apply(data$obs.mat,1,sum) != data$K)){
    Mu <- as.vector(t(kronecker(diag(1, data$K), data$X[i,])) %*% pars$beta)
    f <- which(is.na(data$Y[i,]) == T)
    pars$Z[i,f] <- c(rmvnorm(1,
      mean = Mu[f] + C[f,-f] %*% solve(C[-f,-f]) %*% (pars$Z[i,-f] - Mu[-f]),
      sigma = C[f,f] - C[f,-f] %*% solve(C[-f,-f]) %*% C[-f,f]))
  }
  
  return(pars)
}

drive.lakes <- function(data,pars,priors,iters,print.out){
  # MCMC (Gibbs updates) from posterior of [beta, Sigma, Z.missing | .]
  # Input:
  #   data: list containing data variables (X, K, p, n)
  #   pars: list with initial parameter values (beta, Sigma, Z)
  #   priors: list with priro hyperparameters
  #   iters: number of iterations (i.e., number of simulations)
  #   print.out: print out iteration number (for convenience)
  # Output:
  #   List with posterior samples for beta, Sigma, Z.missing.
  
  # set up structure to hold output
  out=list()
  out$data <- data
  out$priors <- priors
  out$beta <- array(dim=c(data$K,data$p,iters+1))
  out$beta[,,1] <- matrix(pars$beta,ncol=data$p,byrow=T)
  out$Sigma <- array(dim=c(data$K,data$K,iters+1))
  out$Sigma[,,1] <- pars$Sigma
  out$Z <- array(dim=c(data$n,data$K,iters+1))
  out$Z[,,1] <- pars$Z
  
  # Run MCMC for iters
  for(j in 1:iters){
    
    # update Z
    pars <- Z.update(data,pars,priors)
    # update Sigma
    pars <- Sigma.update(data,pars,priors)
    # update beta
    pars <- beta.update(data,pars,priors)
    
    # print out iteration
    if(j%%print.out==0){
      print(j)
    }
    
    # save MCMC updates
    out$beta[,,j+1] <- matrix(pars$beta,ncol=data$p,byrow=T)
    out$Sigma[,,j+1] <- pars$Sigma
    out$Z[,,j+1] <- pars$Z   
  }  
  
  # return output
  return(out)
}

################################################################################
###### 3. Functions used to generate posterior predictive samples
################################################################################

predict.all <- function(run,b,e,iters){
  # Posterior predictive distribution of entire multivariate response.
  # Input:
  #   run: list containing hold-out set data (Xtest) and responses (Y)
  #   b: burn-in indice (avoid sampling from here)
  #   e: end point of indices (should be size of MCMC samples)
  #   iters: number of post. pred. samples
  # Output:
  #   Return a given number of simulations from the posterior 
  # predictive distribution.
  
  # holdout covariates
  newX <- as.matrix(run$data$Xtest)
  # structure for post. predicted values of Y
  Y.pred <- array(dim=c(nrow(newX), ncol(run$data$Y), iters))
  # choose b posterior samples out of the e total, to use for post. pred. samples
  s <- sample(b:e,iters)
  # generate posterior predictive samples
  for(j in 1:iters){
    for(i in 1:nrow(newX)){
      Y.pred[i,,j] <- as.vector(rmvnorm(1,
        as.vector(newX[i,] %*% t(run$beta[,,s[j]])), run$Sigma[,,s[j]]))	
    }
  }
  # save and resurn output
  out <- list()
  out$Y.pred <- Y.pred
  out$Y.true <- run$data$Ytest
  return(out)
}

predict.cond <- function(run,b,e,iters){
  # Posterior predictive distribution of conditional response (i.e., 
  #   predict all covariates except Secchi).
  # Input:
  #   run: list containing hold-out set data (Xtest) and responses (Y)
  #   b: burn-in indice (avoid sampling from here)
  #   e: end point of indices (should be size of MCMC samples)
  #   iters: number of post. pred. samples
  # Output:
  #   Return a given number of simulations from the posterior 
  # predictive distribution of logTP, logTN, logChla cond. on Secchi.
  
  # holdout covariates
  newX <- as.matrix(run$data$Xtest) 
  # structure for post. predicted values of Y
  Y.cond <- array(dim=c(nrow(newX), ncol(run$data$Y), iters))
  # choose b posterior samples out of the e total, to use for post. pred. samples
  s <- sample(b:e,iters)
  
  # post. pred samples of response conditional on secchi depth
  for(i in 1:nrow(newX)){
    if(is.na(run$data$Ytest[i,4]) == T){
      Y.cond[i,,]=NA	
    }
    if(is.na(run$data$Ytest[i,4]) != T){
      for(j in 1:iters){
        f <- 4 # this assumes secchi is the 4th variable
        Y.cond[i,f,j] <- run$data$Ytest[i,f]
        C <- run$Sigma[,,s[j]]
        Mu <- newX[i,] %*% t(run$beta[,,s[j]])
        Y.cond[i,-f,j] <- rmvnorm(1,
          mean = Mu[-f] + C[-f,f] %*% solve(C[f,f]) %*% 
            (run$data$Ytest[i,f] - Mu[f]),
          sigma = C[-f,-f] - C[-f,f] %*% solve(C[f,f]) %*% C[f,-f])	
      }	
    }
  }	
  
  # save output and return
  out <- list()
  out$Y.cond <- Y.cond
  out$Y.true <- run$data$Ytest
  return(out)
}

################################################################################
###### 4. Calculate prediction error
################################################################################

percent.error <- function(pred.out){
  # Calculate median percent error (MPE) for a given prediction output.
  # Input:
  #   pred.out: list with prediction output (either full or conditional)
  # Output:
  #   Matrix of median percent error (MPE) values.
  
  point.est <- apply(exp(pred.out[[1]]), c(1,2), median)
  true <- exp(predOut[[2]])
  true[which(true==0)] <- NA
  #print(head(true))
  p.error <- apply(abs(point.est/true-1), 2, median, na.rm=T)
  return(p.error)	
}	

rmse <- function(pred.out){
  # Calculate the root mean square error (RMSE) for a given prediction output.
  # Input:
  #   pred.out: list with prediction output (either full or conditional)
  # Output:
  #   Matrix of RMSE values.
  
  error <- sqrt(apply((pred.out$Y.true - 
                         apply(pred.out$Y.pred, c(1,2), mean))^2,
                      2, mean, na.rm=T))
  
  error
}
